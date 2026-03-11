use fp2::traits::Fp as FpTrait;
use rand_core::{CryptoRng, RngCore};

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

// TODO: what's the optimum here?
const KARATSUBA_THRESHOLD: usize = 4;

/// Trait for arithmetic for univariate polynomials in Fp[X]
pub trait Poly<Fp: FpTrait>:
    Clone
    + Default
    + Display
    + Index<usize>
    + IndexMut<usize>
    + Sized
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + MulAssign<Self>
    + MulAssign<Fp>
{
    fn new_from_ele(a: &Fp) -> Self;
    fn new_from_slice(a: &[Fp]) -> Self;
    fn set_from_slice(&mut self, a: &[Fp]);

    fn degree(&self) -> Option<usize>;

    fn reverse(&self) -> Self;

    fn scale(&self, a: &Fp) -> Self;

    fn evaluate(&self, a: &Fp) -> Fp;

    fn root_from_quadratic_leaves(leaves: &[[Fp; 3]]) -> Self;

    fn resultant_from_roots_with_tree(&self, eval_tree: &EvalTree<Fp>) -> Fp;
    fn resultant_from_roots_horner(&self, roots: &[Fp]) -> Fp;
}

/// Precomputed data for repeated multipoint evaluation at the same set of roots.
///
/// The product tree is stored in a flat layout: `flat_tree[lyr]` is a single
/// contiguous `Vec<Fp>` of length exactly `n_roots`, holding the non-leading
/// coefficients of every node at that layer.  Because every layer is the product
/// tree of a degree-`n` monic polynomial, the sum of node degrees at each layer
/// is always exactly `n` — so `n` slots suffice regardless of the layer.
///
/// `offsets[lyr]` has `n_nodes_at_lyr + 1` entries.  Node `j` at layer `lyr`
/// occupies `flat_tree[lyr][offsets[lyr][j] .. offsets[lyr][j+1]]`.
/// Its degree is `offsets[lyr][j+1] - offsets[lyr][j]`.
pub struct EvalTree<Fp: fp2::traits::Fp> {
    /// flat_tree[lyr][offsets[lyr][j]..offsets[lyr][j+1]] = non-leading coefficients
    /// of node j at layer lyr.  Length of flat_tree[lyr] is always n_roots.
    pub flat_tree: Vec<Vec<Fp>>,
    pub rev_q_inv: Vec<Fp>,
    pub n_roots: usize,
    /// offsets[lyr][j] = start of node j's coefficient window in flat_tree[lyr].
    pub offsets: Vec<Vec<usize>>,
}

impl<Fp: fp2::traits::Fp> EvalTree<Fp> {
    pub fn new(roots: &[Fp]) -> Self {
        use crate::polynomial_ring::poly::Polynomial;

        let n_roots = roots.len();
        if n_roots == 0 {
            return Self {
                flat_tree: vec![],
                rev_q_inv: vec![],
                n_roots: 0,
                offsets: vec![],
            };
        }

        let (flat_tree, offsets) = Polynomial::<Fp>::product_tree_flat(roots);
        let n_layers = flat_tree.len();
        let n = n_roots;

        // rev(Q_root)^{-1} mod x^n
        //
        // Q_root has degree n with leading coeff 1.  Its non-leading coefficients
        // are stored in flat_tree[n_layers-1][0..n].
        // rev(Q)[i] = Q[n-i], so:
        //   rev_q[0]   = Q[n]   = 1          (the monic leading term, reversed)
        //   rev_q[i]   = Q[n-i] = flat_tree[top][n-i]   for i in 1..=n
        let top = &flat_tree[n_layers - 1];
        let mut rev_q = vec![Fp::ZERO; n + 1];
        rev_q[0] = Fp::ONE;
        for i in 1..=n {
            rev_q[i] = top[n - i];
        }
        let mut rev_q_inv = vec![Fp::ZERO; n];
        Polynomial::<Fp>::inv_mod_xn(&mut rev_q_inv, &rev_q, n);

        Self {
            flat_tree,
            rev_q_inv,
            n_roots,
            offsets,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Polynomial<Fp: FpTrait> {
    coeffs: Vec<Fp>,
}

impl<Fp: FpTrait> Polynomial<Fp> {
    /// Create a polynomial from a finite field element.
    pub fn new_from_ele(a: &Fp) -> Self {
        Self { coeffs: vec![*a] }
    }

    /// Create a polynomial from a slice of finite field elements.
    pub fn new_from_slice(a: &[Fp]) -> Self {
        Self { coeffs: a.to_vec() }
    }

    /// Set the coefficients of a polynomial from a slice
    pub fn set_from_slice(&mut self, a: &[Fp]) {
        self.coeffs.copy_from_slice(a)
    }

    /// The length of the polynomial. TODO: should we trim trailing zeros? If so, how often?
    fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Return the degree of the polynomial.
    // TODO: should we trim trailing zeros and mutate while computing the degree or do the following
    // and allow trailing zeros to remain.
    pub fn degree(&self) -> Option<usize> {
        let mut i = self.len() - 1;
        while i > 0 && self.coeffs[i].is_zero() == u32::MAX {
            i -= 1;
        }
        if i == 0 && self.coeffs[0].is_zero() == u32::MAX {
            None
        } else {
            Some(i)
        }
    }

    /// Return the constant coefficient of the polynomial
    pub fn constant_coefficient(&self) -> Option<Fp> {
        if self.len() == 0 {
            return None;
        }
        Some(self.coeffs[0])
    }

    /// Reverse the coefficients of self in place.
    fn reverse_into(&mut self) {
        self.coeffs.reverse();
    }

    /// Return the polynomial with coefficents reversed.
    pub fn reverse(&self) -> Self {
        let mut r = self.clone();
        r.reverse_into();
        r
    }

    /// Return 0xFFFFFFFF if self and other represent the same polynomial.
    /// Otherwise, return 0x00000000.
    pub fn equals(&self, other: &Self) -> u32 {
        // TODO: do I want this constant time?
        // eg: let mut equals = ct_u32_eq(self.len() as u32, other.len() as u32);
        if self.len() != other.len() {
            return 0;
        }

        let mut equals = u32::MAX;
        for i in 0..self.len() {
            equals &= self.coeffs[i].equals(&other[i]);
        }
        equals
    }

    /// Return 0xFFFFFFFF if self is zero, otherwise, return 0x00000000.
    pub fn is_zero(&self) -> u32 {
        let mut is_zero = u32::MAX;
        for i in 0..self.len() {
            is_zero &= self.coeffs[i].is_zero();
        }
        is_zero
    }

    /// Compute f <-- f + g, assumes that the length of f and g are the same.
    fn add_into(f: &mut [Fp], g: &[Fp]) {
        debug_assert!(f.len() >= g.len());
        for i in 0..g.len() {
            f[i] += g[i];
        }
    }

    /// Compute f <-- f - g, assumes that the length of f and g are the same.
    fn sub_into(f: &mut [Fp], g: &[Fp]) {
        debug_assert!(f.len() >= g.len());
        for i in 0..g.len() {
            f[i] -= g[i];
        }
    }

    /// Compute f * g with O(len(f) * len(g)) Fp multiplications using
    /// schoolbook multiplication. Assumes that fg has enough space for
    /// the result (len(f) + len(g) - 1).
    fn schoolbook_multiplication(fg: &mut [Fp], f: &[Fp], g: &[Fp]) {
        debug_assert!(fg.len() >= f.len() + g.len() - 1);

        for i in 0..f.len() {
            for j in 0..g.len() {
                if i == 0 || j + 1 == g.len() {
                    fg[i + j] = f[i] * g[j]
                } else {
                    fg[i + j] += f[i] * g[j]
                }
            }
        }
    }

    /// Compute f * g with ~O(n^1.5) Fp multiplications using Karastuba multiplication.
    /// Assumes that fg has enough space for the result (len(f) + len(g) - 1).
    fn karatsuba_multiplication(fg: &mut [Fp], f: &[Fp], g: &[Fp]) {
        debug_assert!(fg.len() >= f.len() + g.len() - 1);

        // Ensure that the degree of f is larger or equal to g (for balancing the split later)
        if f.len() < g.len() {
            Self::karatsuba_multiplication(fg, g, f);
            return;
        }

        // If g has length zero, then we set f * g to be zero.
        if g.is_empty() {
            for c in fg.iter_mut() {
                *c = Fp::ZERO;
            }
            return;
        }

        // If g has length one we simply scale all coefficients by g0.
        if g.len() == 1 {
            let g0 = g[0];
            for i in 0..f.len() {
                fg[i] = f[i] * g0;
            }
            return;
        }

        // When f is linear, then we also know g is linear and we can
        // save one multiplication compared to the naive method.
        if f.len() == 2 {
            let t1 = f[0] + f[1];
            let t2 = g[0] + g[1];
            fg[0] = f[0] * g[0];
            fg[2] = f[1] * g[1];
            fg[1] = t1 * t2;
            fg[1] -= fg[0];
            fg[1] -= fg[2];
            return;
        }

        // For small degree f, g we use basic multiplication strategies with
        // O(n^2) operations. TODO: set the right bound for when to fall back
        // to this.
        if f.len() <= KARATSUBA_THRESHOLD {
            Self::schoolbook_multiplication(fg, f, g);
            return;
        }

        // We are now at the point where f.len() >= 4 and we will split f into
        // a high and low part at floor(f.len() / 2) to perform a divide and
        // conquer strategy by Karastuba.
        let nf = f.len() / 2;
        let mf = f.len() - nf;

        // When g is particularly small we cannot split g into two halves, so we
        // have to modify Karatsuba to split f at floor(f.len() / 2) and then
        // multiply f_lo and f_hi by the whole of g.
        if g.len() <= nf {
            // We can compute f_lo * g directly into the first nf + g.len() - 1 elements
            // of the output.
            Self::karatsuba_multiplication(&mut fg[..nf + g.len() - 1], &f[..nf], g);

            // We then compute f_hi * g which will have length mf + g.len() - 1.
            // The bottom g.len() - 1 elements we need to add into fg while we can
            // directly copy the rest of the elements into the top of fg.
            let mut fg_hi = vec![Fp::ZERO; mf + g.len() - 1];
            Self::karatsuba_multiplication(&mut fg_hi, &f[nf..nf + mf], g);
            Self::add_into(&mut fg[nf..nf + g.len() - 1], &fg_hi[..g.len() - 1]);
            fg[nf + g.len() - 1..].copy_from_slice(&fg_hi[g.len() - 1..]);

            return;
        }

        // We are now at the point where we can split f = f_lo + x^nf * f_hi and
        // g = g_lo + x^nf * g_hi without an issue, with deg(f) >= deg(g).
        let mg = g.len() - nf;

        // The idea is that we now perform three calls to karatsuba_multiplication
        // on half-length inputs. Writing f * g = fg_lo + x^fn (fg_mid) + x^2*fn (fg_hi)
        // we need to compute:
        //
        // - fg_lo = f_lo * g_lo
        // - fg_mid = f_lo * g_hi + f_hi * g_lo
        //          = (f_lo + f_hi) * (g_lo + g_hi) - f_lo * g_lo - f_hi * g_hi
        // - fg_hi  = f_hi * g_hi
        //
        // Which means we need to compute only three multiplications: f_lo * g_lo,
        // f_hi * g_hi and (f_lo + f_hi) * (f_lo + f_hi).
        //
        // - fg_lo will have length 2*fn - 1 and fill the bottom of fg[..2*fn - 1].
        // - fg_hi will have length fm + gm - 1 and fill the top (without overlap)
        //   of fg[2*fn..]
        // - fg_mid will have length max(fn, fm) + max(fn, gm) - 1 and will fill
        //   from fg[fn..fn + fg_mid_len], to do this we will add the result into
        //   fg after the copies above.

        // The fg_lo and fg_hi parts of the computation are simple and we can
        // multiply straight into the output buffer.
        Self::karatsuba_multiplication(&mut fg[..nf + nf - 1], &f[..nf], &g[..nf]);
        Self::karatsuba_multiplication(&mut fg[nf + nf..], &f[nf..], &g[nf..]);

        // Now we compute f_lo + f_hi and g_lo + g_hi.

        // As nf is floor(len(f) / 2) then mf will either be nf or nf + 1, so we
        // can fit the sum into mf space and then add nf elements to it.
        // TODO: work with less allocations?
        let mut f_mid = f[nf..].to_vec();
        Self::add_into(&mut f_mid[..nf], &f[..nf]);

        // For g_lo + g_hi we need to be more careful, as len(g) <= len(f) we might
        // have that mg < nf.
        let mut g_mid = vec![Fp::ZERO; nf.max(mg)];
        if mg < nf {
            g_mid.copy_from_slice(&g[..nf]);
            Self::add_into(&mut g_mid[..mg], &g[nf..]);
        } else {
            g_mid.copy_from_slice(&g[nf..]);
            Self::add_into(&mut g_mid[..nf], &g[..nf]);
        }

        // Now we have both pieces, we compute their product into a temp buffer.
        let mut fg_mid = vec![Fp::ZERO; mf + nf.max(mg) - 1];
        Self::karatsuba_multiplication(&mut fg_mid, &f_mid, &g_mid);

        // We then compute (f_lo + f_hi) * (f_lo + f_hi) - f_lo * g_lo - f_hi * g_hi
        // with two subtractions.
        Self::sub_into(&mut fg_mid[..nf + nf - 1], &fg[..nf + nf - 1]);
        Self::sub_into(&mut fg_mid[..mf + mg - 1], &fg[nf + nf..]);

        // We now need to set the remaining window of fg which doesn't overlap with
        // the high and low pieces, and add in the last part. The dumb thing to do
        // is to zero out the single element which has yet to be copied over, but
        // we could save one addition by instead doing one copy at fg[nf + nf - 1]
        // and then add in two windows from fg[nf .. nf + nf - 1] and then fg[nf + nf ...].
        fg[nf + nf - 1] = Fp::ZERO;
        Self::add_into(&mut fg[nf..nf + mf + nf.max(mg) - 1], &fg_mid);
    }

    /// TODO: implement other polynomial multiplication (Karatsuba, probably)
    fn mul_into(fg: &mut [Fp], f: &[Fp], g: &[Fp]) {
        debug_assert!(f.len() + g.len() - 1 <= fg.len());
        Self::karatsuba_multiplication(fg, f, g)
    }

    /// TODO: testing only.
    pub fn basic_mul(&self, other: &Self) -> Self {
        let mut coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];
        Self::schoolbook_multiplication(&mut coeffs, &self.coeffs, &other.coeffs);
        Self { coeffs }
    }

    /// TODO: testing only.
    pub fn mul_low_test(h: &mut [Fp], f: &[Fp], g: &[Fp], n: usize) {
        Self::mul_low(h, f, g, n);
    }

    /// TODO: testing only.
    pub fn mul_middle_test(&self, other: &Self) -> Self {
        let mut coeffs = vec![Fp::ZERO; other.len()];
        Self::mul_middle(&mut coeffs, &self.coeffs, &other.coeffs);
        Self { coeffs }
    }

    /// TODO: testing only.
    pub fn karatsuba_mul(&self, other: &Self) -> Self {
        let mut coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];
        Self::karatsuba_multiplication(&mut coeffs, &self.coeffs, &other.coeffs);
        Self { coeffs }
    }

    /// Set self to it's negative.
    fn set_neg(&mut self) {
        for x in self.coeffs.iter_mut() {
            x.set_neg();
        }
    }

    /// Set self <- self + other
    fn set_add(&mut self, other: &Self) {
        if self.len() < other.len() {
            self.coeffs.resize(other.len(), Fp::ZERO);
        }
        Self::add_into(&mut self.coeffs, &other.coeffs);
    }

    /// Set self <- self - other
    fn set_sub(&mut self, other: &Self) {
        if self.len() < other.len() {
            self.coeffs.resize(other.len(), Fp::ZERO);
        }
        Self::sub_into(&mut self.coeffs, &other.coeffs);
    }

    /// Set self <- self * other
    fn set_mul(&mut self, other: &Self) {
        let mut fg_coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];
        Self::mul_into(&mut fg_coeffs, &self.coeffs, &other.coeffs);
        self.coeffs = fg_coeffs;
    }

    /// Multiply all coefficients of the polynomial by a element of the finite field.
    pub fn scale_into(&mut self, c: &Fp) {
        for x in self.coeffs.iter_mut() {
            *x *= *c;
        }
    }

    /// Return  c * self for some c in the finite field.
    pub fn scale(&self, c: &Fp) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.scale_into(c);
        r
    }

    /// Multiply all coefficients of the polynomial by a small value.
    pub fn scale_small_into(&mut self, k: i32) {
        for x in self.coeffs.iter_mut() {
            x.set_mul_small(k);
        }
    }

    /// Return c * self for some small c
    pub fn scale_small(&self, k: i32) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.scale_small_into(k);
        r
    }

    /// Computation of a product tree given an array of coefficients of
    /// quadratic polynomials as leaves. Avoids performing many memory
    /// allocations, but the code is a little harder to read.
    pub fn product_tree_quadratic_leaves(leaves: &[[Fp; 3]]) -> Vec<Vec<Fp>> {
        // Number of leaves and depth for the tree.
        let n = leaves.len();
        let log_n = usize::BITS - (2 * n - 1).leading_zeros();

        // Each layer of the tree is a vector of Fp elements, the leaves are n
        // polynomials of degree 2, then the next are n/2 polynomials of degree
        // 4 (and potentially one of degree 2 if n is odd) and so on until the
        // root of the tree is a single polynomial of degree 2*n.
        let mut layers: Vec<Vec<Fp>> = Vec::with_capacity(log_n as usize);

        // Store the quadratic polynomials inside an input buffer
        let mut buf_in: Vec<Fp> = leaves
            .iter()
            .flat_map(|poly| poly.iter().copied())
            .collect();
        let mut buf_out = vec![Fp::ZERO; 3 * n];

        // Store the leaves of the tree in layers.
        layers.push(buf_in.clone());

        // Iterate over the remaining layers,
        for i in 1..log_n {
            // At each layer in the tree the degree of each polynomial is at most 2^i.
            let deg = 1 << i;
            let len = deg + 1;
            let out_len = 2 * deg + 1;

            // The majority of the multiplications will be h <- f * g where f and g
            // both have degree `deg`. This will be done k times when (2 * n = 2 * deg * k + r)
            let k = n / deg;
            let r = (n << 1) & ((1 << (i + 1)) - 1);

            // Keep track of the slices for the length deg + 1 polys for input and
            // length 2*degree + 1 polynomials as output for multiplication.
            let mut idx_in = 0;
            let mut idx_out = 0;

            // Compute k full multiplications into the buffer.
            for _ in 0..k {
                Self::mul_into(
                    &mut buf_out[idx_out..idx_out + out_len],
                    &buf_in[idx_in..idx_in + len],
                    &buf_in[idx_in + len..idx_in + 2 * len],
                );
                idx_out += out_len;
                idx_in += 2 * len;
            }

            if r > 0 {
                // By this point, we have consumed pairs of polynomials of degree `deg`
                // If `r` is larger than `deg` then we need to do an unbalanced multiplication
                // of a degree `deg` polynomial with a degree `r - deg` polynomial
                if r > deg {
                    let len_rem = r - deg + 1;
                    Self::mul_into(
                        &mut buf_out[idx_out..idx_out + r + 1],
                        &buf_in[idx_in..idx_in + len],
                        &buf_in[idx_in + len..idx_in + len + len_rem],
                    );
                } else {
                    // Otherwise we want to copy the polynomial from the in buffer to the
                    // out buffer
                    buf_out[idx_out..idx_out + r + 1]
                        .copy_from_slice(&buf_in[idx_in..idx_in + r + 1]);
                }
            }

            // Move the multiplication result into the input buffer and append
            // this into the layers.
            std::mem::swap(&mut buf_in, &mut buf_out);
            layers.push(buf_in.clone());
        }
        debug_assert!(layers.len() == log_n as usize);

        layers
    }

    /// Computes the polynomial product(leaves) assuming a tree given by
    /// product_tree_quadratic_leaves.
    pub fn root_from_quadratic_leaf_tree(tree: &[Vec<Fp>]) -> Self {
        let n = tree[0].len() / 3;
        let root = tree.last().unwrap();
        Self {
            coeffs: root[..2 * n + 1].to_vec(),
        }
    }

    /// Compute the product of many quadratic leaves, we do this without storing the
    /// tree as only the root is required from sqrt velu.
    pub fn root_from_quadratic_leaves(leaves: &[[Fp; 3]]) -> Self {
        let n = leaves.len();
        let log_n = usize::BITS - (2 * n - 1).leading_zeros();

        let mut buf_in: Vec<Fp> = leaves.iter().flat_map(|p| p.iter().copied()).collect();
        let mut buf_out = vec![Fp::ZERO; 3 * n];

        for i in 1..log_n {
            let deg = 1 << i;
            let len = deg + 1;
            let out_len = 2 * deg + 1;
            let k = n / deg;
            let r = (n << 1) & ((1 << (i + 1)) - 1);

            let mut idx_in = 0;
            let mut idx_out = 0;
            for _ in 0..k {
                Self::mul_into(
                    &mut buf_out[idx_out..idx_out + out_len],
                    &buf_in[idx_in..idx_in + len],
                    &buf_in[idx_in + len..idx_in + 2 * len],
                );
                idx_out += out_len;
                idx_in += 2 * len;
            }
            if r > 0 {
                if r > deg {
                    let len_rem = r - deg + 1;
                    Self::mul_into(
                        &mut buf_out[idx_out..idx_out + r + 1],
                        &buf_in[idx_in..idx_in + len],
                        &buf_in[idx_in + len..idx_in + len + len_rem],
                    );
                } else {
                    buf_out[idx_out..idx_out + r + 1]
                        .copy_from_slice(&buf_in[idx_in..idx_in + r + 1]);
                }
            }
            std::mem::swap(&mut buf_in, &mut buf_out);
        }

        Self {
            coeffs: buf_in[..2 * n + 1].to_vec(),
        }
    }

    /// Evaluate a polynomial at a value `a` using Horner's method.
    pub fn evaluate(&self, a: &Fp) -> Fp {
        // Handle degree 0 and 1 cases early.
        if self.len() == 0 {
            return Fp::ZERO;
        } else if self.len() == 1 {
            return self.coeffs[0];
        }

        // Iterate from the last to the first coefficient using Horner's method
        // to evaluate the polynomial.
        let mut bi = Fp::ZERO;
        let deg = self.degree().unwrap();
        for i in 0..=deg {
            bi = *a * bi + self.coeffs[deg - i];
        }
        bi
    }

    /// Set self to a random value
    pub fn set_rand<R: CryptoRng + RngCore>(&mut self, rng: &mut R) {
        for x in self.coeffs.iter_mut() {
            x.set_rand(rng);
        }
    }

    /// Return a new random polynomial with length d
    pub fn rand<R: CryptoRng + RngCore>(rng: &mut R, d: usize) -> Self {
        let mut r = Self {
            coeffs: vec![Fp::ZERO; d],
        };
        r.set_rand(rng);
        r
    }

    /// Compute only the low n coefficients of the product f * g.
    ///
    /// Uses a Karatsuba-style halving recursion with split point m = ceil(n/2):
    ///
    ///   f = f_lo + x^m * f_hi,   g = g_lo + x^m * g_hi
    ///   (f*g)[0..n] = (f_lo*g_lo)[0..2m-1]          full sub-product, kept up to x^{n-1}
    ///               + x^m * (f_lo*g_hi)[0..rem]      cross terms
    ///               + x^m * (f_hi*g_lo)[0..rem]
    ///
    /// The f_hi*g_hi term begins at x^{2m} ≥ x^n and never contributes.
    ///
    /// `h` must have length exactly `n`.
    fn mul_low(h: &mut [Fp], f: &[Fp], g: &[Fp], n: usize) {
        debug_assert_eq!(h.len(), n);

        if n == 0 {
            return;
        }

        // Coefficients at index ≥ n in f or g cannot affect the low n coefficients.
        let f = &f[..f.len().min(n)];
        let g = &g[..g.len().min(n)];

        if f.is_empty() || g.is_empty() {
            for c in h.iter_mut() {
                *c = Fp::ZERO;
            }
            return;
        }

        // Schoolbook for small n — no allocation overhead.
        if n <= KARATSUBA_THRESHOLD || f.len() <= 2 || g.len() <= 2 {
            for c in h.iter_mut() {
                *c = Fp::ZERO;
            }
            for i in 0..f.len() {
                for j in 0..g.len().min(n - i) {
                    h[i + j] += f[i] * g[j];
                }
            }
            return;
        }

        let m = (n + 1) / 2; // ceil(n/2) — the split point
        let rem = n - m;

        let (f_lo, f_hi) = f.split_at(f.len().min(m));
        let (g_lo, g_hi) = g.split_at(g.len().min(m));

        // Full product f_lo * g_lo (length ≤ 2m-1).
        // Low m coefficients → h[0..m].  Overflow (if any) → h[m..n].
        let ll_len = f_lo.len() + g_lo.len() - 1;
        let mut ll = vec![Fp::ZERO; ll_len];
        Self::karatsuba_multiplication(&mut ll, f_lo, g_lo);

        let copy_lo = ll_len.min(m);
        h[..copy_lo].copy_from_slice(&ll[..copy_lo]);
        for c in h[copy_lo..m].iter_mut() {
            *c = Fp::ZERO;
        }

        if rem == 0 {
            return;
        }

        // Seed h[m..n] with overflow from f_lo*g_lo, then add cross terms.
        for c in h[m..].iter_mut() {
            *c = Fp::ZERO;
        }
        if ll_len > m {
            let overflow = (ll_len - m).min(rem);
            h[m..m + overflow].copy_from_slice(&ll[m..m + overflow]);
        }

        // Cross term 1: f_lo * g_hi, low `rem` coefficients, added into h[m..n].
        if !g_hi.is_empty() {
            let mut cross = vec![Fp::ZERO; rem];
            Self::mul_low(&mut cross, f_lo, g_hi, rem);
            Self::add_into(&mut h[m..], &cross);
        }

        // Cross term 2: f_hi * g_lo, low `rem` coefficients, added into h[m..n].
        if !f_hi.is_empty() {
            let mut cross = vec![Fp::ZERO; rem];
            Self::mul_low(&mut cross, f_hi, g_lo, rem);
            Self::add_into(&mut h[m..], &cross);
        }
    }

    /// Compute the middle n coefficients of f * g.
    ///
    /// The result satisfies `h[i] = (f * g)[n-1+i]` for `i` in `0..n`.
    ///
    /// Requires `f.len() == 2*n - 1` and `g.len() == n`.
    ///
    /// Uses the Hanrot-Quercia-Zimmermann (HQZ) 3-recursive-call algorithm.
    ///
    /// With `half = n/2` (floor) and `half_up = n - half` (ceil):
    ///
    ///   tmplo = f  (length 2n-1)
    ///   tmplo[0..n-half_up] += f[half_up..]       (self-overlapping sum)
    ///
    ///   a = mul_middle(tmplo[0..2*half_up-1],           g[half..],  half_up)
    ///   c = mul_middle(tmplo[half_up..half_up+2*half-1], g[..half], half)
    ///   b = mul_middle(f[half_up..3*half_up-1],          g_diff,    half_up)
    ///   where g_diff = g[half..],  g_diff[half_up-half..] -= g[..half]
    ///
    ///   h[0..half_up]  = a - b
    ///   h[half_up..n]  = c + b[0..half]
    fn mul_middle(h: &mut [Fp], f: &[Fp], g: &[Fp]) {
        let n = g.len();
        debug_assert_eq!(h.len(), n);
        debug_assert_eq!(f.len(), 2 * n - 1);

        if n == 0 {
            return;
        }
        if n == 1 {
            h[0] = f[0] * g[0];
            return;
        }
        if n == 2 {
            // Explicit 2-coefficient case.
            // (f*g)[1] = f[0]*g[1] + f[1]*g[0]
            // (f*g)[2] = f[1]*g[1] + f[2]*g[0]
            h[0] = f[0] * g[1] + f[1] * g[0];
            h[1] = f[1] * g[1] + f[2] * g[0];
            return;
        }

        let half = n / 2; // floor(n/2)
        let half_up = n - half; // ceil(n/2)

        // Build tmplo: f with f shifted right by half_up added into its low portion.
        // tmplo[i] = f[i] + f[i+half_up]  for i in 0..(f.len()-half_up)
        // tmplo[i] = f[i]                  for i in (f.len()-half_up)..f.len()
        let mut tmplo = f.to_vec(); // length 2n-1
        let add_len = f.len() - half_up; // = (2n-1) - half_up
        Self::add_into(&mut tmplo[..add_len], &f[half_up..]);

        // a: middle half_up coefficients using the first window of tmplo.
        // Input length: 2*half_up-1  ✓ (≤ 2n-1 since half_up ≤ n)
        let mut a = vec![Fp::ZERO; half_up];
        Self::mul_middle(&mut a, &tmplo[..2 * half_up - 1], &g[half..]);

        // c: middle half coefficients using the second window of tmplo.
        // Window length: 2*half-1  (we use only the part of tmplo we need).
        let mut c = vec![Fp::ZERO; half];
        if half > 0 {
            Self::mul_middle(&mut c, &tmplo[half_up..half_up + 2 * half - 1], &g[..half]);
        }

        // g_diff = g[half..] with g[..half] subtracted from the aligned portion.
        // The offset (half_up - half = n mod 2) handles the odd/even cases:
        //   n even: subtract g[0..half] from g_diff[0..half]  (full overlap)
        //   n odd:  subtract g[0..half] from g_diff[1..half+1] (shifted by 1)
        let mut g_diff = g[half..].to_vec(); // length half_up
        Self::sub_into(&mut g_diff[half_up - half..], &g[..half]);

        // b: middle half_up coefficients using f[half_up..3*half_up-1].
        // Input length: 3*half_up-1 - half_up = 2*half_up-1  ✓
        let mut b = vec![Fp::ZERO; half_up];
        Self::mul_middle(&mut b, &f[half_up..3 * half_up - 1], &g_diff);

        // Combine results.
        for i in 0..half_up {
            h[i] = a[i] - b[i];
        }
        for i in 0..half {
            h[half_up + i] = c[i] + b[i];
        }
    }

    /// Compute the middle `out_n` coefficients of `f * g`, where
    /// `h[i] = (f * g)[in_n - 1 + i]` for `i in 0..out_n`.
    ///
    /// Requires `f.len() == in_n + out_n - 1` and `g.len() == in_n`.
    ///
    /// Unlike the previous version, this function is **allocation-free** in the
    /// recursive case: it recurses directly into the two non-overlapping halves
    /// `h[..half_out]` and `h[half_out..]` of the output slice, eliminating the
    /// two temporary `Vec`s (`h_lo`, `h_hi`) that existed before.
    ///
    /// The only allocation that remains is in the base case (`out_n <= 4`), where
    /// a full product is computed into `scratch` (which is passed in from the
    /// caller and reused across calls).
    fn mul_middle_rect(h: &mut [Fp], f: &[Fp], g: &[Fp], scratch: &mut Vec<Fp>) {
        let out_n = h.len();
        let in_n = g.len();
        debug_assert_eq!(f.len(), in_n + out_n - 1);

        // Balanced base case — delegate to the existing mul_middle.
        if in_n == out_n {
            Self::mul_middle(h, f, g);
            return;
        }

        // Base case: compute full product into scratch, then extract the window.
        // scratch is reused by the caller, so we resize (cheap if already large enough)
        // rather than allocating a fresh Vec each time.
        if out_n <= 4 {
            let full_len = f.len() + g.len() - 1;
            scratch.resize(full_len, Fp::ZERO);
            // Zero the scratch buffer before use (resize only zeros new elements).
            for c in scratch.iter_mut() {
                *c = Fp::ZERO;
            }
            Self::karatsuba_multiplication(scratch, f, g);
            h.copy_from_slice(&scratch[in_n - 1..in_n - 1 + out_n]);
            return;
        }

        // Recursive case: split the output range in half and recurse directly into
        // the two non-overlapping halves of h — zero allocations.
        //
        // f has length in_n + out_n - 1.  Split at half_out:
        //
        //   h[0..half_out]  = mul_middle_rect(f[..in_n+half_out-1], g, half_out)
        //     check: f_lo.len() = in_n + half_out - 1  ✓  (= in_n + half_out - 1)
        //
        //   h[half_out..out_n] = mul_middle_rect(f[half_out..], g, out_hi)
        //     check: f_hi.len() = in_n+out_n-1-half_out = in_n+out_hi-1  ✓
        //
        // Correctness: h_lo[i] = (f*g)[in_n-1+i] for i in 0..half_out uses only
        //   f[0..in_n+half_out-1], which is f_lo.
        //   h_hi[i] = (f*g)[in_n-1+half_out+i] for i in 0..out_hi uses only
        //   f[half_out..in_n+out_n-1], which is f_hi.
        let half_out = out_n / 2;

        // Recurse into h[..half_out] directly — no allocation.
        Self::mul_middle_rect(&mut h[..half_out], &f[..in_n + half_out - 1], g, scratch);

        // Recurse into h[half_out..] directly — no allocation.
        Self::mul_middle_rect(&mut h[half_out..], &f[half_out..], g, scratch);
    }

    /// Compute `h[0..n]` such that `(f * h) ≡ 1 mod x^n`.
    ///
    /// Requires `f[0]` to be invertible. Directly translates the reference Rust
    /// `_inv_mod_xn` (Hanrot-Quercia-Zimmermann §5.1), which works for any `n`:
    ///
    ///   half     = n / 2  (floor)
    ///   half_up  = n - half  (ceil)
    ///
    ///   1. Recurse: compute h_k[0..half_up] s.t. f * h_k ≡ 1 mod x^{half_up}
    ///   2. Error:   e[0..half_up] = mul_middle(f_shifted[0..2*half_up-1], h_k[0..half_up])
    ///              where f_shifted = f[1..] zero-padded to length 2*half_up-1
    ///   3. Correct: h[half_up..n] = -(h_k[0..half] * e[0..half])[0..half]
    ///
    /// The error computation gives `half_up` values; only the first `half` are used
    /// in the correction. This naturally handles both even and odd `n` with no
    /// special-casing — for odd `n`, one extra error value is computed and discarded.
    ///
    /// `h` must have length exactly `n`.
    pub fn inv_mod_xn(h: &mut [Fp], f: &[Fp], n: usize) {
        debug_assert_eq!(h.len(), n);

        if n == 0 {
            return;
        }

        if f[0].equals(&Fp::ONE) == u32::MAX {
            match n {
                1 => {
                    h[0] = Fp::ONE;
                    return;
                }
                2 => {
                    h[0] = Fp::ONE;
                    h[1] = -f[1];
                    return;
                }
                3 => {
                    h[0] = Fp::ONE;
                    h[1] = -f[1];
                    h[2] = f[1].square() - f[2];
                    return;
                }
                _ => {}
            }
        } else if n == 1 {
            h[0] = Fp::ONE / f[0];
            return;
        }

        let half = n / 2; // floor(n/2)
        let half_up = n - half; // ceil(n/2)

        // Step 1: recursively compute h_k[0..half_up] s.t. f*h_k ≡ 1 mod x^{half_up}.
        Self::inv_mod_xn(&mut h[..half_up], &f[..half_up], half_up);

        // Step 2: compute the error e = (f * h_k)[half_up .. n + half_up - 1]
        //   = middle half_up coefficients of (f_shifted * h_k)
        //   where f_shifted[i] = f[i+1]  (zero-padded to length 2*half_up-1).
        //
        // This is a direct translation of the reference's general-case branch:
        //   tmp_p[..n-1] = f[1..n],  tmp_p[n-1..] = 0
        //   _middlemul(tmplo[0..half_up], tmp_p[0..2*half_up-1], z[0..half_up])
        let f_shift_len = 2 * half_up - 1;
        let mut f_shifted = vec![Fp::ZERO; f_shift_len];
        let copy_end = (1 + f_shift_len).min(f.len());
        if copy_end > 1 {
            f_shifted[..copy_end - 1].copy_from_slice(&f[1..copy_end]);
        }

        let mut e = vec![Fp::ZERO; half_up];
        Self::mul_middle(&mut e, &f_shifted, &h[..half_up]);
        // e[0..half_up] now holds (f*h_k)[half_up..n+half_up-1].
        // We only need e[0..half] for the correction (indices half_up..n).

        // Step 3: correction  h[half_up..n] = -(h_k[0..half] * e[0..half])[0..half]
        let mut corr = vec![Fp::ZERO; half];
        Self::mul_low(&mut corr, &h[..half], &e[..half], half);
        for i in 0..half {
            h[half_up + i] = -corr[i];
        }
    }

    /// Build a product tree over `roots` in a flat layout and return it together
    /// with the offset table.
    ///
    /// # Layout
    ///
    /// Every node in the tree is a *monic* polynomial.  We store only its
    /// **non-leading coefficients**: a degree-`d` node contributes exactly `d`
    /// elements.  Because the sum of node degrees at every layer equals `n`
    /// (the total degree of the root), each layer occupies **exactly `n`
    /// contiguous `Fp` elements**.
    ///
    /// `flat_tree[lyr][offsets[lyr][j] .. offsets[lyr][j+1]]` gives the
    /// non-leading coefficients of node `j` at layer `lyr`.
    ///
    /// # Memory
    ///
    /// * `flat_tree`: `n_layers` allocations of `n` elements each.
    /// * `offsets`: `n_layers` small `usize` vectors.
    /// * One shared `scratch` buffer of ≤ `2n + 2` elements (re-used; no
    ///   per-node allocation after that single up-front allocation).
    pub fn product_tree_flat(roots: &[Fp]) -> (Vec<Vec<Fp>>, Vec<Vec<usize>>) {
        let n = roots.len();
        debug_assert!(n > 0);

        // ------------------------------------------------------------------ //
        // Layer 0 (leaves): n linear monic polys  (x - r_i).                 //
        // Non-leading part of (x - r_i) is just [-r_i], a single element.    //
        // So flat_tree[0] = [-r_0, -r_1, ..., -r_{n-1}], length n.           //
        // offsets[0] = [0, 1, 2, ..., n].                                     //
        // ------------------------------------------------------------------ //
        let leaf_layer: Vec<Fp> = roots.iter().map(|r| -*r).collect();
        let leaf_offsets: Vec<usize> = (0..=n).collect();

        let mut flat_tree: Vec<Vec<Fp>> = vec![leaf_layer];
        let mut offsets: Vec<Vec<usize>> = vec![leaf_offsets];

        // Scratch buffer for monic_mul_nonleading.
        //
        // Layout inside monic_mul_nonleading for inputs of degrees n_f and n_g:
        //   [f_nl | 1]  length n_f+1
        //   [g_nl | 1]  length n_g+1
        //   product     length n_f+n_g+1
        // Total = 2*(n_f+n_g)+3.  The worst case is the root layer where
        // n_f+n_g = n, requiring 2n+3 elements.  Allocate 2n+4 so resize()
        // inside monic_mul_nonleading never needs to grow the buffer.
        let mut scratch: Vec<Fp> = vec![Fp::ZERO; 2 * n + 4];

        // ------------------------------------------------------------------ //
        // Build upward until the single root node remains.                    //
        // ------------------------------------------------------------------ //
        loop {
            let prev_off = offsets.last().unwrap();
            let n_prev_nodes = prev_off.len() - 1;
            if n_prev_nodes == 1 {
                break;
            }

            let prev_flat = flat_tree.last().unwrap();

            // Allocate the next layer's flat buffer (always exactly n elements).
            let mut next_flat = vec![Fp::ZERO; n];
            // offsets for this layer: one entry per node + sentinel.
            let mut next_off: Vec<usize> = Vec::with_capacity((n_prev_nodes + 1) / 2 + 1);
            next_off.push(0);

            let mut out_cursor = 0usize;
            let mut i = 0usize;

            while i < n_prev_nodes {
                if i + 1 < n_prev_nodes {
                    // Multiply the pair of monic children.
                    let l_start = prev_off[i];
                    let l_end = prev_off[i + 1];
                    let r_start = prev_off[i + 1];
                    let r_end = prev_off[i + 2];

                    let n_l = l_end - l_start; // degree of left child
                    let n_r = r_end - r_start; // degree of right child
                    let n_out = n_l + n_r; // degree of parent = sum of child degrees

                    Self::monic_mul_nonleading(
                        &mut next_flat[out_cursor..out_cursor + n_out],
                        &prev_flat[l_start..l_end],
                        &mut scratch,
                        n_l,
                        &prev_flat[r_start..r_end],
                        n_r,
                    );

                    out_cursor += n_out;
                    next_off.push(out_cursor);
                    i += 2;
                } else {
                    // Odd node: carry it unchanged into the next layer.
                    let s = prev_off[i];
                    let e = prev_off[i + 1];
                    let len = e - s;
                    next_flat[out_cursor..out_cursor + len].copy_from_slice(&prev_flat[s..e]);
                    out_cursor += len;
                    next_off.push(out_cursor);
                    i += 1;
                }
            }

            debug_assert_eq!(out_cursor, n);

            flat_tree.push(next_flat);
            offsets.push(next_off);
        }

        (flat_tree, offsets)
    }

    /// Multiply two monic polynomials given only their non-leading coefficients,
    /// writing only the non-leading coefficients of the result.
    ///
    /// Given:
    ///   F(x) = x^n_f + f_{n_f-1} x^{n_f-1} + … + f_0    (stored as `f_nl = [f_0..f_{n_f-1}]`)
    ///   G(x) = x^n_g + g_{n_g-1} x^{n_g-1} + … + g_0    (stored as `g_nl = [g_0..g_{n_g-1}]`)
    ///
    /// Computes H = F * G (monic, degree n_f+n_g) and writes
    ///   `out = [h_0, h_1, ..., h_{n_f+n_g-1}]`   (the n_f+n_g non-leading coefficients).
    ///
    /// Uses `scratch` as a temporary buffer for the full product (length n_f+n_g+1);
    /// the buffer is resized if needed but never shrinks, so repeated calls with
    /// similar sizes have zero allocation cost after the first.
    fn monic_mul_nonleading(
        out: &mut [Fp],
        f_nl: &[Fp],
        scratch: &mut Vec<Fp>,
        n_f: usize,
        g_nl: &[Fp],
        n_g: usize,
    ) {
        debug_assert_eq!(f_nl.len(), n_f);
        debug_assert_eq!(g_nl.len(), n_g);
        debug_assert_eq!(out.len(), n_f + n_g);

        // We need to compute (f_nl | 1) * (g_nl | 1) and keep the first n_f+n_g
        // coefficients (discarding the monic leading term at position n_f+n_g).
        //
        // Rather than building temporary full-poly slices, we expand the product
        // symbolically:
        //
        //   H = F * G
        //     = (x^n_f + Σ f_i x^i) * (x^n_g + Σ g_j x^j)
        //
        // The full product has n_f+n_g+1 coefficients.  We compute it into
        // `scratch` using the Karatsuba routine (which works on arbitrary slices)
        // by temporarily appending the leading 1 to both inputs.
        //
        // Strategy: build the full [f_nl | 1] and [g_nl | 1] in scratch's extra
        // space, multiply, then copy the non-leading n_f+n_g coefficients to `out`.
        //
        // To avoid two extra allocations for [f_nl | 1] and [g_nl | 1] we use a
        // single scratch region laid out as:
        //   scratch[0          .. n_f+1]          = f_nl ++ [1]
        //   scratch[n_f+1      .. n_f+1+n_g+1]   = g_nl ++ [1]
        //   scratch[n_f+n_g+2  .. 2*(n_f+n_g)+3] = product (length n_f+n_g+1)
        //
        // Total scratch needed: 3*(n_f+n_g) + 5, well within our 2n+2 budget.

        let f_start = 0;
        let f_end = n_f + 1;
        let g_start = f_end;
        let g_end = g_start + n_g + 1;
        let p_start = g_end;
        let p_end = p_start + n_f + n_g + 1;

        scratch.resize(p_end, Fp::ZERO);

        // Write [f_nl | 1] into scratch[f_start..f_end].
        scratch[f_start..f_start + n_f].copy_from_slice(f_nl);
        scratch[f_start + n_f] = Fp::ONE;

        // Write [g_nl | 1] into scratch[g_start..g_end].
        scratch[g_start..g_start + n_g].copy_from_slice(g_nl);
        scratch[g_start + n_g] = Fp::ONE;

        // Zero the product region before multiplying (karatsuba assumes it).
        for c in scratch[p_start..p_end].iter_mut() {
            *c = Fp::ZERO;
        }

        // We need to split scratch to call karatsuba with non-overlapping slices.
        // Use split_at_mut to get (inputs_region, product_region).
        let (inputs, product_region) = scratch.split_at_mut(p_start);
        Self::karatsuba_multiplication(
            &mut product_region[..n_f + n_g + 1],
            &inputs[f_start..f_end],
            &inputs[g_start..g_end],
        );

        // Copy the non-leading n_f+n_g coefficients into out (discard the leading 1).
        out.copy_from_slice(&scratch[p_start..p_start + n_f + n_g]);
    }

    /// Build a product tree over `roots` — kept for external callers that need
    /// the full `Vec<Vec<Vec<Fp>>>` layout (e.g. tests).  Internally the
    /// EvalTree now uses `product_tree_flat`.
    pub fn product_tree_from_roots(roots: &[Fp]) -> Vec<Vec<Vec<Fp>>> {
        let leaves: Vec<Vec<Fp>> = roots.iter().map(|r| vec![-*r, Fp::ONE]).collect();
        let mut layers: Vec<Vec<Vec<Fp>>> = vec![leaves];

        // Build upward until only one node remains.
        loop {
            let prev = layers.last().unwrap();
            if prev.len() == 1 {
                break;
            }

            let mut next: Vec<Vec<Fp>> = Vec::with_capacity((prev.len() + 1) / 2);
            let mut i = 0;
            while i < prev.len() {
                if i + 1 < prev.len() {
                    // Multiply the pair.
                    let f = &prev[i];
                    let g = &prev[i + 1];
                    let mut fg = vec![Fp::ZERO; f.len() + g.len() - 1];
                    Self::mul_into(&mut fg, f, g);
                    next.push(fg);
                    i += 2;
                } else {
                    // Odd node out — carry it unchanged.
                    next.push(prev[i].clone());
                    i += 1;
                }
            }
            layers.push(next);
        }

        layers
    }

    pub fn multieval_from_tree(&self, eval_tree: &EvalTree<Fp>) -> Vec<Fp> {
        let n_roots = eval_tree.n_roots;
        if n_roots == 0 {
            return vec![];
        }

        let flat_tree = &eval_tree.flat_tree;
        let offsets = &eval_tree.offsets;
        let n_layers = flat_tree.len();
        let n = n_roots;

        // ------------------------------------------------------------------
        // Build the scaled representation  F = revP * revQ^{-1} mod x^n
        // and store it reversed in `node_buf` (our flat "current" buffer).
        // ------------------------------------------------------------------
        let degp_full = match self.degree() {
            Some(d) => d,
            None => return vec![Fp::ZERO; n_roots],
        };

        // Reduce P mod Q_root if deg(P) >= n so the scaled representation is valid.
        // In VeluSqrt deg(P) = 2*size_J < size_I = n, so this branch is dead in practice.
        let owned: Option<Vec<Fp>> = if degp_full >= n {
            // Q_root non-leading coeffs are in flat_tree[top][0..n]; leading coeff is 1.
            let top = &flat_tree[n_layers - 1];
            let mut rem = self.coeffs.clone();
            rem.resize(degp_full + 1, Fp::ZERO);
            for i in (n..=degp_full).rev() {
                if rem[i].is_zero() != u32::MAX {
                    let c = rem[i];
                    // Subtract c * x^{i-n} * Q_root.  Q_root = x^n + top[n-1]*x^{n-1} + … + top[0].
                    // x^{i-n} * Q_root contributes to positions (i-n)..(i).
                    // Non-leading part: positions (i-n)..(i-1) from top[0..n].
                    // Leading part: position i (coefficient 1).
                    for j in 0..n {
                        rem[i - n + j] -= c * top[j];
                    }
                    // rem[i] -= c * 1 = c, but rem[i] is c so it becomes 0.
                    rem[i] = Fp::ZERO;
                }
            }
            rem.truncate(n);
            Some(rem)
        } else {
            None
        };
        let p_coeffs: &[Fp] = owned.as_deref().unwrap_or(&self.coeffs[..=degp_full]);

        let degp = p_coeffs
            .iter()
            .enumerate()
            .rev()
            .find(|(_, c)| c.is_zero() != u32::MAX)
            .map(|(i, _)| i)
            .unwrap_or(0);

        if p_coeffs[degp].is_zero() == u32::MAX {
            return vec![Fp::ZERO; n_roots];
        }

        // rev_p[i] = p_coeffs[degp - i], zero-padded to length n.
        let mut rev_p = vec![Fp::ZERO; n];
        for i in 0..=degp {
            rev_p[i] = p_coeffs[degp - i];
        }

        // series = rev(P) * rev(Q)^{-1} mod x^n.
        let mut series = vec![Fp::ZERO; n];
        Self::mul_low(&mut series, &rev_p, &eval_tree.rev_q_inv, n);

        // node_buf[i] = series[degp - i]: reverse series back to polynomial order.
        let mut node_buf = vec![Fp::ZERO; n];
        for i in 0..=degp {
            node_buf[i] = series[degp - i];
        }

        // dst_buf: ping-pong target.
        let mut dst_buf = vec![Fp::ZERO; n];

        // Shared scratch buffer for mul_middle_rect — one allocation, reused every call.
        let mut mmr_scratch: Vec<Fp> = Vec::new();

        // ------------------------------------------------------------------
        // Descend the tree.
        //
        // At layer lyr+1 (parent), node j occupies node_buf[off_p[j]..off_p[j+1]].
        // Its two children at layer lyr are at indices child_idx and child_idx+1.
        //
        // For each genuine parent (non-carry):
        //   parent slice  = node_buf[off_p[j] .. off_p[j+1]]    length n_l + n_r
        //   parent[1..]   = parent_tail                          length n_l + n_r - 1
        //   q_l_nl = flat_tree[lyr][off_c[ci]   .. off_c[ci+1]]   length n_l
        //   q_r_nl = flat_tree[lyr][off_c[ci+1] .. off_c[ci+2]]   length n_r
        //
        //   dst[off_c[ci]  ..off_c[ci+1]] = mul_middle_rect(parent_tail, q_r_nl) + parent[..n_l]
        //   dst[off_c[ci+1]..off_c[ci+2]] = mul_middle_rect(parent_tail, q_l_nl) + parent[..n_r]
        //
        // For a carry node: copy the parent slice straight into the child slot.
        // ------------------------------------------------------------------
        for lyr in (0..n_layers - 1).rev() {
            let off_p = &offsets[lyr + 1]; // parent offsets
            let off_c = &offsets[lyr]; // child offsets
            let n_children = off_c.len() - 1;
            let lyr_flat = &flat_tree[lyr];

            let mut child_idx = 0usize;
            for parent_j in 0..off_p.len() - 1 {
                let ps = off_p[parent_j];
                let pe = off_p[parent_j + 1];
                let parent = &node_buf[ps..pe];

                if child_idx + 1 < n_children {
                    let ci = child_idx;
                    let n_l = off_c[ci + 1] - off_c[ci];
                    let n_r = off_c[ci + 2] - off_c[ci + 1];
                    let q_l_nl = &lyr_flat[off_c[ci]..off_c[ci + 1]];
                    let q_r_nl = &lyr_flat[off_c[ci + 1]..off_c[ci + 2]];
                    let parent_tail = &parent[1..]; // length n_l + n_r - 1

                    // Left child result.
                    {
                        let out = &mut dst_buf[off_c[ci]..off_c[ci + 1]];
                        if n_r > 0 {
                            Self::mul_middle_rect(out, parent_tail, q_r_nl, &mut mmr_scratch);
                        } else {
                            out.fill(Fp::ZERO);
                        }
                        Self::add_into(out, &parent[..n_l]);
                    }
                    // Right child result.
                    {
                        let out = &mut dst_buf[off_c[ci + 1]..off_c[ci + 2]];
                        if n_l > 0 {
                            Self::mul_middle_rect(out, parent_tail, q_l_nl, &mut mmr_scratch);
                        } else {
                            out.fill(Fp::ZERO);
                        }
                        Self::add_into(out, &parent[..n_r]);
                    }
                    child_idx += 2;
                } else {
                    // Carry: copy parent slice into the single child slot.
                    let ci = child_idx;
                    dst_buf[off_c[ci]..off_c[ci + 1]].copy_from_slice(&node_buf[ps..pe]);
                    child_idx += 1;
                }
            }

            std::mem::swap(&mut node_buf, &mut dst_buf);
            // No need to zero dst_buf here: every slot will be overwritten before
            // being read at the next layer (mul_middle_rect writes before add_into reads,
            // and carry nodes copy unconditionally).
        }

        // At layer 0 each node has exactly 1 coefficient = P(r_i).
        node_buf[..n_roots].to_vec()
    }

    pub fn multieval(&self, roots: &[Fp]) -> Vec<Fp> {
        let eval_tree = EvalTree::new(roots);
        self.multieval_from_tree(&eval_tree)
    }

    pub fn multieval_horner(&self, roots: &[Fp]) -> Vec<Fp> {
        let mut res = vec![Fp::ZERO; roots.len()];
        for (i, r) in roots.iter().enumerate() {
            res[i] = self.evaluate(&r)
        }
        res
    }

    /// Compute the product self(a_i) for all a_i, given a pre-built eval tree.
    pub fn resultant_from_roots_with_tree(&self, eval_tree: &EvalTree<Fp>) -> Fp {
        let mut res = Fp::ONE;
        for a in self.multieval_from_tree(eval_tree).iter() {
            res *= *a;
        }
        res
    }

    pub fn resultant_from_roots(&self, ai: &[Fp]) -> Fp {
        let eval_tree = EvalTree::new(ai);
        self.resultant_from_roots_with_tree(&eval_tree)
    }

    pub fn resultant_from_roots_horner(&self, ai: &[Fp]) -> Fp {
        let mut res = Fp::ONE;
        for a in ai {
            res *= self.evaluate(a)
        }
        res
    }
}

impl<Fp: FpTrait> Poly<Fp> for Polynomial<Fp> {
    fn new_from_ele(a: &Fp) -> Self {
        Self::new_from_ele(a)
    }
    fn new_from_slice(a: &[Fp]) -> Self {
        Self::new_from_slice(a)
    }
    fn set_from_slice(&mut self, a: &[Fp]) {
        self.set_from_slice(a)
    }

    fn degree(&self) -> Option<usize> {
        self.degree()
    }

    fn reverse(&self) -> Self {
        self.reverse()
    }

    fn scale(&self, a: &Fp) -> Self {
        self.scale(a)
    }

    fn evaluate(&self, a: &Fp) -> Fp {
        self.evaluate(a)
    }

    fn root_from_quadratic_leaves(leaves: &[[Fp; 3]]) -> Self {
        Self::root_from_quadratic_leaves(leaves)
    }

    fn resultant_from_roots_with_tree(&self, eval_tree: &EvalTree<Fp>) -> Fp {
        self.resultant_from_roots_with_tree(eval_tree)
    }

    fn resultant_from_roots_horner(&self, roots: &[Fp]) -> Fp {
        self.resultant_from_roots_horner(roots)
    }
}

impl<Fp: FpTrait> Index<usize> for Polynomial<Fp> {
    type Output = Fp;

    fn index(&self, index: usize) -> &Self::Output {
        self.coeffs.get(index).expect("Index out of bounds")
    }
}

impl<Fp: FpTrait> IndexMut<usize> for Polynomial<Fp> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.coeffs.get_mut(index).expect("Index out of bounds")
    }
}

impl<Fp: FpTrait> Default for Polynomial<Fp> {
    fn default() -> Self {
        Self {
            coeffs: vec![Fp::ONE; 1],
        }
    }
}

impl<Fp: FpTrait> Neg for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn neg(self) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_neg();
        r
    }
}

impl<Fp: FpTrait> Add<Polynomial<Fp>> for Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn add(self, other: Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_add(&other);
        r
    }
}

impl<Fp: FpTrait> Add<&Polynomial<Fp>> for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn add(self, other: &Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_add(other);
        r
    }
}

impl<Fp: FpTrait> AddAssign<Polynomial<Fp>> for Polynomial<Fp> {
    #[inline(always)]
    fn add_assign(&mut self, other: Polynomial<Fp>) {
        self.set_add(&other);
    }
}

impl<Fp: FpTrait> AddAssign<&Polynomial<Fp>> for Polynomial<Fp> {
    #[inline(always)]
    fn add_assign(&mut self, other: &Polynomial<Fp>) {
        self.set_add(other);
    }
}

impl<Fp: FpTrait> Sub<Polynomial<Fp>> for Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn sub(self, other: Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_sub(&other);
        r
    }
}

impl<Fp: FpTrait> Sub<&Polynomial<Fp>> for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn sub(self, other: &Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_sub(other);
        r
    }
}

impl<Fp: FpTrait> SubAssign<&Polynomial<Fp>> for Polynomial<Fp> {
    #[inline(always)]
    fn sub_assign(&mut self, other: &Polynomial<Fp>) {
        self.set_sub(other);
    }
}

impl<Fp: FpTrait> Mul<&Polynomial<Fp>> for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn mul(self, other: &Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_mul(other);
        r
    }
}

// TODO: ideally we don't want to do the moving here?
impl<Fp: FpTrait> MulAssign<Polynomial<Fp>> for Polynomial<Fp> {
    #[inline(always)]
    fn mul_assign(&mut self, other: Polynomial<Fp>) {
        self.set_mul(&other);
    }
}

impl<Fp: FpTrait> MulAssign<&Polynomial<Fp>> for Polynomial<Fp> {
    #[inline(always)]
    fn mul_assign(&mut self, other: &Polynomial<Fp>) {
        self.set_mul(other);
    }
}

impl<Fp: FpTrait> MulAssign<Fp> for Polynomial<Fp> {
    #[inline(always)]
    fn mul_assign(&mut self, other: Fp) {
        self.scale_into(&other);
    }
}

impl<Fp: FpTrait> MulAssign<&Fp> for Polynomial<Fp> {
    #[inline(always)]
    fn mul_assign(&mut self, other: &Fp) {
        self.scale_into(other);
    }
}

impl<Fp: FpTrait> ::std::fmt::Display for Polynomial<Fp> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        for (i, c) in self.coeffs.iter().enumerate().rev() {
            if i == 0 {
                write!(f, "({c})")?
            } else {
                write!(f, "({c})*x^{i} + ")?
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test_poly_multieval {
    use crate::fields::sqisign::SqiField248Base as Fp;
    use crate::{
        polynomial_ring::poly::{EvalTree, Polynomial},
        utilities::test_utils::drng::DRNG,
    };

    type PR = Polynomial<Fp>;

    /// mul_low(f, g, n) must equal the first n coefficients of schoolbook(f, g).
    #[test]
    fn test_mul_low_matches_schoolbook() {
        let mut rng = DRNG::from_seed("mul_low".as_bytes());
        for _ in 0..20 {
            let f = PR::rand(&mut rng, 12);
            let g = PR::rand(&mut rng, 13);
            let full = f.basic_mul(&g); // reference schoolbook, length 24

            for &n in &[1usize, 2, 3, 5, 8, 13, 20, 24] {
                let mut h = vec![Fp::ZERO; n];
                Polynomial::mul_low(&mut h, &f.coeffs, &g.coeffs, n);
                for i in 0..n {
                    assert!(
                        h[i].equals(&full[i]) == u32::MAX,
                        "mul_low: coeff {i} wrong for n={n}"
                    );
                }
            }
        }
    }

    /// Larger inputs exercise the full recursive Karatsuba path.
    #[test]
    fn test_mul_low_large() {
        let mut rng = DRNG::from_seed("mul_low_large".as_bytes());
        for i in 0..5 {
            let f = PR::rand(&mut rng, 60 - i);
            let g = PR::rand(&mut rng, 60 - 2 * i);
            let full = f.basic_mul(&g);
            let n = 80;
            let mut h = vec![Fp::ZERO; n];
            Polynomial::mul_low(&mut h, &f.coeffs, &g.coeffs, n);
            for i in 0..n {
                assert!(
                    h[i].equals(&full[i]) == u32::MAX,
                    "mul_low large: coeff {i}"
                );
            }
        }
    }

    /// mul_middle(f, g) where f.len()==2n-1 must give (f*g)[n-1..2n-2].
    #[test]
    fn test_mul_middle_matches_schoolbook() {
        let mut rng = DRNG::from_seed("mul_middle".as_bytes());
        for _ in 0..20 {
            for &n in &[1usize, 2, 3, 4, 5, 8, 16, 32] {
                let f = PR::rand(&mut rng, 2 * n - 1);
                let g = PR::rand(&mut rng, n);
                let full = f.basic_mul(&g); // length 3n-2

                let mut h = vec![Fp::ZERO; n];
                Polynomial::mul_middle(&mut h, &f.coeffs, &g.coeffs);

                for i in 0..n {
                    assert!(
                        h[i].equals(&full[n - 1 + i]) == u32::MAX,
                        "mul_middle: coeff {i} wrong for n={n}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_mul_middle_large() {
        let mut rng = DRNG::from_seed("mul_middle_large".as_bytes());
        for &n in &[50usize, 64, 100, 128] {
            let f = PR::rand(&mut rng, 2 * n - 1);
            let g = PR::rand(&mut rng, n);
            let full = f.basic_mul(&g);
            let mut h = vec![Fp::ZERO; n];
            Polynomial::mul_middle(&mut h, &f.coeffs, &g.coeffs);
            for i in 0..n {
                assert!(h[i].equals(&full[n - 1 + i]) == u32::MAX, "n={n} coeff {i}");
            }
        }
    }

    #[test]
    fn test_mul_middle_rect_balanced() {
        // When in_n == out_n it must agree with mul_middle.
        let mut rng = DRNG::from_seed("mmr_balanced".as_bytes());
        for &n in &[1usize, 2, 3, 4, 8, 16, 32] {
            for _ in 0..5 {
                let f = PR::rand(&mut rng, 2 * n - 1);
                let g = PR::rand(&mut rng, n);
                let full = f.basic_mul(&g);

                let mut h_rect = vec![Fp::ZERO; n];
                let mut scratch = Vec::new();
                PR::mul_middle_rect(&mut h_rect, &f.coeffs, &g.coeffs, &mut scratch);

                for i in 0..n {
                    assert!(
                        h_rect[i].equals(&full[n - 1 + i]) == u32::MAX,
                        "mmr balanced n={n} coeff {i}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_mul_middle_rect_unbalanced() {
        // Various (in_n, out_n) pairs with in_n != out_n.
        let mut rng = DRNG::from_seed("mmr_unbalanced".as_bytes());
        for &(in_n, out_n) in &[
            (1usize, 2usize),
            (2, 1),
            (3, 5),
            (5, 3),
            (7, 8),
            (8, 7),
            (15, 16),
            (16, 15),
            (31, 32),
            (32, 31),
            (32, 64),
            (64, 32),
        ] {
            for _ in 0..5 {
                let f_len = in_n + out_n - 1;
                let f = PR::rand(&mut rng, f_len);
                let g = PR::rand(&mut rng, in_n);
                let full = f.basic_mul(&g);

                let mut h = vec![Fp::ZERO; out_n];
                let mut scratch = Vec::new();
                PR::mul_middle_rect(&mut h, &f.coeffs, &g.coeffs, &mut scratch);

                for i in 0..out_n {
                    assert!(
                        h[i].equals(&full[in_n - 1 + i]) == u32::MAX,
                        "mmr in_n={in_n} out_n={out_n} coeff {i}"
                    );
                }
            }
        }
    }

    /// (f * h) mod x^n == 1.
    #[test]
    fn test_inv_mod_xn() {
        let mut rng = DRNG::from_seed("inv_mod_xn".as_bytes());
        for _ in 0..20 {
            for &n in &[1usize, 2, 3, 4, 5, 8, 16, 32, 64] {
                let f = PR::rand(&mut rng, n + 3); // extra coefficients are clamped internally
                let mut h = vec![Fp::ZERO; n];
                Polynomial::inv_mod_xn(&mut h, &f.coeffs, n);

                let mut fh = vec![Fp::ZERO; n];
                Polynomial::mul_low(&mut fh, &f.coeffs, &h, n);

                assert!(
                    fh[0].equals(&Fp::ONE) == u32::MAX,
                    "inv: constant != 1, n={n}"
                );
                for i in 1..n {
                    assert!(fh[i].is_zero() == u32::MAX, "inv: coeff {i} nonzero, n={n}");
                }
            }
        }
    }

    /// Test that product_tree_flat produces the same root polynomial as
    /// product_tree_from_roots.
    #[test]
    fn test_product_tree_flat_matches_old() {
        let mut rng = DRNG::from_seed("flat_tree".as_bytes());
        for &n in &[1usize, 2, 3, 4, 5, 7, 8, 13, 16, 31, 32, 33, 64] {
            let roots: Vec<Fp> = (0..n).map(|_| Fp::rand(&mut rng)).collect();

            let old_tree = PR::product_tree_from_roots(&roots);
            let old_root = old_tree.last().unwrap()[0].clone();

            let (flat_tree, offsets) = PR::product_tree_flat(&roots);
            let n_layers = flat_tree.len();

            // The flat root layer has exactly n non-leading coefficients.
            assert_eq!(flat_tree[n_layers - 1].len(), n, "root layer length n={n}");
            assert_eq!(offsets[n_layers - 1], vec![0, n], "root offsets n={n}");

            // Compare non-leading coefficients against old tree root.
            // old_root = [c_0, c_1, ..., c_{n-1}, 1] (length n+1, monic)
            // flat root = [c_0, c_1, ..., c_{n-1}]   (length n)
            for i in 0..n {
                assert!(
                    flat_tree[n_layers - 1][i].equals(&old_root[i]) == u32::MAX,
                    "flat root coeff {i} wrong for n={n}"
                );
            }
        }
    }

    /// multieval must agree with Horner evaluation at every point.
    #[test]
    fn test_multieval_matches_horner() {
        let mut rng = DRNG::from_seed("multieval".as_bytes());
        for _ in 0..10 {
            for &(deg, n_roots) in &[(4usize, 8usize), (10, 16), (20, 32), (50, 64)] {
                let f = PR::rand(&mut rng, deg + 1);
                let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
                let multi_eval = f.multieval(&roots);
                let horner = f.multieval_horner(&roots);
                assert_eq!(multi_eval.len(), n_roots);
                for i in 0..n_roots {
                    assert!(
                        multi_eval[i].equals(&horner[i]) == u32::MAX,
                        "multieval: root {i} wrong (deg={deg}, n_roots={n_roots})"
                    );
                }
            }
        }
    }

    #[test]
    fn test_multieval_at_zero_and_one() {
        // f(0) = f[0], f(1) = sum of all coefficients — easy to verify.
        let mut rng = DRNG::from_seed("multieval_special".as_bytes());
        let f = PR::rand(&mut rng, 8);
        let vals = f.multieval(&[Fp::ZERO, Fp::ONE]);
        assert!(vals[0].equals(&f[0]) == u32::MAX, "multieval at 0");
        let sum = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
        assert!(vals[1].equals(&sum) == u32::MAX, "multieval at 1");
    }

    #[test]
    fn test_eval_tree_reuse() {
        // The primary VeluSqrt use case: one tree, four polynomials.
        let mut rng = DRNG::from_seed("etree_reuse".as_bytes());
        let n_roots = 64;
        let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
        let et = EvalTree::new(&roots);

        for _ in 0..4 {
            let f = PR::rand(&mut rng, 65);
            let fast = f.multieval_from_tree(&et);
            for (i, &r) in roots.iter().enumerate() {
                assert!(
                    fast[i].equals(&f.evaluate(&r)) == u32::MAX,
                    "reuse root {i} wrong"
                );
            }
        }
    }

    #[test]
    fn test_multieval_odd_counts() {
        // Non-power-of-two sizes exercise every carry path.
        let mut rng = DRNG::from_seed("fast_odd".as_bytes());
        for &n in &[3usize, 5, 6, 7, 9, 11, 13, 15, 17, 31, 33, 63, 65] {
            let f = PR::rand(&mut rng, n + 1);
            let roots: Vec<Fp> = (0..n).map(|_| Fp::rand(&mut rng)).collect();
            let et = EvalTree::new(&roots);
            let fast = f.multieval_from_tree(&et);
            for (i, &r) in roots.iter().enumerate() {
                assert!(
                    fast[i].equals(&f.evaluate(&r)) == u32::MAX,
                    "odd n={n} root {i} wrong"
                );
            }
        }
    }

    #[test]
    fn test_multieval_large() {
        let mut rng = DRNG::from_seed("fast_large".as_bytes());
        let f = PR::rand(&mut rng, 201);
        let roots: Vec<Fp> = (0..256).map(|_| Fp::rand(&mut rng)).collect();
        let et = EvalTree::new(&roots);
        let fast = f.multieval_from_tree(&et);
        for (i, &r) in roots.iter().enumerate() {
            assert!(
                fast[i].equals(&f.evaluate(&r)) == u32::MAX,
                "large root {i}"
            );
        }
    }

    #[test]
    fn test_multieval_special_values() {
        let mut rng = DRNG::from_seed("fast_special".as_bytes());
        let f = PR::rand(&mut rng, 8);
        let et = EvalTree::new(&[Fp::ZERO, Fp::ONE]);
        let vals = f.multieval_from_tree(&et);
        assert!(vals[0].equals(&f[0]) == u32::MAX, "f(0) wrong");
        let sum = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
        assert!(vals[1].equals(&sum) == u32::MAX, "f(1) wrong");
    }

    #[test]
    fn test_resultant_matches_horner() {
        let mut rng = DRNG::from_seed("resultant".as_bytes());
        for _ in 0..10 {
            for &(deg, n_roots) in &[(4usize, 8usize), (10, 16), (20, 32), (50, 64)] {
                let f = PR::rand(&mut rng, deg + 1);
                let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
                let res = f.resultant_from_roots(&roots);
                let horner = f.resultant_from_roots_horner(&roots);
                assert_eq!(res.equals(&horner), u32::MAX);
            }
        }
    }
}
