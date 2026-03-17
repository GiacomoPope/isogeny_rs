use fp2::traits::Fp as FpTrait;
use rand_core::{CryptoRng, RngCore};

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

// TODO: what's the optimum here? It seems to be surprisingly low!
const KARATSUBA_THRESHOLD: usize = 8;

/// Trait for arithmetic for univariate polynomials in Fp[X].
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
/// Stores a product tree in a flat layout so that multiple polynomials can be
/// evaluated at the same set of roots without rebuilding the tree each time.
///
/// tree[layer] is a contiguous Vec<Fp> of length exactly n_roots holding the
/// non-leading coefficients of every node at that layer. Because every layer is
/// the product tree of a degree-n monic polynomial, the sum of node degrees at
/// each layer is always exactly n, so n slots suffice for every layer.
///
/// offsets[layer] has n_nodes_at_layer + 1 entries. Node j at layer `layer`
/// occupies tree[layer][offsets[layer][j] .. offsets[layer][j+1]], so its
/// degree is offsets[layer][j+1] - offsets[layer][j].
pub struct EvalTree<Fp: fp2::traits::Fp> {
    /// tree[lyr][offsets[lyr][j]..offsets[lyr][j+1]] = non-leading coefficients
    /// of node j at layer lyr. Length of tree[lyr] is always n_roots.
    pub tree: Vec<Vec<Fp>>,
    pub rev_q_inv: Vec<Fp>,
    pub n_roots: usize,
    /// offsets[lyr][j] = start of node j's coefficient window in tree[lyr].
    pub offsets: Vec<Vec<usize>>,
}

impl<Fp: fp2::traits::Fp> EvalTree<Fp> {
    /// Build a product tree over `roots` and precompute rev(Q)^{-1} mod x^n,
    /// where Q = product(x - r_i). The inverse is used once per multieval call
    /// via a single mul_low, amortising the cost across all polynomial evaluations
    /// that share this tree.
    pub fn new(roots: &[Fp]) -> Self {
        use crate::polynomial_ring::poly::Polynomial;

        let n_roots = roots.len();
        if n_roots == 0 {
            return Self {
                tree: vec![],
                rev_q_inv: vec![],
                n_roots: 0,
                offsets: vec![],
            };
        }

        let (tree, offsets) = Polynomial::<Fp>::product_tree_from_roots(roots);
        let n_layers = tree.len();
        let n = n_roots;

        // Compute rev(Q)^{-1} mod x^n.
        // Q_root is monic of degree n; its non-leading coefficients are stored in
        // tree[n_layers-1][0..n]. rev(Q)[i] = Q[n-i], so:
        //   rev_q[0]   = Q[n]   = 1   (the leading term, reversed to position 0)
        //   rev_q[i]   = Q[n-i]       for i in 1..=n
        let root = &tree[n_layers - 1];
        let mut rev_q = vec![Fp::ZERO; n + 1];
        rev_q[0] = Fp::ONE;
        for i in 1..=n {
            rev_q[i] = root[n - i];
        }
        let mut rev_q_inv = vec![Fp::ZERO; n];
        Polynomial::<Fp>::inv_mod_xn(&mut rev_q_inv, &rev_q, n);

        Self {
            tree,
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

    /// Set the coefficients of a polynomial from a slice.
    pub fn set_from_slice(&mut self, a: &[Fp]) {
        self.coeffs.copy_from_slice(a)
    }

    /// The length of the coefficient vector. Note that trailing zeros are not
    /// trimmed automatically, so this may exceed the true degree plus one.
    fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Return the degree of the polynomial, or None if the polynomial is zero.
    ///
    /// Scans from the highest index downward to skip trailing zero coefficients.
    /// Trailing zeros in the coefficient vector are tolerated and do not mutate
    /// self.
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

    /// Return the constant coefficient, or None if the polynomial is empty.
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

    /// Return a new polynomial with the coefficients of self reversed.
    pub fn reverse(&self) -> Self {
        let mut r = self.clone();
        r.reverse_into();
        r
    }

    /// Return 0xFFFFFFFF if self and other represent the same polynomial,
    /// otherwise return 0x00000000.
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

    /// Return 0xFFFFFFFF if every coefficient is zero, otherwise 0x00000000.
    pub fn is_zero(&self) -> u32 {
        let mut is_zero = u32::MAX;
        for i in 0..self.len() {
            is_zero &= self.coeffs[i].is_zero();
        }
        is_zero
    }

    // Simple arithmetic functions ===============================================

    /// Compute f <-- f + g. Requires f.len() >= g.len().
    #[inline]
    fn add_into(f: &mut [Fp], g: &[Fp]) {
        debug_assert!(f.len() >= g.len());
        for i in 0..g.len() {
            f[i] += g[i];
        }
    }

    /// Compute f <-- f - g. Requires f.len() >= g.len().
    #[inline]
    fn sub_into(f: &mut [Fp], g: &[Fp]) {
        debug_assert!(f.len() >= g.len());
        for i in 0..g.len() {
            f[i] -= g[i];
        }
    }

    /// Compute f * g using schoolbook multiplication in O(n^2) Fp multiplications.
    /// Callers must zero fg before calling and ensure fg.len() >= f.len() + g.len() - 1.
    fn schoolbook_multiplication(fg: &mut [Fp], f: &[Fp], g: &[Fp]) {
        debug_assert!(fg.len() >= f.len() + g.len() - 1);
        for i in 0..f.len() {
            for j in 0..g.len() {
                // I could do this if I knew everything was zeroed in fg
                // fg[i + j] += f[i] * g[j];
                if i == 0 || j + 1 == g.len() {
                    fg[i + j] = f[i] * g[j]
                } else {
                    fg[i + j] += f[i] * g[j]
                }
            }
        }
    }

    /// Return the minimum scratch buffer length required by
    /// `karatsuba_multiplication_impl` for inputs of total length `n_f + n_g`
    /// where `n_f >= n_g`.
    ///
    /// The recurrence mirrors the scratch layout used by the algorithm:
    ///
    ///   T(n_f, n_g) = 2 * half + T(mf, max(nf, mg))
    ///
    /// where `half = n_f / 2`, `mf = n_f - half`.  For the common balanced case
    /// this resolves to roughly `2n`.
    pub fn karatsuba_scratch_size(n_f: usize, n_g: usize) -> usize {
        // Base cases that don't recurse into the scratch-buffer path.
        if n_g == 0 || n_g == 1 || n_f == 2 || n_f <= KARATSUBA_THRESHOLD {
            return 0;
        }

        // Ensure n_f >= n_g (mirrors the swap at the top of the impl).
        let (n_f, n_g) = if n_f < n_g { (n_g, n_f) } else { (n_f, n_g) };

        let nf = n_f / 2;
        let mf = n_f - nf;

        if n_g <= nf {
            // Unbalanced branch: recurse with (nf, n_g) and (mf, n_g).
            // The two recursive calls are sequential so they share the same
            // scratch space; take the max.
            let s1 = Self::karatsuba_scratch_size(nf, n_g);
            let s2 = Self::karatsuba_scratch_size(mf, n_g);
            return s1.max(s2);
        }

        let mg = n_g - nf;

        // Balanced branch scratch layout (see karatsuba_multiplication_impl):
        //   [0       .. mf]         f_mid  (f_lo + f_hi, length mf)
        //   [mf      .. mf+mg_max]  g_mid  (g_lo + g_hi, length max(nf,mg))
        //   [mf+mg_max ..]          scratch for the three recursive sub-calls
        //
        // The three sub-calls are (nf+nf-1, nf), (mf+mg-1, ?), (mf+mg_max-1, mg_max)
        // and they run sequentially so they share the suffix region.
        let mg_max = nf.max(mg);
        let header = mf + mg_max;
        let s_lo = Self::karatsuba_scratch_size(nf, nf);
        let s_hi = Self::karatsuba_scratch_size(mf + mg - 1, mf.max(mg)); // approx
        let s_mid = Self::karatsuba_scratch_size(mf, mg_max);
        let s_rec = s_lo.max(s_hi).max(s_mid);
        header + s_rec
    }

    /// Inner Karatsuba implementation that uses a caller-supplied scratch buffer
    /// instead of allocating.
    ///
    /// `tmp` must have length >= `karatsuba_scratch_size(f.len(), g.len())`.
    ///
    /// Scratch layout for the balanced recursive case (f.len() >= g.len(),
    /// nf = f.len()/2, mf = f.len()-nf, mg = g.len()-nf):
    ///
    ///   tmp[0      .. mf]         f_mid  = f_hi + f_lo[..mf] (length mf)
    ///   tmp[mf     .. mf+mg_max]  g_mid  = g_hi + g_lo[..mg_max] (length mg_max)
    ///   tmp[mf+mg_max ..]         scratch forwarded to recursive calls
    fn karatsuba_multiplication(fg: &mut [Fp], f: &[Fp], g: &[Fp], tmp: &mut [Fp]) {
        debug_assert!(fg.len() >= f.len() + g.len() - 1);

        // Ensure that the degree of f is larger or equal to g
        if f.len() < g.len() {
            Self::karatsuba_multiplication(fg, g, f, tmp);
            return;
        }

        if g.is_empty() {
            for c in fg.iter_mut() {
                *c = Fp::ZERO;
            }
            return;
        }

        if g.len() == 1 {
            let g0 = g[0];
            for i in 0..f.len() {
                fg[i] = f[i] * g0;
            }
            return;
        }

        if f.len() == 2 {
            // f and g are both linear
            let t1 = f[0] + f[1];
            let t2 = g[0] + g[1];
            fg[0] = f[0] * g[0];
            fg[2] = f[1] * g[1];
            fg[1] = t1 * t2;
            fg[1] -= fg[0];
            fg[1] -= fg[2];
            return;
        }

        if f.len() <= KARATSUBA_THRESHOLD {
            Self::schoolbook_multiplication(fg, f, g);
            return;
        }

        let nf = f.len() / 2;
        let mf = f.len() - nf;

        if g.len() <= nf {
            // Unbalanced: f is much longer than g.
            // Recurse with the same scratch buffer (calls are sequential).
            Self::karatsuba_multiplication(&mut fg[..nf + g.len() - 1], &f[..nf], g, tmp);

            let mut fg_hi = vec![Fp::ZERO; mf + g.len() - 1];
            Self::karatsuba_multiplication(&mut fg_hi, &f[nf..nf + mf], g, tmp);
            Self::add_into(&mut fg[nf..nf + g.len() - 1], &fg_hi[..g.len() - 1]);
            fg[nf + g.len() - 1..].copy_from_slice(&fg_hi[g.len() - 1..]);

            return;
        }

        let mg = g.len() - nf;
        let mg_max = nf.max(mg);

        // Scratch layout:
        //   tmp[0      .. mf]         f_mid
        //   tmp[mf     .. mf+mg_max]  g_mid
        //   tmp[mf+mg_max ..]         forwarded scratch for recursive calls
        debug_assert!(tmp.len() >= mf + mg_max);
        let (f_mid, rest) = tmp.split_at_mut(mf);
        let (g_mid, rec_scratch) = rest.split_at_mut(mg_max);

        // f_mid = f_hi (length mf), then add f_lo (length nf <= mf) into it.
        f_mid.copy_from_slice(&f[nf..]);
        Self::add_into(&mut f_mid[..nf], &f[..nf]);

        // g_mid = g_hi or g_lo depending on which is longer, then add the other.
        if mg < nf {
            g_mid[..nf].copy_from_slice(&g[..nf]);
            Self::add_into(&mut g_mid[..mg], &g[nf..]);
        } else {
            g_mid[..mg].copy_from_slice(&g[nf..]);
            Self::add_into(&mut g_mid[..nf], &g[..nf]);
        }

        // fg_lo = f_lo * g_lo  (written directly into fg[..2*nf-1])
        Self::karatsuba_multiplication(&mut fg[..nf + nf - 1], &f[..nf], &g[..nf], rec_scratch);

        // fg_hi = f_hi * g_hi  (written directly into fg[2*nf..])
        Self::karatsuba_multiplication(&mut fg[nf + nf..], &f[nf..], &g[nf..], rec_scratch);

        // fg_mid = f_mid * g_mid  (into a temporary on the stack, then added in)
        // We need mf + mg_max - 1 elements for this product.  We can't write it
        // directly into fg because the window overlaps both fg_lo and fg_hi.
        // We reuse rec_scratch prefix as the output buffer, but we need to save
        // and restore fg's overlap region first.
        //
        // Instead, allocate a small Vec here (this is the ONE remaining allocation
        // per top-level call, bounded in size and not inside a hot inner loop once
        // KARATSUBA_THRESHOLD kicks in the schoolbook path).
        let fg_mid_len = mf + mg_max - 1;
        let mut fg_mid = vec![Fp::ZERO; fg_mid_len];
        Self::karatsuba_multiplication(&mut fg_mid, f_mid, &g_mid[..mg_max], rec_scratch);

        // fg_mid -= fg_lo
        Self::sub_into(&mut fg_mid[..nf + nf - 1], &fg[..nf + nf - 1]);
        // fg_mid -= fg_hi
        Self::sub_into(&mut fg_mid[..mf + mg - 1], &fg[nf + nf..]);

        // Accumulate fg_mid into fg[nf..].
        fg[nf + nf - 1] = Fp::ZERO;
        Self::add_into(&mut fg[nf..nf + fg_mid_len], &fg_mid);
    }

    /// Compute f * g with ~O(n^1.5) Fp multiplications using Karatsuba multiplication.
    /// Assumes that fg has enough space for the result (len(f) + len(g) - 1).

    /// Compute the product of f with g into fg. Requires fg.len() >= f.len() + g.len() - 1.
    fn mul_into(fg: &mut [Fp], f: &[Fp], g: &[Fp]) {
        debug_assert!(f.len() + g.len() - 1 <= fg.len());

        let n_f = f.len();
        let n_g = g.len();

        // No scratch buffer needed
        if n_f.max(n_g) < KARATSUBA_THRESHOLD {
            Self::karatsuba_multiplication(fg, f, g, &mut []);
            return;
        }

        let scratch_len = Self::karatsuba_scratch_size(n_f.max(n_g), n_f.min(n_g));
        let mut tmp = vec![Fp::ZERO; scratch_len];
        Self::karatsuba_multiplication(fg, f, g, &mut tmp)
    }

    /// Set self to its negative.
    #[inline(always)]
    fn set_neg(&mut self) {
        for x in self.coeffs.iter_mut() {
            x.set_neg();
        }
    }

    /// Set self <-- self + other, extending self if necessary.
    fn set_add(&mut self, other: &Self) {
        if self.len() < other.len() {
            self.coeffs.resize(other.len(), Fp::ZERO);
        }
        Self::add_into(&mut self.coeffs, &other.coeffs);
    }

    /// Set self <-- self - other, extending self if necessary.
    fn set_sub(&mut self, other: &Self) {
        if self.len() < other.len() {
            self.coeffs.resize(other.len(), Fp::ZERO);
        }
        Self::sub_into(&mut self.coeffs, &other.coeffs);
    }

    /// Set self <-- self * other.
    fn set_mul(&mut self, other: &Self) {
        let mut fg_coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];
        Self::mul_into(&mut fg_coeffs, &self.coeffs, &other.coeffs);
        self.coeffs = fg_coeffs;
    }

    /// Multiply all coefficients of the polynomial by a field element, in place.
    pub fn scale_into(&mut self, c: &Fp) {
        for x in self.coeffs.iter_mut() {
            *x *= *c;
        }
    }

    /// Return c * self for some c in the finite field.
    pub fn scale(&self, c: &Fp) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.scale_into(c);
        r
    }

    /// Multiply all coefficients by a small integer, in place.
    pub fn scale_small_into(&mut self, k: i32) {
        for x in self.coeffs.iter_mut() {
            x.set_mul_small(k);
        }
    }

    /// Return k * self for some small integer k.
    pub fn scale_small(&self, k: i32) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.scale_small_into(k);
        r
    }

    /// Evaluate the polynomial at a using Horner's method in O(deg) multiplications.
    pub fn evaluate(&self, a: &Fp) -> Fp {
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

    /// Set all coefficients to random field elements.
    pub fn set_rand<R: CryptoRng + RngCore>(&mut self, rng: &mut R) {
        for x in self.coeffs.iter_mut() {
            x.set_rand(rng);
        }
    }

    /// Return a new polynomial of length d with random coefficients.
    pub fn rand<R: CryptoRng + RngCore>(rng: &mut R, d: usize) -> Self {
        let mut r = Self {
            coeffs: vec![Fp::ZERO; d],
        };
        r.set_rand(rng);
        r
    }

    // Specialised arithmetic functions =================================

    /// Multiply two monic polynomials, storing only the non-leading coefficients
    /// of the result.
    ///
    /// Given monic polynomials F of degree n_f and G of degree n_g, represented
    /// by their non-leading coefficients f_nl and g_nl, compute H = F * G and
    /// write its n_f + n_g non-leading coefficients into out.
    ///
    /// Uses scratch as a temporary buffer of length at least 2*(n_f+n_g)+3.
    /// The caller is responsible for sizing it; product_tree_from_roots
    /// allocates a single buffer of 2n+4 elements and reuses it for every node.
    fn monic_mul_nonleading(
        out: &mut [Fp],
        f_nl: &[Fp],
        n_f: usize,
        g_nl: &[Fp],
        n_g: usize,
        tmp: &mut [Fp],
    ) {
        debug_assert_eq!(f_nl.len(), n_f);
        debug_assert_eq!(g_nl.len(), n_g);
        debug_assert_eq!(out.len(), n_f + n_g);
        debug_assert!(tmp.len() >= 2 * (n_f + n_g) + 3);

        // F, G are linear: F = x + f0, G = x + g0.
        // H = x^2 + (f0+g0)*x + f0*g0.
        if n_f == 1 && n_g == 1 {
            out[0] = f_nl[0] * g_nl[0];
            out[1] = f_nl[0] + g_nl[0];
            return;
        }

        // F quadratic, G linear: F = x^2 + f1*x + f0, G = x + g0.
        // H = x^3 + (f1+g0)*x^2 + (f0+f1*g0)*x + f0*g0.
        // out = [f0*g0, f0+f1*g0, f1+g0].  Two multiplications.
        if n_f == 2 && n_g == 1 {
            let t = f_nl[1] * g_nl[0]; // f1*g0
            out[0] = f_nl[0] * g_nl[0]; // f0*g0
            out[1] = f_nl[0] + t; // f0 + f1*g0
            out[2] = f_nl[1] + g_nl[0]; // f1 + g0
            return;
        }

        // F, Q quadratic: F = x^2 + f1*x + f0, G = x^2 + g1*x + g0.
        //
        // H = x^4 + (f1+g1)*x^3 + (f0+g0+f1*g1)*x^2 + (f0*g1+f1*g0)*x + f0*g0.
        // out = [f0*g0, f0*g1+f1*g0, f0+g0+f1*g1, f1+g1].
        // Three multiplications via Karatsuba on the non-leading parts.
        if n_f == 2 && n_g == 2 {
            let ll = f_nl[0] * g_nl[0]; // f0*g0
            let hh = f_nl[1] * g_nl[1]; // f1*g1
            let mid = (f_nl[0] + f_nl[1]) * (g_nl[0] + g_nl[1]) - ll - hh;
            out[0] = ll;
            out[1] = mid;
            out[2] = f_nl[0] + g_nl[0] + hh; // f0 + g0 + f1*g1
            out[3] = f_nl[1] + g_nl[1];
            return;
        }

        // Build [f_nl | 1] and [g_nl | 1] in scratch, then multiply into scratch,
        // then copy the first n_f + n_g coefficients to out (discarding the
        // leading 1 of the product at position n_f + n_g).
        //
        // Scratch layout:
        //   [0          .. n_f+1]          f_nl ++ [1]
        //   [n_f+1      .. n_f+n_g+2]      g_nl ++ [1]
        //   [n_f+n_g+2  .. 2*(n_f+n_g)+3]  product (length n_f+n_g+1)

        let f_start = 0;
        let f_end = n_f + 1;
        let g_start = f_end;
        let g_end = g_start + n_g + 1;
        let p_start = g_end;
        let p_end = p_start + n_f + n_g + 1;

        tmp[f_start..f_start + n_f].copy_from_slice(f_nl);
        tmp[f_start + n_f] = Fp::ONE;

        tmp[g_start..g_start + n_g].copy_from_slice(g_nl);
        tmp[g_start + n_g] = Fp::ONE;

        for c in tmp[p_start..p_end].iter_mut() {
            *c = Fp::ZERO;
        }

        let (inputs, product_region) = tmp.split_at_mut(p_start);
        // TODO: remove the allocation here
        Self::mul_into(
            &mut product_region[..n_f + n_g + 1],
            &inputs[f_start..f_end],
            &inputs[g_start..g_end],
        );

        out.copy_from_slice(&tmp[p_start..p_start + n_f + n_g]);
    }

    /// Return the minimum scratch buffer length required by `mul_low_impl` for
    /// inputs clamped to length n.
    ///
    /// Scratch layout per recursive level (m = ceil(n/2), rem = n-m,
    /// ll_len = 2m-1 in the worst case):
    ///
    ///   [0        .. 2m-1]       ll    – full product of f_lo * g_lo
    ///   [2m-1     .. 2m-1+rem]   cross – reused for both cross terms sequentially
    ///   [2m-1+rem ..]            forwarded to recursive sub-calls
    ///
    /// Recurrence: T(n) = (2m-1) + rem + T(m),  T(n) = 0 when n <= THRESHOLD.
    pub fn mul_low_scratch_size(n: usize) -> usize {
        if n <= KARATSUBA_THRESHOLD {
            return 0;
        }
        let m = (n + 1) / 2;
        let rem = n - m;
        let ll_max = 2 * m - 1; // worst-case length of f_lo*g_lo
        ll_max + rem + Self::mul_low_scratch_size(m)
    }

    /// Inner implementation of mul_low that uses a caller-supplied scratch buffer.
    ///
    /// `tmp` must have length >= `mul_low_scratch_size(n)`.
    fn mul_low_impl(h: &mut [Fp], f: &[Fp], g: &[Fp], n: usize, tmp: &mut [Fp]) {
        debug_assert_eq!(h.len(), n);

        if n == 0 {
            return;
        }

        let f = &f[..f.len().min(n)];
        let g = &g[..g.len().min(n)];

        if f.is_empty() || g.is_empty() {
            for c in h.iter_mut() {
                *c = Fp::ZERO;
            }
            return;
        }

        // Schoolbook base case – no scratch needed.
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

        let m = (n + 1) / 2;
        let rem = n - m;

        let (f_lo, f_hi) = f.split_at(f.len().min(m));
        let (g_lo, g_hi) = g.split_at(g.len().min(m));

        // Carve scratch:
        //   ll    = tmp[0     .. ll_len]          (f_lo * g_lo, worst-case 2m-1)
        //   cross = tmp[2m-1  .. 2m-1+rem]        (reused for both cross terms)
        //   rec   = tmp[2m-1+rem ..]              (forwarded to recursive calls)
        let ll_len = f_lo.len() + g_lo.len() - 1;
        let ll_slot = 2 * m - 1; // fixed-size slot so cross starts at a stable offset
        let (ll_and_cross, rec_scratch) = tmp.split_at_mut(ll_slot + rem);
        let (ll, cross) = ll_and_cross.split_at_mut(ll_slot);

        // f_lo * g_lo into ll[..ll_len]; the rest of the ll slot is unused padding.
        for c in ll[..ll_len].iter_mut() {
            *c = Fp::ZERO;
        }
        Self::karatsuba_multiplication(&mut ll[..ll_len], f_lo, g_lo, rec_scratch);

        // Distribute ll into h.
        let copy_lo = ll_len.min(m);
        h[..copy_lo].copy_from_slice(&ll[..copy_lo]);
        for c in h[copy_lo..m].iter_mut() {
            *c = Fp::ZERO;
        }

        if rem == 0 {
            return;
        }

        for c in h[m..].iter_mut() {
            *c = Fp::ZERO;
        }
        if ll_len > m {
            let overflow = (ll_len - m).min(rem);
            h[m..m + overflow].copy_from_slice(&ll[m..m + overflow]);
        }

        // Two cross terms, computed sequentially into the same `cross` buffer.
        if !g_hi.is_empty() {
            for c in cross.iter_mut() {
                *c = Fp::ZERO;
            }
            Self::mul_low_impl(cross, f_lo, g_hi, rem, rec_scratch);
            Self::add_into(&mut h[m..], cross);
        }

        if !f_hi.is_empty() {
            for c in cross.iter_mut() {
                *c = Fp::ZERO;
            }
            Self::mul_low_impl(cross, f_hi, g_lo, rem, rec_scratch);
            Self::add_into(&mut h[m..], cross);
        }
    }

    /// Compute only the low n coefficients of the product f * g.
    ///
    /// Uses a Karatsuba-style halving recursion with split point m = ceil(n/2):
    ///
    ///   f = f_lo + x^m * f_hi,   g = g_lo + x^m * g_hi
    ///
    /// The low n coefficients of f * g are:
    ///
    ///   (f_lo * g_lo)[0..n]               full sub-product, truncated at x^(n-1)
    ///   + x^m * (f_lo * g_hi)[0..rem]     cross term
    ///   + x^m * (f_hi * g_lo)[0..rem]     cross term
    ///
    /// The f_hi * g_hi term starts at x^(2m) >= x^n and is discarded.
    ///
    /// h must have length exactly n.
    fn mul_low(h: &mut [Fp], f: &[Fp], g: &[Fp], n: usize) {
        if n <= KARATSUBA_THRESHOLD {
            // Fast path: schoolbook, no scratch needed.
            Self::mul_low_impl(h, f, g, n, &mut []);
            return;
        }
        let scratch_len = Self::mul_low_scratch_size(n);
        let mut tmp = vec![Fp::ZERO; scratch_len];
        Self::mul_low_impl(h, f, g, n, &mut tmp);
    }

    /// Compute the middle n coefficients of f * g, i.e. h[i] = (f*g)[n-1+i].
    ///
    /// Requires f.len() == 2*n - 1 and g.len() == n.
    ///
    /// Uses the Hanrot-Quercia-Zimmermann (HQZ) 3-way recursive algorithm.
    /// Splitting at half = floor(n/2) and half_up = ceil(n/2):
    ///
    ///   a = middle(f + x^{half_up}*f,             g[half..])          length half_up
    ///   c = middle(shifted window of the same sum, g[..half])          length half
    ///   b = middle(f[half_up..3*half_up-1],        g[half..]-g[..half]) length half_up
    ///
    ///   h[0..half_up]  = a - b
    ///   h[half_up..n]  = c + b[0..half]
    ///
    /// The tmp slice is scratch space that avoids heap allocation inside the
    /// recursion. It must have length at least mul_middle_scratch_size(n).
    ///
    /// Scratch layout per call level, carved with split_at_mut:
    ///
    ///   [0                .. 2n-1]          tmplo  -- modified copy of f
    ///   [2n-1             .. 2n-1+half_up]  g_diff -- g[half..] - g[..half]
    ///   [2n-1+half_up     .. 2n-1+2*half_up] b     -- output of the b sub-call
    ///   [2n-1+2*half_up   ..]               recursive scratch for all sub-calls
    ///
    /// The a and c sub-calls share the region from 2n-1+half_up onward as
    /// scratch, running sequentially so they reuse the same memory. The b
    /// sub-call carves off half_up elements for its output first, then uses
    /// the remainder for its own recursive scratch.
    fn mul_middle(h: &mut [Fp], f: &[Fp], g: &[Fp], tmp: &mut [Fp]) {
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
            // Explicit 2-coefficient case avoids the recursion overhead.
            // (f*g)[1] = f[0]*g[1] + f[1]*g[0]
            // (f*g)[2] = f[1]*g[1] + f[2]*g[0]
            h[0] = f[0] * g[1] + f[1] * g[0];
            h[1] = f[1] * g[1] + f[2] * g[0];
            return;
        }

        let half = n / 2; // floor(n/2)
        let half_up = n - half; // ceil(n/2)

        // Carve the scratch buffer into regions for this call level.
        // tmphi is the region passed as scratch to all recursive sub-calls.
        let (tmplo, rest) = tmp.split_at_mut(2 * n - 1);
        let (g_diff, tmphi) = rest.split_at_mut(half_up);

        // Build tmplo = f + x^{half_up} * f, i.e. tmplo[i] += f[i + half_up].
        // This constructs the two shifted sums (P2*P1 + P1*P0) and (P3*P2 + P2*P1)
        // from HQZ in a single pass.
        tmplo.copy_from_slice(f);
        let add_len = f.len() - half_up;
        Self::add_into(&mut tmplo[..add_len], &f[half_up..]);

        // Zero h before writing a and c into it.
        for x in h.iter_mut() {
            *x = Fp::ZERO;
        }

        // Compute a into h[0..half_up]. tmphi is reused as scratch for a and c
        // since they run sequentially.
        Self::mul_middle(
            &mut h[..half_up],
            &tmplo[..2 * half_up - 1],
            &g[half..],
            tmphi,
        );

        // Compute c into h[half_up..n].
        if half > 0 {
            Self::mul_middle(
                &mut h[half_up..],
                &tmplo[half_up..half_up + 2 * half - 1],
                &g[..half],
                tmphi,
            );
        }

        // Build g_diff = g[half..] with g[..half] subtracted from the aligned
        // portion. For even n, g_diff[i] = g[half+i] - g[i]. For odd n, the
        // subtraction is shifted right by one: g_diff[0] = g[half], then
        // g_diff[i] = g[half+i] - g[i-1] for i in 1..half_up.
        g_diff.copy_from_slice(&g[half..]);
        Self::sub_into(&mut g_diff[half_up - half..], &g[..half]);

        // Compute b into tmphi[0..half_up] so its output survives past the
        // recursive call. The remainder of tmphi is used as scratch for b's
        // own sub-calls.
        let (b, tmphi2) = tmphi.split_at_mut(half_up);
        for x in b.iter_mut() {
            *x = Fp::ZERO;
        }
        Self::mul_middle(b, &f[half_up..3 * half_up - 1], g_diff, tmphi2);

        // Combine: h[0..half_up] = a - b,  h[half_up..n] = c + b[0..half].
        for i in 0..half_up {
            h[i] -= b[i];
        }
        for i in 0..half {
            h[half_up + i] += b[i];
        }
    }

    /// Return the minimum scratch buffer length required by mul_middle for an
    /// input g of length n (equivalently, f of length 2n-1).
    ///
    /// The recurrence reflects the scratch layout:
    ///
    ///   T(n) = (2n-1) + 2*half_up + T(half_up),  T(1) = T(2) = 0
    ///
    /// where (2n-1) is for tmplo, the first half_up is for g_diff, the second
    /// half_up is for the b output buffer, and T(half_up) is the recursive
    /// scratch needed by the deepest sub-call.
    pub fn mul_middle_scratch_size(n: usize) -> usize {
        if n <= 2 {
            // Even for small n, the base case in mul_middle_rect may compute a
            // full product of two polynomials of combined length up to 2n-1+n-1
            // = 3n-2. Return that as a safe floor.
            return if n == 0 { 0 } else { 3 * n };
        }
        let half_up = n - n / 2;
        (2 * n - 1) + 2 * half_up + Self::mul_middle_scratch_size(half_up)
    }

    /// Compute the middle out_n coefficients of f * g, i.e. h[i] = (f*g)[in_n-1+i].
    ///
    /// Requires f.len() == in_n + out_n - 1 and g.len() == in_n.
    ///
    /// When in_n == out_n this delegates directly to mul_middle. Otherwise the
    /// larger of the two dimensions is halved each call, so in_n and out_n
    /// converge to equality and mul_middle is always reached.
    ///
    /// When out_n > in_n we split out_n:
    ///   h[0..half]     = MP(f[..in_n+half-1], g)   middle offset in_n-1
    ///   h[half..out_n] = MP(f[half..],         g)   middle offset in_n-1
    ///
    /// When out_n < in_n we split in_n at half = floor(in_n/2):
    ///   g = g_lo(len=half) + x^half * g_hi(len=half_up)
    ///   (f*g)[in_n-1+i] = (f*g_lo)[half_up+half-1+i] + (f*g_hi)[half_up-1+i]
    ///   g_hi: f[..half_up+out_n-1], in_n=half_up  (into h)
    ///   g_lo: f[half_up..half_up+half+out_n-1], in_n=half  (accumulated into h)
    ///
    /// The scratch slice must have length at least mul_middle_scratch_size(max(in_n,out_n)).
    fn mul_middle_rect(h: &mut [Fp], f: &[Fp], g: &[Fp], scratch: &mut [Fp]) {
        let out_n = h.len();
        let in_n = g.len();
        debug_assert_eq!(f.len(), in_n + out_n - 1);

        if in_n == 0 || out_n == 0 {
            return;
        }

        // Balanced: hand off to the optimised HQZ middle product.
        if in_n == out_n {
            Self::mul_middle(h, f, g, scratch);
            return;
        }

        // Schoolbook base case: when either dimension is small enough,
        // a direct loop over the relevant coefficients is cheapest.
        // We only compute the out_n coefficients at offset in_n-1.
        if in_n <= KARATSUBA_THRESHOLD || out_n <= KARATSUBA_THRESHOLD {
            for c in h.iter_mut() {
                *c = Fp::ZERO;
            }
            // h[k] = sum_{i+j = in_n-1+k} f[i]*g[j]
            //      = sum_{j=0}^{in_n-1} f[in_n-1+k-j] * g[j]   for k in 0..out_n
            for k in 0..out_n {
                let offset = in_n - 1 + k;
                // j ranges so that both f[offset-j] and g[j] are in bounds
                let j_lo = offset.saturating_sub(f.len() - 1);
                let j_hi = offset.min(in_n - 1);
                for j in j_lo..=j_hi {
                    h[k] += f[offset - j] * g[j];
                }
            }
            return;
        }

        if out_n > in_n {
            let half = out_n / 2;
            Self::mul_middle_rect(&mut h[..half], &f[..in_n + half - 1], g, scratch);
            Self::mul_middle_rect(&mut h[half..], &f[half..], g, scratch);
        } else {
            let half = in_n / 2;
            let half_up = in_n - half;

            let (tmp, rec_scratch) = scratch.split_at_mut(out_n);

            for c in h.iter_mut() {
                *c = Fp::ZERO;
            }
            Self::mul_middle_rect(h, &f[..half_up + out_n - 1], &g[half..], rec_scratch);

            for c in tmp.iter_mut() {
                *c = Fp::ZERO;
            }
            Self::mul_middle_rect(
                tmp,
                &f[half_up..half_up + half + out_n - 1],
                &g[..half],
                rec_scratch,
            );
            Self::add_into(h, tmp);
        }
    }

    pub fn mul_middle_rect_scratch_size(in_n: usize, out_n: usize) -> usize {
        if in_n == 0 || out_n == 0 {
            return 0;
        }
        if in_n == out_n {
            return Self::mul_middle_scratch_size(in_n);
        }
        // Base case: schoolbook, no scratch needed.
        if in_n <= KARATSUBA_THRESHOLD || out_n <= KARATSUBA_THRESHOLD {
            return 0;
        }
        if out_n > in_n {
            let half = out_n / 2;
            let s1 = Self::mul_middle_rect_scratch_size(in_n, half);
            let s2 = Self::mul_middle_rect_scratch_size(in_n, out_n - half);
            s1.max(s2)
        } else {
            let half = in_n / 2;
            let half_up = in_n - half;
            let s1 = Self::mul_middle_rect_scratch_size(half_up, out_n);
            let s2 = Self::mul_middle_rect_scratch_size(half, out_n);
            out_n + s1.max(s2)
        }
    }

    /// Return the minimum scratch buffer length required by `inv_mod_xn` for
    /// an input of length n.
    ///
    /// Scratch layout per recursive level (half_up = ceil(n/2), half = floor(n/2)):
    ///
    ///   [0                    .. 2*half_up-1]           f_shifted
    ///   [2*half_up-1          .. 3*half_up-1]           e
    ///   [3*half_up-1          .. 3*half_up-1+half]      corr
    ///   [3*half_up-1+half     ..]                       mm_scratch (mul_middle)
    ///                                                   + recursive scratch for
    ///                                                     inv_mod_xn and mul_low
    ///
    /// The recursive inv_mod_xn call operates on half_up elements, and the
    /// mul_middle and mul_low calls also operate on at most half_up elements,
    /// so they all share the suffix region sequentially.
    pub fn inv_mod_xn_scratch_size(n: usize) -> usize {
        if n <= 1 {
            return 0; // n==0 or n==1 never reaches scratch carving
        }
        let half = n / 2;
        let half_up = n - half;
        let header = (2 * half_up - 1) + half_up + half; // f_shifted + e + corr
        let mm = Self::mul_middle_scratch_size(half_up);
        let rec_inv = Self::inv_mod_xn_scratch_size(half_up);
        let rec_low = Self::mul_low_scratch_size(half);
        let suffix = mm.max(rec_inv).max(rec_low);
        header + suffix
    }

    /// Compute h[0..n] such that (f * h) == 1 mod x^n.
    ///
    /// Requires f[0] to be invertible. Uses the Hanrot-Quercia-Zimmermann
    /// lifting algorithm (section 5.1): compute an approximation to half the
    /// precision, then lift it to full precision in one step.
    ///
    ///   half     = floor(n/2)
    ///   half_up  = ceil(n/2)
    ///
    ///   1. Recurse to get h_k[0..half_up] with f * h_k == 1 mod x^{half_up}.
    ///   2. Compute the error e = (f * h_k)[half_up..n+half_up-1] via a middle
    ///      product of f[1..] (zero-padded to length 2*half_up-1) against h_k.
    ///   3. Correct: h[half_up..n] = -(h_k[0..half] * e[0..half])[0..half].
    ///
    /// h must have length exactly n.
    pub fn inv_mod_xn(h: &mut [Fp], f: &[Fp], n: usize) {
        debug_assert_eq!(h.len(), n);
        if n == 0 {
            return;
        }
        let scratch_len = Self::inv_mod_xn_scratch_size(n);
        let mut tmp = vec![Fp::ZERO; scratch_len];
        Self::inv_mod_xn_impl(h, f, n, &mut tmp);
    }

    fn inv_mod_xn_impl(h: &mut [Fp], f: &[Fp], n: usize, tmp: &mut [Fp]) {
        debug_assert_eq!(h.len(), n);

        if n == 0 {
            return;
        }

        // Base case: always correct regardless of f[0].
        if n == 1 {
            h[0] = Fp::ONE / f[0];
            return;
        }

        // Fast paths for n=2,3 when f[0] == 1, avoiding the general path which
        // requires scratch. When f[0] != 1 we fall through to the general path;
        // inv_mod_xn_scratch_size accounts for scratch at n=2 and n=3.
        if f[0].equals(&Fp::ONE) == u32::MAX {
            match n {
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
        }

        let half = n / 2;
        let half_up = n - half;

        // Scratch layout:
        //   f_shifted = tmp[0               .. 2*half_up-1]
        //   e         = tmp[2*half_up-1     .. 3*half_up-1]
        //   corr      = tmp[3*half_up-1     .. 3*half_up-1+half]
        //   suffix    = tmp[3*half_up-1+half ..]   shared by mm_scratch / rec calls
        let f_shift_end = 2 * half_up - 1;
        let e_end = f_shift_end + half_up;
        let corr_end = e_end + half;
        let (f_shifted, rest) = tmp.split_at_mut(f_shift_end);
        let (e, rest2) = rest.split_at_mut(half_up);
        let (corr, suffix) = rest2.split_at_mut(half);

        // Step 1: recurse to get h_k mod x^{half_up}.
        Self::inv_mod_xn_impl(&mut h[..half_up], &f[..half_up], half_up, suffix);

        // Step 2: build f_shifted = f[1..1+f_shift_len], zero-padded.
        let f_shift_len = 2 * half_up - 1;
        for c in f_shifted.iter_mut() {
            *c = Fp::ZERO;
        }
        let copy_end = (1 + f_shift_len).min(f.len());
        if copy_end > 1 {
            f_shifted[..copy_end - 1].copy_from_slice(&f[1..copy_end]);
        }

        // Step 2 cont: middle product to get the error e.
        for c in e.iter_mut() {
            *c = Fp::ZERO;
        }
        Self::mul_middle(e, f_shifted, &h[..half_up], suffix);

        // Step 3: correction via mul_low, result negated into h[half_up..n].
        for c in corr.iter_mut() {
            *c = Fp::ZERO;
        }
        Self::mul_low_impl(corr, &h[..half], &e[..half], half, suffix);
        for i in 0..half {
            h[half_up + i] = -corr[i];
        }

        // suppress unused-variable warnings from the size calculation
        let _ = corr_end;
    }

    /// Compute the product of all the quadratic polynomials in leaves, returning
    /// the root of the product tree without storing the intermediate layers.
    ///
    /// Each leaf is stored as [c0, c1, c2] representing c0 + c1*x + c2*x^2.
    /// The result has degree 2 * leaves.len() and length 2 * leaves.len() + 1.
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

        Self::new_from_slice(&buf_in[..2 * n + 1])
    }

    /// Build a product tree over roots using a simple node-per-Vec layout.
    ///
    /// Kept for testing purposes. Prefer product_tree_from_roots for production
    /// use as this version performs significantly more heap allocation.
    pub fn product_tree_from_roots_simple(roots: &[Fp]) -> Vec<Vec<Vec<Fp>>> {
        let leaves: Vec<Vec<Fp>> = roots.iter().map(|r| vec![-*r, Fp::ONE]).collect();
        let mut layers: Vec<Vec<Vec<Fp>>> = vec![leaves];

        loop {
            let prev = layers.last().unwrap();
            if prev.len() == 1 {
                break;
            }

            let mut next: Vec<Vec<Fp>> = Vec::with_capacity((prev.len() + 1) / 2);
            let mut i = 0;
            while i < prev.len() {
                if i + 1 < prev.len() {
                    let f = &prev[i];
                    let g = &prev[i + 1];
                    let mut fg = vec![Fp::ZERO; f.len() + g.len() - 1];
                    Self::mul_into(&mut fg, f, g);
                    next.push(fg);
                    i += 2;
                } else {
                    next.push(prev[i].clone());
                    i += 1;
                }
            }
            layers.push(next);
        }

        layers
    }

    /// Build a product tree over roots in a flat memory layout.
    ///
    /// Every node is a monic polynomial stored by its non-leading coefficients
    /// only. Since all node degrees at a given layer sum to n (the total number
    /// of roots), each layer fits in exactly n contiguous Fp elements.
    ///
    /// Returns (flat_tree, offsets) where:
    ///   flat_tree[lyr][offsets[lyr][j] .. offsets[lyr][j+1]]
    /// holds the non-leading coefficients of node j at layer lyr.
    ///
    /// A single scratch buffer of size <= 2n + 4 is allocated once and reused
    /// across all node multiplications, so there is at most one allocation per
    /// layer beyond the layer buffer itself.
    pub fn product_tree_from_roots(roots: &[Fp]) -> (Vec<Vec<Fp>>, Vec<Vec<usize>>) {
        let n = roots.len();
        debug_assert!(n > 0);

        // Layer 0: n linear monic polynomials (x - r_i).
        // Non-leading part of (x - r_i) is [-r_i], a single element.
        // flat_tree[0] = [-r_0, -r_1, ..., -r_{n-1}], length n.
        // offsets[0]   = [0, 1, 2, ..., n].
        let leaf_layer: Vec<Fp> = roots.iter().map(|r| -*r).collect();
        let leaf_offsets: Vec<usize> = (0..=n).collect();

        let mut flat_tree: Vec<Vec<Fp>> = vec![leaf_layer];
        let mut offsets: Vec<Vec<usize>> = vec![leaf_offsets];

        // Scratch buffer for monic_mul_nonleading. The worst case is at the root
        // layer where the two children have combined degree n, requiring
        // 2*(n_f+n_g)+3 = 2n+3 elements. Allocate 2n+4 as a safe upper bound.
        let mut scratch: Vec<Fp> = vec![Fp::ZERO; 2 * n + 4];

        loop {
            let prev_off = offsets.last().unwrap();
            let n_prev_nodes = prev_off.len() - 1;
            if n_prev_nodes == 1 {
                break;
            }

            let prev_flat = flat_tree.last().unwrap();

            let mut next_flat = vec![Fp::ZERO; n];
            let mut next_off: Vec<usize> = Vec::with_capacity((n_prev_nodes + 1) / 2 + 1);
            next_off.push(0);

            let mut out_cursor = 0usize;
            let mut i = 0usize;

            while i < n_prev_nodes {
                if i + 1 < n_prev_nodes {
                    let l_start = prev_off[i];
                    let l_end = prev_off[i + 1];
                    let r_start = prev_off[i + 1];
                    let r_end = prev_off[i + 2];

                    let n_l = l_end - l_start;
                    let n_r = r_end - r_start;
                    let n_out = n_l + n_r;

                    Self::monic_mul_nonleading(
                        &mut next_flat[out_cursor..out_cursor + n_out],
                        &prev_flat[l_start..l_end],
                        n_l,
                        &prev_flat[r_start..r_end],
                        n_r,
                        &mut scratch,
                    );

                    out_cursor += n_out;
                    next_off.push(out_cursor);
                    i += 2;
                } else {
                    // Odd node out: carry it unchanged into the next layer.
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

    /// Evaluate self at each root in eval_tree using a scaled remainder tree.
    ///
    /// The algorithm computes the Laurent series F = rev(P) * rev(Q)^{-1} mod x^n
    /// and then descends the product tree, splitting F at each node into the two
    /// child remainders via middle products. Each leaf holds P(r_i).
    ///
    /// Building an EvalTree once and reusing it via this method amortises the
    /// cost of the product tree construction and the power series inversion when
    /// the same set of roots is used multiple times.
    pub fn multieval_from_tree(&self, eval_tree: &EvalTree<Fp>) -> Vec<Fp> {
        let n_roots = eval_tree.n_roots;
        if n_roots == 0 {
            return vec![];
        }

        let flat_tree = &eval_tree.tree;
        let offsets = &eval_tree.offsets;
        let n_layers = flat_tree.len();
        let n = n_roots;

        let degp_full = match self.degree() {
            Some(d) => d,
            None => return vec![Fp::ZERO; n_roots],
        };

        // Reduce P mod Q_root if deg(P) >= n so the scaled representation is valid.
        // In VeluSqrt deg(P) = 2*size_J < size_I = n so this branch is dead in practice.
        let owned: Option<Vec<Fp>> = if degp_full >= n {
            let top = &flat_tree[n_layers - 1];
            let mut rem = self.coeffs.clone();
            rem.resize(degp_full + 1, Fp::ZERO);
            for i in (n..=degp_full).rev() {
                if rem[i].is_zero() != u32::MAX {
                    let c = rem[i];
                    for j in 0..n {
                        rem[i - n + j] -= c * top[j];
                    }
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

        // Compute the scaled series: rev(P) * rev(Q)^{-1} mod x^n.
        let mut series = vec![Fp::ZERO; n];
        Self::mul_low(&mut series, &rev_p, &eval_tree.rev_q_inv, n);

        // Reverse the series back to polynomial order for the tree descent.
        let mut node_buf = vec![Fp::ZERO; n];
        for i in 0..=degp {
            node_buf[i] = series[degp - i];
        }

        let mut dst_buf = vec![Fp::ZERO; n];

        // Allocate scratch for mul_middle_rect once and reuse across all calls.
        let scratch_size = Self::mul_middle_rect_scratch_size(n, n);
        let mut mmr_scratch: Vec<Fp> = vec![Fp::ZERO; scratch_size];

        // Descend the tree layer by layer.
        //
        // At layer lyr+1 (parent), node j occupies node_buf[off_p[j]..off_p[j+1]].
        // For each genuine pair of children (non-carry):
        //
        //   parent_tail = parent[1..]              length n_l + n_r - 1
        //   q_l_nl = non-leading coeffs of left child
        //   q_r_nl = non-leading coeffs of right child
        //
        //   left  child = mul_middle_rect(parent_tail, q_r_nl) + parent[..n_l]
        //   right child = mul_middle_rect(parent_tail, q_l_nl) + parent[..n_r]
        //
        // For a carry node, copy the parent slice directly into the child slot.
        for lyr in (0..n_layers - 1).rev() {
            let off_p = &offsets[lyr + 1];
            let off_c = &offsets[lyr];
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
                    let parent_tail = &parent[1..];

                    {
                        let out = &mut dst_buf[off_c[ci]..off_c[ci + 1]];
                        if n_r > 0 {
                            Self::mul_middle_rect(out, parent_tail, q_r_nl, &mut mmr_scratch);
                        }
                        Self::add_into(out, &parent[..n_l]);
                    }
                    {
                        let out = &mut dst_buf[off_c[ci + 1]..off_c[ci + 2]];
                        if n_l > 0 {
                            Self::mul_middle_rect(out, parent_tail, q_l_nl, &mut mmr_scratch);
                        }
                        Self::add_into(out, &parent[..n_r]);
                    }
                    child_idx += 2;
                } else {
                    let ci = child_idx;
                    dst_buf[off_c[ci]..off_c[ci + 1]].copy_from_slice(&node_buf[ps..pe]);
                    child_idx += 1;
                }
            }

            std::mem::swap(&mut node_buf, &mut dst_buf);
            // dst_buf does not need to be zeroed here: every slot is overwritten
            // before being read at the next layer (mul_middle_rect writes before
            // add_into reads, and carry nodes copy unconditionally).
        }

        // At layer 0 each node has exactly 1 coefficient equal to P(r_i).
        node_buf[..n_roots].to_vec()
    }

    // TODO this could be refactored into the above potentially, there's a lot of code duplication
    /// Evaluate self at a perfectly balanced power-of-2 set of roots.
    ///
    /// Uses a simplified descent that exploits the uniform stride layout:
    /// at layer lyr, node j occupies tree[lyr][j*(1<<lyr) .. (j+1)*(1<<lyr)].
    /// This avoids all offset lookups and calls mul_middle directly (never
    /// mul_middle_rect) since all sibling pairs have equal size.
    ///
    /// Requires eval_tree to have been built over exactly n roots where
    /// n is a power of 2, using product_tree_from_roots.
    pub fn multieval_from_balanced_tree(&self, eval_tree: &EvalTree<Fp>) -> Vec<Fp> {
        let n = eval_tree.n_roots;
        debug_assert!(n.is_power_of_two());
        let log_n = n.trailing_zeros() as usize;

        let degp_full = match self.degree() {
            Some(d) => d,
            None => return vec![Fp::ZERO; n],
        };

        // Reduce P mod Q if necessary.
        let n_layers = eval_tree.tree.len();
        let owned: Option<Vec<Fp>> = if degp_full >= n {
            let top = &eval_tree.tree[n_layers - 1];
            let mut rem = self.coeffs.clone();
            rem.resize(degp_full + 1, Fp::ZERO);
            for i in (n..=degp_full).rev() {
                if rem[i].is_zero() != u32::MAX {
                    let c = rem[i];
                    for j in 0..n {
                        rem[i - n + j] -= c * top[j];
                    }
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
            return vec![Fp::ZERO; n];
        }

        // Build rev_p zero-padded to length n.
        let mut rev_p = vec![Fp::ZERO; n];
        for i in 0..=degp {
            rev_p[i] = p_coeffs[degp - i];
        }

        // Scaled series: rev(P) * rev(Q)^{-1} mod x^n.
        let mut series = vec![Fp::ZERO; n];
        Self::mul_low(&mut series, &rev_p, &eval_tree.rev_q_inv, n);

        // Reverse back to polynomial order into the working buffer.
        // dst_buf[0..n] holds the root-level node.
        let mut dst_buf = vec![Fp::ZERO; n];
        for i in 0..=degp {
            dst_buf[i] = series[degp - i];
        }

        // Scratch buffer: largest mul_middle call is at the first descent step
        // with child_size = n/2.
        let max_scratch = Self::mul_middle_scratch_size(n / 2);
        let mut scratch = vec![Fp::ZERO; max_scratch];

        // dst buffer reused each layer.
        let mut dst = vec![Fp::ZERO; n];

        // Descend from the root toward the leaves.
        // After processing layer lyr (0-indexed from leaves), we split
        // nodes of size (1 << (lyr+1)) into pairs of size (1 << lyr).
        //
        // Loop variable: node_size is the size of the CURRENT (parent) nodes.
        // child_size = node_size / 2.
        // n_nodes = n / node_size  (number of parent nodes at this step).
        //
        // We go from node_size=n (1 parent) down to node_size=2 (n/2 parents,
        // each producing 2 children of size 1).
        let mut node_size = n;

        for lyr in (0..log_n).rev() {
            // lyr is the child layer index (0 = leaves).
            // child_size = 1 << lyr.
            // parent layer index = lyr + 1, node_size = 1 << (lyr + 1).
            let child_size = 1_usize << lyr;
            debug_assert_eq!(node_size, child_size * 2);

            let n_parents = n / node_size;

            // At child layer lyr, node j has non-leading coefficients at:
            // eval_tree.tree[lyr][j * child_size .. (j+1) * child_size]
            let lyr_flat = &eval_tree.tree[lyr];

            for j in 0..n_parents {
                // Parent node j occupies dst_buf[j*node_size .. (j+1)*node_size].
                let p_start = j * node_size;
                let parent = &dst_buf[p_start..p_start + node_size];

                // parent_tail = parent[1..], length = node_size - 1 = 2*child_size - 1.
                // This is the `f` argument to mul_middle (length 2*child_size - 1).
                let parent_tail = &parent[1..]; // length 2*child_size - 1

                // Left child (index 2j) non-leading coeffs.
                let q_l = &lyr_flat[2 * j * child_size..(2 * j + 1) * child_size];
                // Right child (index 2j+1) non-leading coeffs.
                let q_r = &lyr_flat[(2 * j + 1) * child_size..(2 * j + 2) * child_size];

                // Both q_l and q_r have length child_size.
                // parent_tail has length 2*child_size - 1.
                // So mul_middle(out, parent_tail, q_x) is exactly the balanced case.
                debug_assert_eq!(parent_tail.len(), 2 * child_size - 1);
                debug_assert_eq!(q_l.len(), child_size);
                debug_assert_eq!(q_r.len(), child_size);

                // Left child remainder = MP(parent_tail, q_r) + parent[..child_size]
                {
                    let out = &mut dst[2 * j * child_size..(2 * j + 1) * child_size];
                    Self::mul_middle(out, parent_tail, q_r, &mut scratch);
                    Self::add_into(out, &parent[..child_size]);
                }

                // Right child remainder = MP(parent_tail, q_l) + parent[..child_size]
                {
                    let out = &mut dst[(2 * j + 1) * child_size..(2 * j + 2) * child_size];
                    Self::mul_middle(out, parent_tail, q_l, &mut scratch);
                    Self::add_into(out, &parent[..child_size]);
                }
            }

            std::mem::swap(&mut dst_buf, &mut dst);
            node_size = child_size;
        }

        // dst_buf[i] now holds P(roots[i]) for i in 0..n.
        dst_buf
    }

    // Multi point evaluation methods ===========================================

    /// Evaluate self at each point in ai using a scaled remainder tree.
    ///
    /// Builds an EvalTree internally. If the same set of evaluation points will
    /// be used more than once, build an EvalTree explicitly and call
    /// multieval_from_tree to avoid recomputing the product tree each time.
    pub fn multieval(&self, ai: &[Fp]) -> Vec<Fp> {
        let eval_tree = EvalTree::new(ai);
        self.multieval_from_tree(&eval_tree)
    }

    /// Evaluate self at each point in roots using n calls to Horner's method.
    pub fn multieval_horner(&self, roots: &[Fp]) -> Vec<Fp> {
        let mut res = vec![Fp::ZERO; roots.len()];
        for (i, r) in roots.iter().enumerate() {
            res[i] = self.evaluate(&r)
        }
        res
    }

    /// Compute the product of self evaluated at each root in eval_tree, using a
    /// pre-built EvalTree. Equivalent to the resultant of self against the
    /// polynomial whose roots are those in the tree.
    pub fn resultant_from_roots_with_tree(&self, eval_tree: &EvalTree<Fp>) -> Fp {
        let mut res = Fp::ONE;
        for a in self.multieval_from_tree(eval_tree).iter() {
            res *= *a;
        }
        res
    }

    /// Compute the product of self evaluated at each root using a scaled
    /// remainder tree. Builds an EvalTree internally.
    pub fn resultant_from_roots(&self, roots: &[Fp]) -> Fp {
        let eval_tree = EvalTree::new(roots);
        self.resultant_from_roots_with_tree(&eval_tree)
    }

    /// Compute the product of self evaluated at each root using Horner's method.
    pub fn resultant_from_roots_horner(&self, ai: &[Fp]) -> Fp {
        let mut res = Fp::ONE;
        for a in ai {
            res *= self.evaluate(a)
        }
        res
    }

    /// Experimental strategy where we take 2^n elements from roots for a balanced product
    /// tree, and then the remaining elements are multiplied in with Horner
    pub fn resultant_mixed_strategy_with_tree(
        &self,
        balanced_tree: &EvalTree<Fp>,
        rem_roots: &[Fp],
    ) -> Fp {
        // Compute the resultant of 2^n elements with a balanced tree
        let mut resultant = Fp::ONE;
        for r in self.multieval_from_balanced_tree(&balanced_tree) {
            resultant *= r
        }

        // Use Horner for the rest, if I got things faster we could use an unbalanced tree?
        for a in rem_roots {
            resultant *= self.evaluate(a)
        }

        resultant
    }

    /// Experimental strategy where we take 2^n elements from roots for a balanced product
    /// tree, and then the remaining elements are multiplied in with Horner
    pub fn resultant_mixed_strategy(&self, ai: &[Fp]) -> Fp {
        // For the first 2^n elements, use a product tree
        let n_roots = ai.len();
        let n = usize::BITS - n_roots.leading_zeros() - 1;
        let two_n = 1 << n;
        let balanced_tree = EvalTree::new(&ai[..two_n]);

        // For the rest of it, use Horner, but we could make two trees, one balanced and one unbalanced?
        let mut resultant = Fp::ONE;
        for r in self.multieval_from_balanced_tree(&balanced_tree) {
            resultant *= r
        }

        for a in &ai[two_n..] {
            resultant *= self.evaluate(a)
        }

        resultant
    }

    // Extra functions for benchmarking and testing ===========================================

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
        let n = other.len();
        let mut coeffs = vec![Fp::ZERO; n];
        let mut scratch = vec![Fp::ZERO; Self::mul_middle_scratch_size(n)];
        Self::mul_middle(&mut coeffs, &self.coeffs, &other.coeffs, &mut scratch);
        Self { coeffs }
    }

    /// TODO: testing only.
    pub fn karatsuba_mul(&self, other: &Self) -> Self {
        let mut coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];

        let n_f = self.len();
        let n_g = other.len();

        let scratch_len = Self::karatsuba_scratch_size(n_f.max(n_g), n_f.min(n_g));
        let mut tmp = vec![Fp::ZERO; scratch_len];
        Self::karatsuba_multiplication(&mut coeffs, &self.coeffs, &other.coeffs, &mut tmp);
        Self { coeffs }
    }
}

// Trait Methods implementations ===========================================================

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

// End of Trait Methods ===========================================================

#[cfg(test)]
mod test_poly {
    use crate::fields::sqisign::SqiField248Base as Fp;
    use crate::{
        polynomial_ring::poly::{EvalTree, Polynomial},
        utilities::test_utils::drng::DRNG,
    };

    type PR = Polynomial<Fp>;

    #[test]
    fn test_addition() {
        let mut rng = DRNG::from_seed("polynomial_addition".as_bytes());

        let f = PR::rand(&mut rng, 4);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 6);

        let t1 = &(&f + &g) + &h;
        let t2 = &f + &(&g + &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 += &g;
        t1 += &h;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = PR::new_from_ele(&Fp::ZERO);
        let t1 = &f + &zero;
        assert!(t1.equals(&f) == u32::MAX);

        let f_neg = -&f;
        assert!((&f + &f_neg).is_zero() == u32::MAX);
    }

    #[test]
    fn test_subtraction() {
        let mut rng = DRNG::from_seed("polynomial_subtraction".as_bytes());

        let f = PR::rand(&mut rng, 4);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 6);

        let t1 = &(&f - &g) - &h;
        let t2 = &f - &(&g + &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 -= &g;
        t1 -= &h;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = &f - &f;
        assert!(zero.is_zero() == u32::MAX);
    }

    #[test]
    fn test_multiplication() {
        let mut rng = DRNG::from_seed("polynomial_multiplication".as_bytes());

        let f = PR::rand(&mut rng, 4);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 6);

        let t1 = &(&f * &g) * &h;
        let t2 = &f * &(&g * &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 *= &g;
        t1 *= &h;
        assert!(t1.equals(&t2) == u32::MAX);

        let one = PR::new_from_ele(&Fp::ONE);
        t1 *= &one;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = PR::new_from_ele(&Fp::ZERO);
        t1 *= &zero;
        assert!(t1.is_zero() == u32::MAX);
    }

    #[test]
    fn test_multiplication_large() {
        let mut rng = DRNG::from_seed("polynomial_multiplication_large".as_bytes());

        for i in 0..5 {
            let f = PR::rand(&mut rng, 100 - i);
            let g = PR::rand(&mut rng, 150 - 20 * i);
            let h = PR::rand(&mut rng, 200 - 15 * i);

            let t1 = f.basic_mul(&g).basic_mul(&h);
            let t2 = f.karatsuba_mul(&g).karatsuba_mul(&h);

            assert!(t1.equals(&t2) == u32::MAX);
        }
    }

    #[test]
    fn test_squaring() {
        // TODO: do specalised squaring method?
        let mut rng = DRNG::from_seed("polynomial_squaring".as_bytes());

        let f = PR::rand(&mut rng, 5);

        let t1 = &f * &f;
        let mut t2 = f.clone();
        t2 *= &f;

        assert!(t1.equals(&t2) == u32::MAX);
    }

    #[test]
    fn test_scaling() {
        let mut rng = DRNG::from_seed("polynomial_scaling".as_bytes());

        let f = PR::rand(&mut rng, 5);

        // Compute [4]f in different ways
        let mut t1 = f.clone();
        t1 += &f;
        t1 += &f;
        t1 += &f;

        let t2 = f.scale(&Fp::FOUR);
        assert!(t1.equals(&t2) == u32::MAX);

        let mut t2 = f.clone();
        t2 *= &Fp::FOUR;
        assert!(t1.equals(&t2) == u32::MAX);

        let t2 = f.scale_small(4);
        assert!(t1.equals(&t2) == u32::MAX);

        let g = f.scale(&Fp::ONE);
        assert!(f.equals(&g) == u32::MAX);

        let g = f.scale_small(1);
        assert!(f.equals(&g) == u32::MAX);

        let g = f.scale(&Fp::ZERO);
        assert!(g.is_zero() == u32::MAX);

        let g = f.scale_small(0);
        assert!(g.is_zero() == u32::MAX);
    }

    #[test]
    fn test_evaluation() {
        let mut rng = DRNG::from_seed("polynomial_eval".as_bytes());

        // Evaluating a constant polynomial should always give the same result
        let a = Fp::rand(&mut rng);
        let b = Fp::rand(&mut rng);
        let f = PR::rand(&mut rng, 1);
        let ea = f.evaluate(&a);
        let eb = f.evaluate(&b);
        assert!(ea.equals(&eb) == u32::MAX);

        // Random degree four polynomial
        let f = PR::rand(&mut rng, 5);

        // Evaluating at zero should give you the constant coefficient
        let e0 = f.evaluate(&Fp::ZERO);
        let c0 = f.constant_coefficient().unwrap();
        assert!(e0.equals(&c0) == u32::MAX);

        // Evaluating at one should give you the sum of the coefficients
        let e1 = f.evaluate(&Fp::ONE);
        let c1 = f[0] + f[1] + f[2] + f[3] + f[4];
        assert!(e1.equals(&c1) == u32::MAX);

        // Dumb evaluation to compare against the type method.
        let a = Fp::rand(&mut rng);
        let ea = f.evaluate(&a);
        let mut test_ea = f[0];
        test_ea += f[1] * a;
        test_ea += f[2] * (a * a);
        test_ea += f[3] * (a * a * a);
        test_ea += f[4] * (a * a * a * a);

        assert!(ea.equals(&test_ea) == u32::MAX);
    }

    #[test]
    fn test_degree() {
        let f = PR::new_from_slice(&[Fp::ZERO]);
        assert!(f.degree() == None);

        let f = PR::new_from_slice(&[Fp::ONE]);
        assert!(f.degree().unwrap() == 0);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ONE]);
        assert!(f.degree().unwrap() == 1);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ONE, Fp::ZERO]);
        assert!(f.degree().unwrap() == 4);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO]);
        assert!(f.degree() == None);

        let f = PR::new_from_slice(&[Fp::ONE, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO]);
        assert!(f.degree().unwrap() == 0);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ONE, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO]);
        assert!(f.degree().unwrap() == 1);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ONE, Fp::ZERO]);
        assert!(f.degree().unwrap() == 4);
    }

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
                let mut scratch = vec![Fp::ZERO; Polynomial::<Fp>::mul_middle_scratch_size(n)];
                Polynomial::mul_middle(&mut h, &f.coeffs, &g.coeffs, &mut scratch);

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
            let mut scratch = vec![Fp::ZERO; Polynomial::<Fp>::mul_middle_scratch_size(n)];
            Polynomial::mul_middle(&mut h, &f.coeffs, &g.coeffs, &mut scratch);
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
                let mut scratch = vec![Fp::ZERO; Polynomial::<Fp>::mul_middle_scratch_size(n)];
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
                let scratch_size = Polynomial::<Fp>::mul_middle_scratch_size(in_n.max(out_n));
                let mut scratch = vec![Fp::ZERO; scratch_size];
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

            let old_tree = PR::product_tree_from_roots_simple(&roots);
            let old_root = old_tree.last().unwrap()[0].clone();

            let (flat_tree, offsets) = PR::product_tree_from_roots(&roots);
            let n_layers = flat_tree.len();

            assert_eq!(flat_tree[n_layers - 1].len(), n, "root layer length n={n}");
            assert_eq!(offsets[n_layers - 1], vec![0, n], "root offsets n={n}");

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
        // f(0) = f[0], f(1) = sum of all coefficients -- easy to verify.
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

    #[test]
    fn test_resultant_mixed_strategy_matches_horner() {
        let mut rng = DRNG::from_seed("resultant".as_bytes());
        for _ in 0..10 {
            for &(deg, n_roots) in &[(5usize, 5usize), (10, 10), (20, 20), (50, 50)] {
                let f = PR::rand(&mut rng, deg + 1);
                let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
                let res = f.resultant_mixed_strategy(&roots);
                let horner = f.resultant_from_roots_horner(&roots);
                assert_eq!(res.equals(&horner), u32::MAX);
            }
        }
    }
}
