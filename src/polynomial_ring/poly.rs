#![allow(dead_code)] // for now

use fp2::traits::Fp as FpTrait;
use rand_core::{CryptoRng, RngCore};

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

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

    fn product_tree_root(v: &[Self]) -> Self;
    fn product_tree_quadratic_leaves(leaves: &[[Fp; 3]]) -> Vec<Vec<Fp>>;
    fn product_from_tree(tree: &[Vec<Fp>]) -> Self;

    fn resultant_from_roots(&self, ai: &[Fp]) -> Fp;
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

    /// Truncate leading zeros from the polynomial
    fn truncate_leading_zeros(&mut self) {
        // find the index of the first non-zero element
        let mut i = self.len() - 1;
        while i > 0 && self.coeffs[i].is_zero() == u32::MAX {
            i -= 1;
        }
        self.coeffs.truncate(i);
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

        // For small degree f, g we use basic multiplication strategies with
        // O(n^2) operations. TODO: set the right bound for when to fall back
        // to this.
        if f.len() <= 4 {
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
        assert!(f.len() + g.len() - 1 <= fg.len());
        Self::karatsuba_multiplication(fg, f, g)
    }

    /// TODO: testing only.
    pub fn basic_mul(&self, other: &Self) -> Self {
        let mut coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];
        Self::schoolbook_multiplication(&mut coeffs, &self.coeffs, &other.coeffs);
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

    /// Computes the root of a product tree given a slice of the leaves.
    pub fn product_tree_root(v: &[Self]) -> Self {
        if v.len() == 1 {
            return v[0].clone();
        }
        let half = v.len() >> 1;
        &Self::product_tree_root(&v[..half]) * &Self::product_tree_root(&v[half..])
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
        // 4 (and potenitally one of degree 2 if n is odd) and so on until the
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
    pub fn product_from_tree(tree: &[Vec<Fp>]) -> Self {
        let n = tree[0].len() / 3;
        let root = tree.last().unwrap();
        Self {
            coeffs: root[..2 * n + 1].to_vec(),
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

        // Iterate from the last to the first coefficicent using Horner's method
        // to evaluate the polynomial.
        let mut bi = Fp::ZERO;
        let deg = self.degree().unwrap();
        for i in 0..=deg {
            bi = *a * bi + self.coeffs[deg - i];
        }
        bi
    }

    /// Compute the resultant of self with a polynomial g = \prod {x - ai}
    /// given the roots ai.
    // TODO: this is a very slow and stupid method, but speed comes later and
    // I want to sketch sqrt velu.
    pub fn resultant_from_roots(&self, ai: &[Fp]) -> Fp {
        let mut res = Fp::ONE;
        for a in ai.iter() {
            res *= self.evaluate(a);
        }
        res
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

    fn product_tree_root(v: &[Self]) -> Self {
        Self::product_tree_root(v)
    }

    fn product_tree_quadratic_leaves(leaves: &[[Fp; 3]]) -> Vec<Vec<Fp>> {
        Self::product_tree_quadratic_leaves(leaves)
    }

    fn product_from_tree(tree: &[Vec<Fp>]) -> Self {
        Self::product_from_tree(tree)
    }

    fn resultant_from_roots(&self, ai: &[Fp]) -> Fp {
        self.resultant_from_roots(ai)
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
                write!(f, "({})", c)?
            } else {
                write!(f, "({})*x^{} + ", c, i)?
            }
        }
        Ok(())
    }
}
