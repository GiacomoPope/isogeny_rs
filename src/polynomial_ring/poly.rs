#![allow(dead_code)] // for now

use fp2::traits::Fp as FpTrait;
use rand_core::{CryptoRng, RngCore};

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

/// Trait for arithmetic for univariate polynomials in Fp[X]
pub trait Poly:
    Index<usize>
    + IndexMut<usize>
    + Sized
    + Neg<Output = Self>
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Display
{
    // TODO
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

    /// The length of the polynomial. TODO: should we trim trailing zeros? If so, how often?
    fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// The degree of the polynomial. TODO: should we trim trailing zeros? If so, how often?
    /// We return None for the zero polynomial (instead of -inf which is a bit too big for my
    /// computer...)
    fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.len() - 1)
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
    pub fn reverse(&self) -> Polynomial<Fp> {
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

    /// Set self to it's negative.
    fn set_neg(&mut self) {
        for x in self.coeffs.iter_mut() {
            x.set_neg();
        }
    }

    /// Set self <- self + other
    // TODO: should we assume other has length smaller or equal to
    // self, or dynamically extend?
    fn set_add(&mut self, other: &Self) {
        let k = self.len().min(other.len());
        for i in 0..k {
            self.coeffs[i] += other[i];
        }
    }

    /// Set self <- self - other
    // TODO: should we assume other has length smaller or equal to
    // self, or dynamically extend?
    fn set_sub(&mut self, other: &Self) {
        let k = self.len().min(other.len());
        for i in 0..k {
            self.coeffs[i] -= other[i];
        }
    }

    /// Compute f * g with O(len(f) * len(g)) Fq multiplications using
    /// schoolbook multiplication. Assumes that fg has enough space for
    /// the result (len(f) + len(g) - 1).
    fn schoolbook_multiplication(fg: &mut [Fp], f: &[Fp], g: &[Fp]) {
        for i in 0..f.len() {
            for j in 0..g.len() {
                // TODO: this could be sped up when we know c[i + j] is zero
                // which happens when i = 0 or when j + 1 = len(other)
                fg[i + j] += f[i] * g[j]
            }
        }
    }

    /// Set self <- self * other
    fn set_mul(&mut self, other: &Self) {
        let mut fg_coeffs = vec![Fp::ZERO; self.len() + other.len() - 1];
        Self::schoolbook_multiplication(&mut fg_coeffs, &self.coeffs, &other.coeffs);
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

    /// Evaluate a polynomial at a value `a`
    pub fn evaluate(&self, a: &Fp) -> Fp {
        // Handle degree 0 and 1 cases early.
        if self.len() == 0 {
            return Fp::ZERO;
        } else if self.len() == 1 {
            return self.coeffs[0];
        } else if self.len() == 2 {
            return self.coeffs[0] + *a * self.coeffs[1];
        }

        // Otherwise compute (c0 + c1 * a) for the linear piece
        let mut res = self.coeffs[0] + *a * self.coeffs[1];
        // Precompute a^2 for the quadratic piece
        let mut a_n = a.square();
        for i in 2..self.len() {
            // Compute res += c_i * a^i for each step.
            res += self.coeffs[i] * a_n;

            // For the last step we can skip multiplying by a
            if i != self.len() - 1 {
                a_n *= *a;
            }
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

impl<Fp: FpTrait> Poly for Polynomial<Fp> {}

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

impl<Fp: FpTrait> Neg for Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn neg(self) -> Polynomial<Fp> {
        let mut r = self;
        r.set_neg();
        r
    }
}

impl<Fp: FpTrait> Add for Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn add(self, other: Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self;
        r.set_add(&other);
        r
    }
}

impl<Fp: FpTrait> Add for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn add(self, other: &Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_add(&other);
        r
    }
}

impl<Fp: FpTrait> AddAssign for Polynomial<Fp> {
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

impl<Fp: FpTrait> Sub for Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn sub(self, other: Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self;
        r.set_sub(&other);
        r
    }
}

impl<Fp: FpTrait> Sub for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn sub(self, other: &Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_sub(&other);
        r
    }
}

impl<Fp: FpTrait> SubAssign for Polynomial<Fp> {
    #[inline(always)]
    fn sub_assign(&mut self, other: Polynomial<Fp>) {
        self.set_sub(&other);
    }
}

impl<Fp: FpTrait> SubAssign<&Polynomial<Fp>> for Polynomial<Fp> {
    #[inline(always)]
    fn sub_assign(&mut self, other: &Polynomial<Fp>) {
        self.set_sub(other);
    }
}

impl<Fp: FpTrait> Mul for Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn mul(self, other: Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self;
        r.set_mul(&other);
        r
    }
}

impl<Fp: FpTrait> Mul for &Polynomial<Fp> {
    type Output = Polynomial<Fp>;

    #[inline(always)]
    fn mul(self, other: &Polynomial<Fp>) -> Polynomial<Fp> {
        let mut r = self.clone();
        r.set_mul(&other);
        r
    }
}

impl<Fp: FpTrait> MulAssign for Polynomial<Fp> {
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

impl<Fq: FpTrait> ::std::fmt::Display for Polynomial<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "TODO")
    }
}
