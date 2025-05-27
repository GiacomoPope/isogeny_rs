#![allow(dead_code)] // for now

use fp2::traits::Fp as FpTrait;

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

    /// Reverse the coefficients of self in place.
    fn reverse_into(&mut self) {
        self.coeffs.reverse();
    }

    /// Return the polynomial with coefficents reversed.
    pub fn reverse(self) -> Polynomial<Fp> {
        let mut r = self;
        r.reverse_into();
        r
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
    /// the result (len(f) + len(g) - 2).
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
        let mut fg_coeffs = vec![Fp::ZERO; self.len() + other.len() - 2];
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
    pub fn scale(self, c: &Fp) -> Polynomial<Fp> {
        let mut r = self;
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
    pub fn scale_small(self, k: i32) -> Polynomial<Fp> {
        let mut r = self;
        r.scale_small_into(k);
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

impl<Fp: FpTrait> AddAssign for Polynomial<Fp> {
    #[inline(always)]
    fn add_assign(&mut self, other: Polynomial<Fp>) {
        self.set_add(&other);
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

impl<Fp: FpTrait> SubAssign for Polynomial<Fp> {
    #[inline(always)]
    fn sub_assign(&mut self, other: Polynomial<Fp>) {
        self.set_sub(&other);
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

impl<Fp: FpTrait> MulAssign for Polynomial<Fp> {
    #[inline(always)]
    fn mul_assign(&mut self, other: Polynomial<Fp>) {
        self.set_mul(&other);
    }
}

impl<Fq: FpTrait> ::std::fmt::Display for Polynomial<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "TODO")
    }
}
