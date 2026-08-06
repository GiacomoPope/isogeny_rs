use core::ops::{Add, Div, Mul, Neg, Sub};
use core::fmt::Debug;
use core::marker::PhantomData;
use crate::bigint::BigIntAlg;
use crate::quaternion::hnf_ext::HnfModExt;


// ========================================================================
// Algebra Configuration (Zero-Sized Types)
// ========================================================================

/// Trait to securely bind a static prime `p` to a specific algebra type
pub trait QuatConfig<T: BigIntAlg> {
    fn p() -> &'static T;
}

// ========================================================================
// Integer Quaternions (Numerator Only)
// ========================================================================

pub struct IntQuat<T: BigIntAlg, P: QuatConfig<T>> {
    pub coords: [T; 4],
    _marker: PhantomData<P>,
}

impl<T: BigIntAlg, P: QuatConfig<T>> Clone for IntQuat<T, P> {
    fn clone(&self) -> Self {
        Self {
            coords: self.coords.clone(),
            _marker: PhantomData,
        }
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> Debug for IntQuat<T, P> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("IntQuat").field("coords", &self.coords).finish()
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> IntQuat<T, P> {
    pub fn new(c0: T, c1: T, c2: T, c3: T) -> Self {
        Self {
            coords: [c0, c1, c2, c3],
            _marker: PhantomData,
        }
    }

    pub fn new_i32(c0: i32, c1: i32, c2: i32, c3: i32) -> Self {
        Self::new(T::from_i32(c0), T::from_i32(c1), T::from_i32(c2), T::from_i32(c3))
    }

    pub fn zero() -> Self {
        Self::new_i32(0, 0, 0, 0)
    }

    pub fn ct_eq(&self, other: &Self) -> u32 {
        self.coords[0].ct_eq(&other.coords[0])
            & self.coords[1].ct_eq(&other.coords[1])
            & self.coords[2].ct_eq(&other.coords[2])
            & self.coords[3].ct_eq(&other.coords[3])
    }

    pub fn is_zero(&self) -> bool {
        // TODO: make that constant time
        self.coords[0].is_zero()
            && self.coords[1].is_zero()
            && self.coords[2].is_zero()
            && self.coords[3].is_zero()
    }

    pub fn conj(&self) -> Self {
        Self {
            coords: [
                self.coords[0].clone(),
                -self.coords[1].clone(),
                -self.coords[2].clone(),
                -self.coords[3].clone(),
            ],
            _marker: PhantomData,
        }
    }

    pub fn norm(&self) -> T {
        let conj_x = self.conj();
        // Multiplication implicitly uses P::p() via the overloaded operator
        // TODO: faster, direct computation
        let norm_elem = self * &conj_x;
        norm_elem.coords[0].clone().abs()
    }


    // HNF additions
    /// Computes the GCD of the 4 coordinates of this quaternion
    pub fn coords_gcd(&self) -> T {
        let mut g = self.coords[0].clone();
        for i in 1..4 {
            g = g.gcd(&self.coords[i]);
        }
        g
    }

    /// Modifies the quaternion in place by applying centered_mod to each coordinate
    pub fn centered_mod_mut(&mut self, m: &T) {
        for i in 0..4 {
            self.coords[i] = self.coords[i].centered_mod(m);
        }
    }
}

// --- IntQuat Standard Operators ---

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Add<&'b IntQuat<T, P>> for &'a IntQuat<T, P> {
    type Output = IntQuat<T, P>;
    fn add(self, rhs: &'b IntQuat<T, P>) -> Self::Output {
        IntQuat {
            coords: [
                self.coords[0].clone() + rhs.coords[0].clone(),
                self.coords[1].clone() + rhs.coords[1].clone(),
                self.coords[2].clone() + rhs.coords[2].clone(),
                self.coords[3].clone() + rhs.coords[3].clone(),
            ],
            _marker: PhantomData,
        }
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Sub<&'b IntQuat<T, P>> for &'a IntQuat<T, P> {
    type Output = IntQuat<T, P>;
    fn sub(self, rhs: &'b IntQuat<T, P>) -> Self::Output {
        IntQuat {
            coords: [
                self.coords[0].clone() - rhs.coords[0].clone(),
                self.coords[1].clone() - rhs.coords[1].clone(),
                self.coords[2].clone() - rhs.coords[2].clone(),
                self.coords[3].clone() - rhs.coords[3].clone(),
            ],
            _marker: PhantomData,
        }
    }
}

// Scalar Multiplication: &IntQuat * &T
impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Mul<&'b T> for &'a IntQuat<T, P> {
    type Output = IntQuat<T, P>;

    fn mul(self, scalar: &'b T) -> Self::Output {
        IntQuat {
            coords: [
                self.coords[0].clone() * scalar.clone(),
                self.coords[1].clone() * scalar.clone(),
                self.coords[2].clone() * scalar.clone(),
                self.coords[3].clone() * scalar.clone(),
            ],
            _marker: PhantomData,
        }
    }
}

// Quaternion Multiplication: &IntQuat * &IntQuat
impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Mul<&'b IntQuat<T, P>> for &'a IntQuat<T, P> {
    type Output = IntQuat<T, P>;
    fn mul(self, rhs: &'b IntQuat<T, P>) -> Self::Output {
        let p = P::p();
        let a = &self.coords;
        let b = &rhs.coords;

        // Quaternion multiplication in 12 multiplications
        // Using Gauss' optimization:
        //    consider a = A + Bj    (A = a₀+a₁i, B = a₂+a₃i)
        //             b = C + Dj    (C = b₀+b₁i, D = b₂+b₃i)
        //    we then have:
        //       (A+Bj)(C+Dj) = (AC - pBD̅) + (AD + BC̅)j
        //    we can use Karatsuba to perform each complex multiplication
        //    in 3 element multiplication instead of 4.

        let m1 = a[0].clone() * b[0].clone();
        let m2 = a[1].clone() * b[1].clone();
        let m3 = (a[0].clone() + a[1].clone()) * (b[0].clone() + b[1].clone());

        let m4 = a[2].clone() * b[2].clone();
        let m5 = a[3].clone() * b[3].clone();
        let m6 = (a[2].clone() + a[3].clone()) * (b[2].clone() - b[3].clone());

        let m7 = a[0].clone() * b[2].clone();
        let m8 = a[1].clone() * b[3].clone();
        let m9 = (a[0].clone() + a[1].clone()) * (b[2].clone() + b[3].clone());

        let m10 = a[2].clone() * b[0].clone();
        let m11 = a[3].clone() * b[1].clone();
        let m12 = (a[2].clone() + a[3].clone()) * (b[0].clone() - b[1].clone());

        let c0 = (m1.clone() - m2.clone()) - p.clone() * (m4.clone() + m5.clone());
        let c1 = (m3 - m1.clone() - m2.clone()) - p.clone() * (m6 - m4 + m5);
        let c2 = (m7.clone() - m8.clone()) + (m10.clone() + m11.clone());
        let c3 = (m9 - m7 - m8) + (m12 - m10 + m11);

        IntQuat {
            coords: [c0, c1, c2, c3],
            _marker: PhantomData,
        }
    }
}

impl<'a, T: BigIntAlg, P: QuatConfig<T>> Neg for &'a IntQuat<T, P> {
    type Output = IntQuat<T, P>;
    fn neg(self) -> Self::Output {
        IntQuat {
            coords: [
                -self.coords[0].clone(),
                -self.coords[1].clone(),
                -self.coords[2].clone(),
                -self.coords[3].clone(),
            ],
            _marker: PhantomData,
        }
    }
}

// Scalar Division: &IntQuat / &T (Assumes exact division)
impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Div<&'b T> for &'a IntQuat<T, P> {
    type Output = IntQuat<T, P>;

    fn div(self, scalar: &'b T) -> Self::Output {
        IntQuat {
            coords: [
                self.coords[0].clone() / scalar.clone(),
                self.coords[1].clone() / scalar.clone(),
                self.coords[2].clone() / scalar.clone(),
                self.coords[3].clone() / scalar.clone(),
            ],
            _marker: PhantomData,
        }
    }
}

// --- Macros to auto-implement Val/Ref combinations for +, -, * ---
macro_rules! impl_binop {
    ($trait:ident, $method:ident) => {
        impl<T: BigIntAlg, P: QuatConfig<T>> $trait<IntQuat<T, P>> for IntQuat<T, P> {
            type Output = IntQuat<T, P>;
            fn $method(self, rhs: IntQuat<T, P>) -> Self::Output { (&self).$method(&rhs) }
        }
        impl<'a, T: BigIntAlg, P: QuatConfig<T>> $trait<&'a IntQuat<T, P>> for IntQuat<T, P> {
            type Output = IntQuat<T, P>;
            fn $method(self, rhs: &'a IntQuat<T, P>) -> Self::Output { (&self).$method(rhs) }
        }
        impl<'a, T: BigIntAlg, P: QuatConfig<T>> $trait<IntQuat<T, P>> for &'a IntQuat<T, P> {
            type Output = IntQuat<T, P>;
            fn $method(self, rhs: IntQuat<T, P>) -> Self::Output { self.$method(&rhs) }
        }
    };
}
impl_binop!(Add, add);
impl_binop!(Sub, sub);
impl_binop!(Mul, mul);

macro_rules! impl_scalar_mul {
    () => {
        impl<T: BigIntAlg, P: QuatConfig<T>> Mul<T> for IntQuat<T, P> {
            type Output = IntQuat<T, P>;
            fn mul(self, rhs: T) -> Self::Output { (&self).mul(&rhs) }
        }
        impl<'a, T: BigIntAlg, P: QuatConfig<T>> Mul<&'a T> for IntQuat<T, P> {
            type Output = IntQuat<T, P>;
            fn mul(self, rhs: &'a T) -> Self::Output { (&self).mul(rhs) }
        }
        impl<'a, T: BigIntAlg, P: QuatConfig<T>> Mul<T> for &'a IntQuat<T, P> {
            type Output = IntQuat<T, P>;
            fn mul(self, rhs: T) -> Self::Output { self.mul(&rhs) }
        }
    }
}
impl_scalar_mul!();

// ========================================================================
// Rational Quaternions (Integer + Denominator)
// ========================================================================

#[derive(Debug)]
pub struct RatQuat<T: BigIntAlg, P: QuatConfig<T>> {
    pub num: IntQuat<T, P>,
    pub denom: T,
}

impl<T: BigIntAlg, P: QuatConfig<T>> Clone for RatQuat<T, P> {
    fn clone(&self) -> Self {
        Self {
            num: self.num.clone(),
            denom: self.denom.clone(),
        }
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> RatQuat<T, P> {
    pub fn new(num: IntQuat<T, P>, denom: T) -> Self {
        debug_assert!(!denom.is_zero(), "RatQuat initialized with zero denominator");
        Self { num, denom }
    }

    pub fn new_i32(denom: i32, c0: i32, c1: i32, c2: i32, c3: i32) -> Self {
        debug_assert!(denom != 0, "RatQuat initialized with zero denominator");
        Self {
            num: IntQuat::new_i32(c0, c1, c2, c3),
            denom: T::from_i32(denom),
        }
    }

    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }

    pub fn ct_eq(&self, other: &Self) -> u32 {
        self.num.ct_eq(&other.num) & self.denom.ct_eq(&other.denom)
    }

    pub fn conj(&self) -> Self {
        Self {
            num: self.num.conj(),
            denom: self.denom.clone(),
        }
    }

    pub fn norm(&self) -> (T, T) {
        let num_norm = self.num.norm();
        let denom_sq = self.denom.clone() * self.denom.clone();

        let g = num_norm.gcd(&denom_sq);
        let res_num = (num_norm / g.clone()).abs();
        let res_denom = (denom_sq / g).abs();

        (res_num, res_denom)
    }

    pub fn normalize(&mut self) {
        debug_assert!(!self.denom.is_zero(), "RatQuat::normalize called on zero denominator");

        let mut g = self.num.coords[0].clone();
        for i in 1..4 {
            g = g.gcd(&self.num.coords[i]);
        }
        g = g.gcd(&self.denom);

        let mut sign = T::one();
        if self.denom < T::zero() {
            sign = T::from_i32(-1);
        }

        let divisor = g * sign;
        debug_assert!(!divisor.is_zero(), "RatQuat::normalize encountered zero divisor during reduction");

        for i in 0..4 {
            self.num.coords[i] = self.num.coords[i].clone() / divisor.clone();
        }
        self.denom = self.denom.clone() / divisor;
    }

    // --- Unnormalized (Lazy) Operations ---

    pub fn add_lazy(&self, other: &Self) -> Self {
        let gcd = self.denom.gcd(&other.denom);
        let a_mult = other.denom.clone() / gcd.clone();
        let b_mult = self.denom.clone() / gcd;

        let a_scaled = &self.num * &a_mult;
        let b_scaled = &other.num * &b_mult;

        Self {
            num: &a_scaled + &b_scaled,
            denom: self.denom.clone() * a_mult,
        }
    }

    pub fn sub_lazy(&self, other: &Self) -> Self {
        let gcd = self.denom.gcd(&other.denom);
        let a_mult = other.denom.clone() / gcd.clone();
        let b_mult = self.denom.clone() / gcd;

        let a_scaled = &self.num * &a_mult;
        let b_scaled = &other.num * &b_mult;

        Self {
            num: &a_scaled - &b_scaled,
            denom: self.denom.clone() * a_mult,
        }
    }

    pub fn mul_lazy(&self, other: &Self) -> Self {
        Self {
            num: &self.num * &other.num,
            denom: self.denom.clone() * other.denom.clone(),
        }
    }
}

// --- RatQuat Standard Operators (Implicitly Normalized) ---

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Add<&'b RatQuat<T, P>> for &'a RatQuat<T, P> {
    type Output = RatQuat<T, P>;
    fn add(self, rhs: &'b RatQuat<T, P>) -> Self::Output {
        // TODO: is it possible to perform lcm without division for marginal
        // speed gain?
        let mut res = self.add_lazy(rhs);
        res.normalize();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Sub<&'b RatQuat<T, P>> for &'a RatQuat<T, P> {
    type Output = RatQuat<T, P>;
    fn sub(self, rhs: &'b RatQuat<T, P>) -> Self::Output {
        let mut res = self.sub_lazy(rhs);
        res.normalize();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Mul<&'b RatQuat<T, P>> for &'a RatQuat<T, P> {
    type Output = RatQuat<T, P>;
    fn mul(self, rhs: &'b RatQuat<T, P>) -> Self::Output {
        let mut res = self.mul_lazy(rhs);
        res.normalize();
        res
    }
}
