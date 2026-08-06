pub mod rat;

use core::fmt::Debug;
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};

// ========================================================================
// Generic BigInt Interface
// ========================================================================

pub trait BigIntAlg:
    Clone
    + Debug
    + PartialEq
    + Eq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Neg<Output = Self>
    + Sized
    + 'static
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_i32(val: i32) -> Self;
    fn gcd(&self, other: &Self) -> Self;
    fn xgcd(&self, other: &Self) -> (Self, Self, Self);
    fn abs(&self) -> Self;
    fn is_zero(&self) -> bool;
    fn sqrt(&self) -> Self;

    /// Returns u32::MAX (0xFFFFFFFF) if equal, 0 otherwise
    fn ct_eq(&self, other: &Self) -> u32;
}

// ========================================================================
// Backend Implementations
// ========================================================================

pub mod backends {
    use super::BigIntAlg;

    // --- num_bigint Backend ---
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{Signed, Zero, One};

    impl BigIntAlg for BigInt {
        fn zero() -> Self { <BigInt as Zero>::zero() }
        fn one() -> Self { <BigInt as One>::one() }
        fn from_i32(val: i32) -> Self { BigInt::from(val) }
        fn gcd(&self, other: &Self) -> Self { self.extended_gcd(other).gcd }
        fn xgcd(&self, other: &Self) -> (Self, Self, Self) {
            let res = self.extended_gcd(other);
            (res.gcd, res.x, res.y)
        }
        fn abs(&self) -> Self { <BigInt as Signed>::abs(self) }
        fn is_zero(&self) -> bool { <BigInt as Zero>::is_zero(self) }
        fn sqrt(&self) -> Self { self.sqrt() }
        fn ct_eq(&self, other: &Self) -> u32 {
            ((self == other) as u32).wrapping_neg()
        }
    }

    // --- crypto_bigint Backend ---
    use crypto_bigint::{Int, NonZero};
    use core::ops::{Add, Div, Mul, Neg, Rem, Sub};

    #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct CryptoInt<const LIMBS: usize>(pub Int<LIMBS>);

    impl<const LIMBS: usize> Add for CryptoInt<LIMBS> {
        type Output = Self;
        fn add(self, rhs: Self) -> Self::Output { CryptoInt(self.0 + rhs.0) }
    }

    impl<const LIMBS: usize> Sub for CryptoInt<LIMBS> {
        type Output = Self;
        fn sub(self, rhs: Self) -> Self::Output { CryptoInt(self.0 - rhs.0) }
    }

    impl<const LIMBS: usize> Mul for CryptoInt<LIMBS> {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self::Output { CryptoInt(self.0 * rhs.0) }
    }

    impl<const LIMBS: usize> Div for CryptoInt<LIMBS> {
        type Output = Self;
        fn div(self, rhs: Self) -> Self::Output {
            let nz: Option<NonZero<Int<LIMBS>>> = NonZero::new(rhs.0).into();
            let res: Option<Int<LIMBS>> = (self.0 / nz.expect("Division by zero in CryptoInt::div")).into();
            CryptoInt(res.expect("Division overflowed"))
        }
    }

    impl<const LIMBS: usize> Rem for CryptoInt<LIMBS> {
        type Output = Self;
        fn rem(self, rhs: Self) -> Self::Output {
            let nz: Option<NonZero<Int<LIMBS>>> = NonZero::new(rhs.0).into();
            let res: Option<Int<LIMBS>> = (self.0 % nz.expect("Division by zero in CryptoInt::rem")).into();
            CryptoInt(res.expect("Remainder overflowed"))
        }
    }

    impl<const LIMBS: usize> Neg for CryptoInt<LIMBS> {
        type Output = Self;
        fn neg(self) -> Self::Output { CryptoInt(Int::ZERO - self.0) }
    }

    impl<const LIMBS: usize> BigIntAlg for CryptoInt<LIMBS> {
        fn zero() -> Self { CryptoInt(Int::ZERO) }
        fn one() -> Self { CryptoInt(Int::ONE) }
        fn from_i32(val: i32) -> Self { CryptoInt(Int::from_i32(val)) }
        fn gcd(&self, other: &Self) -> Self {
            let g_uint = self.0.abs().gcd_vartime(&other.0.abs());
            CryptoInt(Int::new(g_uint.into()))
        }
        fn xgcd(&self, other: &Self) -> (Self, Self, Self) {
            let out = self.0.xgcd(&other.0);
            (CryptoInt(Int::new(out.gcd.into())), CryptoInt(out.x), CryptoInt(out.y))
        }
        fn abs(&self) -> Self { CryptoInt(Int::new(self.0.abs().into())) }
        fn is_zero(&self) -> bool { self.0.is_zero().into() }
        fn sqrt(&self) -> Self {
            let root_uint = self.0.abs().floor_sqrt_vartime();
            CryptoInt(Int::new(root_uint.into()))
        }

        fn ct_eq(&self, other: &Self) -> u32 {
            ((self.0 == other.0) as u32).wrapping_neg()
        }
    }
}
