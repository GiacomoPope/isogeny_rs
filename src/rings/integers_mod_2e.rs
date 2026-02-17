use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Trait for ring operations
pub trait Ring:
    Copy
    + Neg<Output = Self>
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Mul<Output = Self>
    + MulAssign
    + Div<Output = Self>
    + DivAssign
{
    /// The length of the encoded representation of the ring element in bytes
    const ENCODED_LENGTH: usize;

    /// Predefined constant element representing the value 0.
    const ZERO: Self;

    /// Return the value x for a given integer x of type `i32`.
    fn from_i32(x: i32) -> Self;

    /// Return the value x for a given integer x of type `u32`.
    fn from_u32(x: u32) -> Self;

    /// Return the value x for a given integer x of type `i64`.
    fn from_i64(x: i64) -> Self;

    /// Return the value x for a given integer x of type `u64`.
    fn from_u64(x: u64) -> Self;

    /// Return 0xFFFFFFFF if this value is invertible in the ring,
    /// or 0x00000000 otherwise.
    fn is_invertible(self) -> u32;

    /// Return `0xFFFFFFFF` if this value is zero,
    /// or `0x00000000` otherwise.
    fn is_zero(self) -> u32;

    /// Return `0xFFFFFFFF` if this value is equal to rhs,
    /// or `0x00000000`otherwise.
    fn equals(self, rhs: &Self) -> u32;

    /// Negate this value.
    fn set_neg(&mut self);

    /// Replace this value with its square.
    fn set_square(&mut self);

    /// If invertible, replace this value with its inverse otherwise set to zero.
    fn set_invert(&mut self);

    /// Compute the square of this value.
    fn square(self) -> Self;

    /// Compute the inverse of this value
    fn invert(self) -> Self;

    /// Set this value to `rhs` if `ctl` is `0xFFFFFFFF`; leave it unchanged if
    /// `ctl` is `0x00000000`.
    /// The value of `ctl` MUST be either `0x00000000` or `0xFFFFFFFF`.
    fn set_cond(&mut self, rhs: &Self, ctl: u32);

    // Encode this value into bytes. Encoding uses little-endian, has
    // a fixed size (for a given field), and is canonical.
    fn encode(self) -> [u8; Self::ENCODED_LENGTH];

    /// Decode the provided bytes into a field element. The source slice
    /// can have arbitrary length; the bytes are interpreted with the
    /// unsigned little-endian convention (no sign bit), with the first half
    /// of the bytes corresponding to x0 and the latter half to x1. For each
    /// resulting integer, the result is reduced modulo the field modulus p.
    /// By definition, this function does not enforce canonicality of the source
    /// value.
    fn decode_reduce(buf: &[u8]) -> Self;
}

/// A macro to define the ring Z / (2^e) Z
///
/// Macro expectations:
/// - A typename for the ring generated
/// - The exponent e such that the ring is Z mod 2^e
///
/// Why a macro?
///
/// This macro almost wasn't needed, and we could have ideally just made
/// `Zmod2e<const E: usize>`. But the issue then is we can't use `Self::N`
///  to define the struct, so then you need redundant Zmod2e<const E: usize, const N: usize>`
/// but this also wasn't enough when integrating the type with the Ring trait because
/// `fn encode(self) -> [u8; Self::ENCODED_LENGTH];` caused further issues.
/// Once I understand generic const expressions better maybe this macro can be
/// removed.
#[macro_export]
macro_rules! define_z_mod_2e_core {
        (
        typename = $typename:ident,
        exponent = $e:literal,
    ) => {

        // ========================================================================
        // Core type implementation
        // ========================================================================

        /// The ring of integers modulo 2^e
        ///
        /// A ring element x is encoded into bytes using the unsigned little-endian convention
        /// over the unique element x in the range [0, ..., 2^e - 1].
        #[derive(Clone, Copy, Debug)]
        pub struct $typename([u64; $typename::N]);

        impl $typename {
            // Number of words needed to represent the ring element.
            const N: usize = ($e + 63) >> 6;

            // Mask used on the top word to ensure e bits are used.
            const HI_MASK: u64 = 0xFFFFFFFFFFFFFFFFu64 >> (63 - (($e - 1) & 63));

            /// Length of a ring element encoded into bytes.
            pub const ENCODED_LENGTH: usize = ($e + 7) >> 3;

            /// The element zero in the ring.
            pub const ZERO: Self = Self([0u64; Self::N]);

            /// Create an element by converting the provided integer.
            pub fn from_i32(x: i32) -> Self {
                Self::from_i64(x as i64)
            }

            /// Create an element by converting the provided integer.
            pub fn from_i64(x: i64) -> Self {
                let sx = (x >> 63) as u64;
                let ax = ((x as u64) ^ sx).wrapping_sub(sx);
                let mut r = Self::from_u64(ax);
                let v = -r;
                r.set_cond(&v, sx as u32);
                r
            }

            /// Create an element by converting the provided integer.
            pub fn from_u32(x: u32) -> Self {
                Self::from_u64(x as u64)
            }

            /// Create an element by converting the provided integer.
            pub fn from_u64(x: u64) -> Self {
                let mut r = Self::ZERO;
                r.0[0] = x;
                r
            }

            /// Set this value to rhs if ctl is 0xFFFFFFFF; leave it unchanged if
            /// ctl is 0x00000000.
            /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
            pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
                let c = (ctl as u64) | ((ctl as u64) << 32);
                for i in 0..Self::N {
                    let wa = self.0[i];
                    let wb = rhs.0[i];
                    self.0[i] = wa ^ (c & (wa ^ wb));
                }
            }

            /// Return 0xFFFFFFFF if this value is zero, or 0x00000000 otherwise.
            pub fn is_zero(self) -> u32 {
                let mut r = self.0[0];
                for i in 1..Self::N {
                    r |= self.0[i];
                }
                (((r | r.wrapping_neg()) >> 63) as u32).wrapping_sub(1)
            }

            /// Return 0xFFFFFFFF if this value is invertible (odd), or 0x00000000 otherwise.
            pub fn is_invertible(self) -> u32 {
                let r = self.0[0] & 1;
                (r as u32).wrapping_neg()
            }

            /// Return 0xFFFFFFFF if values are equal, or 0x00000000 otherwise.
            pub fn equals(self, rhs: &Self) -> u32 {
                let mut r = self.0[0] ^ rhs.0[0];
                for i in 1..Self::N {
                    r |= self.0[i] ^ rhs.0[i];
                }
                (((r | r.wrapping_neg()) >> 63) as u32).wrapping_sub(1)
            }

            /// Add `rhs` to this value.
            fn set_add(&mut self, rhs: &Self) {
                let (d, mut cc) = fp2::utils64::addcarry_u64(self.0[0], rhs.0[0], 0);
                self.0[0] = d;
                for i in 1..Self::N {
                    let (d, ee) = fp2::utils64::addcarry_u64(self.0[i], rhs.0[i], cc);
                    self.0[i] = d;
                    cc = ee;
                }
                self.0[Self::N - 1] &= Self::HI_MASK;
            }

            /// Subtract `rhs` from this value.
            fn set_sub(&mut self, rhs: &Self) {
                let (d, mut cc) = fp2::utils64::subborrow_u64(self.0[0], rhs.0[0], 0);
                self.0[0] = d;
                for i in 1..Self::N {
                    let (d, ee) = fp2::utils64::subborrow_u64(self.0[i], rhs.0[i], cc);
                    self.0[i] = d;
                    cc = ee;
                }
                self.0[Self::N - 1] &= Self::HI_MASK;
            }

            /// Negate this value
            fn set_neg(&mut self) {
                let mut cc = 1;
                for i in 0..Self::N {
                    let (d, ee) = fp2::utils64::addcarry_u64(!self.0[i], 0, cc);
                    self.0[i] = d;
                    cc = ee;
                }
                self.0[Self::N - 1] &= Self::HI_MASK;
            }

            /// Multiply this value by `rhs`.
            fn set_mul(&mut self, rhs: &Self) {
                let mut t = [0u64; Self::N];
                for i in 0..Self::N {
                    let f = rhs.0[i];
                    let (d, mut cc) = fp2::utils64::umull_add(f, self.0[0], t[i]);
                    t[i] = d;
                    for j in 1..(Self::N - i) {
                        let (d, hi) = fp2::utils64::umull_add2(f, self.0[j], t[i + j], cc);
                        t[i + j] = d;
                        cc = hi;
                    }
                }
                t[Self::N - 1] &= Self::HI_MASK;
                self.0[..].copy_from_slice(&t);
            }

            /// Square this value
            pub fn set_square(&mut self) {
                let mut t = [0u64; Self::N];
                for i in 0..(Self::N >> 1) {
                    let f = self.0[i];
                    let (d, mut cc) = fp2::utils64::umull_add(f, self.0[i + 1], t[(i << 1) + 1]);
                    t[(i << 1) + 1] = d;
                    for j in (i + 2)..(Self::N - i) {
                        let (d, hi) = fp2::utils64::umull_add2(f, self.0[j], t[i + j], cc);
                        t[i + j] = d;
                        cc = hi;
                    }
                }
                let mut cc = 0;
                for i in 1..Self::N {
                    let w = t[i];
                    let ee = w >> 63;
                    t[i] = (w << 1) | cc;
                    cc = ee;
                }
                let mut cc = 0;
                for i in 0..(Self::N >> 1) {
                    let (lo, hi) = fp2::utils64::umull(self.0[i], self.0[i]);
                    let (d0, ee) = fp2::utils64::addcarry_u64(lo, t[i << 1], cc);
                    let (d1, ee) = fp2::utils64::addcarry_u64(hi, t[(i << 1) + 1], ee);
                    t[i << 1] = d0;
                    t[(i << 1) + 1] = d1;
                    cc = ee;
                }
                if (Self::N & 1) != 0 {
                    let w = self.0[Self::N >> 1];
                    let w = w.wrapping_mul(w);
                    let (d, _) = fp2::utils64::addcarry_u64(w, t[Self::N - 1], cc);
                    t[Self::N - 1] = d;
                }
                t[Self::N - 1] &= Self::HI_MASK;
                self.0[..].copy_from_slice(&t);
            }

            /// Return the square of this value
            pub fn square(self) -> Self {
                let mut r = self;
                r.set_square();
                r
            }

            /// Invert this value, if this value has no inverse (it is even) then the value
            /// is set to zero.
            pub fn set_invert(&mut self) {
                // If a value is even we can not compute the inverse and zero is returned
                let is_even = !self.is_invertible();

                // Invert modulo 2^64. This uses Hensel lifting:
                // if x*y = 1 mod r, then x*(y*(2 - y*x)) = 1 mod r^2
                let m0 = self.0[0];
                let y = 2u64.wrapping_sub(m0);
                let y = y.wrapping_mul(2u64.wrapping_sub(y.wrapping_mul(m0)));
                let y = y.wrapping_mul(2u64.wrapping_sub(y.wrapping_mul(m0)));
                let y = y.wrapping_mul(2u64.wrapping_sub(y.wrapping_mul(m0)));
                let y = y.wrapping_mul(2u64.wrapping_sub(y.wrapping_mul(m0)));
                let y = y.wrapping_mul(2u64.wrapping_sub(y.wrapping_mul(m0)));

                // Special case when using only one limb.
                if Self::N == 1 {
                    self.0[0] = y & Self::HI_MASK;
                    return;
                }

                let mut yy = [0u64; Self::N];
                yy[0] = y;
                let mut n = 1;
                while n < Self::N {
                    // Compute new length.
                    let mut n2 = n << 1;
                    if n2 > Self::N {
                        n2 = Self::N;
                    }
                    n = n2;

                    // t <- y*x
                    let mut tt = [0u64; Self::N];
                    for i in 0..n {
                        let f = self.0[i];
                        let (d, mut cc) = fp2::utils64::umull_add(f, yy[0], tt[i]);
                        tt[i] = d;
                        for j in 1..(n - i) {
                            let (d, hi) = fp2::utils64::umull_add2(f, yy[j], tt[i + j], cc);
                            tt[i + j] = d;
                            cc = hi;
                        }
                    }

                    // t <- 2 - t
                    let (d, mut cc) = fp2::utils64::subborrow_u64(2u64, tt[0], 0);
                    tt[0] = d;
                    for i in 1..n {
                        let (d, ee) = fp2::utils64::subborrow_u64(0u64, tt[i], cc);
                        tt[i] = d;
                        cc = ee;
                    }

                    // y <- y*t
                    let mut uu = [0u64; Self::N];
                    for i in 0..n {
                        let f = yy[i];
                        let (d, mut cc) = fp2::utils64::umull_add(f, tt[0], uu[i]);
                        uu[i] = d;
                        for j in 1..(n - i) {
                            let (d, hi) = fp2::utils64::umull_add2(f, tt[j], uu[i + j], cc);
                            uu[i + j] = d;
                            cc = hi;
                        }
                    }
                    yy[..n].copy_from_slice(&uu[..n]);
                }

                yy[Self::N - 1] &= Self::HI_MASK;
                self.0[..].copy_from_slice(&yy);
                self.set_cond(&Self::ZERO, is_even);
            }

            /// Return the inverse of this value, if `self` is even, then zero is returned
            pub fn invert(self) -> Self {
                let mut r = self;
                r.set_invert();
                r
            }

            /// Divide this value by `rhs`. If `rhs` is zero, then this sets this value
            /// to zero.
            fn set_div(&mut self, rhs: &Self) {
                self.set_mul(&rhs.invert());
            }

            /// Encode this value into bytes. Encoding uses little-endian, has
            /// a fixed size (for a given field), and is canonical.
            pub fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
                let mut r = [0u8; Self::ENCODED_LENGTH];
                for i in 0..(Self::N - 1) {
                    r[(8 * i)..(8 * (i + 1))].copy_from_slice(&self.0[i].to_le_bytes());
                }
                r[(8 * (Self::N - 1))..]
                    .copy_from_slice(&self.0[Self::N - 1].to_le_bytes()[..(Self::ENCODED_LENGTH - 8 * (Self::N - 1))]);
                r
            }

            /// Decode the provided bytes into a field element. The source slice
            /// can have arbitrary length; the bytes are interpreted with the
            /// unsigned little-endian convention (no sign bit), and the resulting
            /// integer is reduced modulo the field modulus 2^e. By definition, this
            /// function does not enforce canonicality of the source value.
            pub fn decode_reduce(buf: &[u8]) -> Self {
                let mut r = Self::ZERO;
                let mut i = 0;
                while i < Self::N && (8 * (i + 1)) < buf.len() {
                    r.0[i] =
                        u64::from_le_bytes(*<&[u8; 8]>::try_from(&buf[(8 * i)..(8 * (i + 1))]).unwrap());
                    i += 1;
                }
                if i < Self::N {
                    let mut j = i * 8;
                    while j < buf.len() {
                        r.0[i] |= (buf[j] as u64) << (8 * (j & 7));
                        j += 1;
                    }
                }
                r.0[Self::N - 1] &= Self::HI_MASK;
                r
            }
        }

        // ========================================================================
        // Operator trait implementations
        // ========================================================================

        impl ::core::ops::Add<$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn add(self, other: $typename) -> $typename {
                let mut r = self;
                r.set_add(&other);
                r
            }
        }

        impl ::core::ops::Add<&$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn add(self, other: &$typename) -> $typename {
                let mut r = self;
                r.set_add(other);
                r
            }
        }

        impl ::core::ops::Add<$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn add(self, other: $typename) -> $typename {
                let mut r = *self;
                r.set_add(&other);
                r
            }
        }

        impl ::core::ops::Add<&$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn add(self, other: &$typename) -> $typename {
                let mut r = *self;
                r.set_add(other);
                r
            }
        }

        impl ::core::ops::AddAssign<$typename> for $typename {
            #[inline(always)]
            fn add_assign(&mut self, other: $typename) {
                self.set_add(&other);
            }
        }

        impl ::core::ops::AddAssign<&$typename> for $typename {
            #[inline(always)]
            fn add_assign(&mut self, other: &$typename) {
                self.set_add(other);
            }
        }

        impl ::core::ops::Div<$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn div(self, other: $typename) -> $typename {
                let mut r = self;
                r.set_div(&other);
                r
            }
        }

        impl ::core::ops::Div<&$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn div(self, other: &$typename) -> $typename {
                let mut r = self;
                r.set_div(other);
                r
            }
        }

        impl ::core::ops::Div<$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn div(self, other: $typename) -> $typename {
                let mut r = *self;
                r.set_div(&other);
                r
            }
        }

        impl ::core::ops::Div<&$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn div(self, other: &$typename) -> $typename {
                let mut r = *self;
                r.set_div(other);
                r
            }
        }

        impl ::core::ops::DivAssign<$typename> for $typename {
            #[inline(always)]
            fn div_assign(&mut self, other: $typename) {
                self.set_div(&other);
            }
        }

        impl ::core::ops::DivAssign<&$typename> for $typename {
            #[inline(always)]
            fn div_assign(&mut self, other: &$typename) {
                self.set_div(other);
            }
        }

        impl ::core::ops::Mul<$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn mul(self, other: $typename) -> $typename {
                let mut r = self;
                r.set_mul(&other);
                r
            }
        }

        impl ::core::ops::Mul<&$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn mul(self, other: &$typename) -> $typename {
                let mut r = self;
                r.set_mul(other);
                r
            }
        }

        impl ::core::ops::Mul<$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn mul(self, other: $typename) -> $typename {
                let mut r = *self;
                r.set_mul(&other);
                r
            }
        }

        impl ::core::ops::Mul<&$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn mul(self, other: &$typename) -> $typename {
                let mut r = *self;
                r.set_mul(other);
                r
            }
        }

        impl ::core::ops::MulAssign<$typename> for $typename {
            #[inline(always)]
            fn mul_assign(&mut self, other: $typename) {
                self.set_mul(&other);
            }
        }

        impl ::core::ops::MulAssign<&$typename> for $typename {
            #[inline(always)]
            fn mul_assign(&mut self, other: &$typename) {
                self.set_mul(other);
            }
        }

        impl ::core::ops::Neg for $typename {
            type Output = $typename;
            #[inline(always)]
            fn neg(self) -> $typename {
                let mut r = self;
                r.set_neg();
                r
            }
        }

        impl ::core::ops::Neg for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn neg(self) -> $typename {
                let mut r = *self;
                r.set_neg();
                r
            }
        }

        impl ::core::ops::Sub<$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn sub(self, other: $typename) -> $typename {
                let mut r = self;
                r.set_sub(&other);
                r
            }
        }

        impl ::core::ops::Sub<&$typename> for $typename {
            type Output = $typename;
            #[inline(always)]
            fn sub(self, other: &$typename) -> $typename {
                let mut r = self;
                r.set_sub(other);
                r
            }
        }

        impl ::core::ops::Sub<$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn sub(self, other: $typename) -> $typename {
                let mut r = *self;
                r.set_sub(&other);
                r
            }
        }

        impl ::core::ops::Sub<&$typename> for &$typename {
            type Output = $typename;
            #[inline(always)]
            fn sub(self, other: &$typename) -> $typename {
                let mut r = *self;
                r.set_sub(other);
                r
            }
        }

        impl ::core::ops::SubAssign<$typename> for $typename {
            #[inline(always)]
            fn sub_assign(&mut self, other: $typename) {
                self.set_sub(&other);
            }
        }

        impl ::core::ops::SubAssign<&$typename> for $typename {
            #[inline(always)]
            fn sub_assign(&mut self, other: &$typename) {
                self.set_sub(other);
            }
        }

        // ========================================================================
        // Implement the generic trait Mod2E
        // ========================================================================

        impl super::Ring for $typename {
            // Reexport constants for base field Trait
            const ENCODED_LENGTH: usize = $typename::ENCODED_LENGTH;
            const ZERO: Self = Self::ZERO;

            fn from_i32(x: i32) -> Self {
                Self::from_i32(x)
            }
            fn from_u32(x: u32) -> Self {
                Self::from_u32(x)
            }
            fn from_i64(x: i64) -> Self {
                Self::from_i64(x)
            }
            fn from_u64(x: u64) -> Self {
                Self::from_u64(x)
            }
            fn is_invertible(self) -> u32 {
                self.is_invertible()
            }
            fn is_zero(self) -> u32 {
                self.is_zero()
            }
            fn equals(self, rhs: &Self) -> u32 {
                self.equals(rhs)
            }
            fn set_neg(&mut self) {
                self.set_neg()
            }
            fn set_square(&mut self) {
                self.set_square()
            }
            fn set_invert(&mut self) {
                self.set_invert()
            }
            fn square(self) -> Self {
                self.square()
            }
            fn invert(self) -> Self {
                self.invert()
            }
            fn set_cond(&mut self, rhs: &Self, ctl: u32) {
                self.set_cond(rhs, ctl)
            }
            fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
                self.encode()
            }
            fn decode_reduce(buf: &[u8]) -> Self {
                Self::decode_reduce(buf)
            }
        }

    }
}

#[cfg(test)]
mod tests {
    use num_bigint::{BigInt, Sign, ToBigInt};
    use sha2::{Digest, Sha256};

    // Number of times each randomised test is performed
    static TEST_ITERATIONS: usize = 100;

    // Type alias for the specific bit size used in tests
    static EXPONENT: usize = 248;
    define_z_mod_2e_core!(typename = ZMod2_248, exponent = 248,);

    /// Helper function to generate test vectors using SHA256
    fn generate_test_vectors(index: usize) -> ([u8; 256], [u8; 256]) {
        let mut va = [0u8; 256];
        let mut vb = [0u8; 256];
        let mut sh = Sha256::new();

        for j in 0..8 {
            sh.update(((256 * index + 8 * j + 0) as u64).to_le_bytes());
            va[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
        }
        for j in 0..8 {
            sh.update(((256 * index + 8 * j + 1) as u64).to_le_bytes());
            vb[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
        }

        (va, vb)
    }

    /// Helper to create BigInt modulus 2^E
    fn zmod_2e_modulus_bigint(e: usize) -> BigInt {
        1u32.to_bigint().unwrap() << e
    }

    /// Helper to decode bytes to BigInt
    fn bytes_to_bigint(bytes: &[u8]) -> BigInt {
        BigInt::from_bytes_le(Sign::Plus, bytes)
    }

    // ========================================================================
    // Encode/Decode Tests
    // ========================================================================

    #[test]
    fn test_encode_decode() {
        for i in 0..TEST_ITERATIONS {
            let (va, _) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let za = bytes_to_bigint(&va);

            let encoded = a.encode();
            let decoded_bigint = bytes_to_bigint(&encoded);
            let expected = &za % &zmod_2e_modulus_bigint(EXPONENT);

            assert_eq!(
                decoded_bigint, expected,
                "Encode/decode mismatch at iteration {}",
                i
            );
        }
    }

    #[test]
    fn test_encode_decode_with_zero() {
        let va = [0u8; 256];
        let a = ZMod2_248::decode_reduce(&va);
        let encoded = a.encode();
        let decoded_bigint = bytes_to_bigint(&encoded);

        assert_eq!(decoded_bigint, BigInt::from(0u32), "Zero encoding failed");
    }

    // ========================================================================
    // Addition Tests
    // ========================================================================

    #[test]
    fn test_addition() {
        for i in 0..TEST_ITERATIONS {
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);
            let za = bytes_to_bigint(&va);
            let zb = bytes_to_bigint(&vb);

            let c = a + b;
            let encoded = c.encode();
            let result = bytes_to_bigint(&encoded);
            let expected = (&za + &zb) % &zmod_2e_modulus_bigint(EXPONENT);

            assert_eq!(result, expected, "Addition failed at iteration {}", i);
        }
    }

    #[test]
    fn test_addition_with_zero() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);
        let zero = ZMod2_248::ZERO;

        let result = a + zero;
        assert_eq!(a.equals(&result), 0xFFFFFFFF, "Addition with zero failed");
    }

    #[test]
    fn test_addition_commutativity() {
        for i in 0..TEST_ITERATIONS {
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);

            let c1 = a + b;
            let c2 = b + a;

            assert_eq!(
                c1.equals(&c2),
                0xFFFFFFFF,
                "Addition is not commutative at iteration {}",
                i
            );
        }
    }

    // ========================================================================
    // Subtraction Tests
    // ========================================================================

    #[test]
    fn test_subtraction() {
        for i in 0..TEST_ITERATIONS {
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);
            let za = bytes_to_bigint(&va);
            let zb = bytes_to_bigint(&vb);

            let c = a - b;
            let encoded = c.encode();
            let result = bytes_to_bigint(&encoded);

            let zp = zmod_2e_modulus_bigint(EXPONENT);
            let zpz = &zp << 64;
            let expected = ((&zpz + &za) - (&zb % &zp)) % &zp;

            assert_eq!(result, expected, "Subtraction failed at iteration {}", i);
        }
    }

    #[test]
    fn test_subtraction_self() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);

        let result = a - a;
        assert_eq!(
            result.is_zero(),
            0xFFFFFFFF,
            "Subtraction of self should be zero"
        );
    }

    // ========================================================================
    // Negation Tests
    // ========================================================================

    #[test]
    fn test_negation() {
        for i in 0..TEST_ITERATIONS {
            let (va, _) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let za = bytes_to_bigint(&va);

            let c = -a;
            let encoded = c.encode();
            let result = bytes_to_bigint(&encoded);

            let zp = zmod_2e_modulus_bigint(EXPONENT);
            let expected = (&zp - (&za % &zp)) % &zp;

            assert_eq!(result, expected, "Negation failed at iteration {}", i);
        }
    }

    #[test]
    fn test_negation_involution() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);

        let neg_a = -a;
        let neg_neg_a = -neg_a;

        assert_eq!(
            a.equals(&neg_neg_a),
            0xFFFFFFFF,
            "Double negation should return original"
        );
    }

    #[test]
    fn test_negation_zero() {
        let zero = ZMod2_248::ZERO;
        let neg_zero = -zero;

        assert_eq!(
            neg_zero.is_zero(),
            0xFFFFFFFF,
            "Negation of zero should be zero"
        );
    }

    // ========================================================================
    // Multiplication Tests
    // ========================================================================

    #[test]
    fn test_multiplication() {
        for i in 0..TEST_ITERATIONS {
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);
            let za = bytes_to_bigint(&va);
            let zb = bytes_to_bigint(&vb);

            let c = a * b;
            let encoded = c.encode();
            let result = bytes_to_bigint(&encoded);
            let expected = (&za * &zb) % &zmod_2e_modulus_bigint(EXPONENT);

            assert_eq!(result, expected, "Multiplication failed at iteration {}", i);
        }
    }

    #[test]
    fn test_multiplication_by_zero() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);
        let zero = ZMod2_248::ZERO;

        let result = a * zero;
        assert_eq!(
            result.is_zero(),
            0xFFFFFFFF,
            "Multiplication by zero should be zero"
        );
    }

    #[test]
    fn test_multiplication_commutativity() {
        for i in 0..TEST_ITERATIONS {
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);

            let c1 = a * b;
            let c2 = b * a;

            assert_eq!(
                c1.equals(&c2),
                0xFFFFFFFF,
                "Multiplication is not commutative at index {}",
                i
            );
        }
    }

    // ========================================================================
    // Squaring Tests
    // ========================================================================

    #[test]
    fn test_square() {
        for i in 0..TEST_ITERATIONS {
            let (va, _) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let za = bytes_to_bigint(&va);

            let c = a.square();
            let encoded = c.encode();
            let result = bytes_to_bigint(&encoded);
            let expected = (&za * &za) % &zmod_2e_modulus_bigint(EXPONENT);

            assert_eq!(result, expected, "Square failed at iteration {}", i);
        }
    }

    #[test]
    fn test_square_equals_multiply_self() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);

        let squared = a.square();
        let multiplied = a * a;

        assert_eq!(
            squared.equals(&multiplied),
            0xFFFFFFFF,
            "Square should equal self-multiplication"
        );
    }

    // ========================================================================
    // Invertibility tests
    // ========================================================================

    #[test]
    fn test_is_invertible() {
        for i in 0..TEST_ITERATIONS {
            let (va, _) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);

            if (va[0] & 1) == 0 {
                assert_eq!(a.is_invertible(), 0);
            } else {
                assert_eq!(a.is_invertible(), 0xFFFFFFFF);
            }
        }
    }

    // ========================================================================
    // Division Tests
    // ========================================================================

    #[test]
    fn test_division() {
        for i in 0..TEST_ITERATIONS {
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);

            let c = a / b;
            let d = c * b;

            // Branch depending on whether b is invertible or not
            if (vb[0] & 1) == 0 {
                assert_eq!(
                    c.is_zero(),
                    0xFFFFFFFF,
                    "Inversion by even value should be zero"
                );
            } else {
                assert_eq!(a.equals(&d), 0xFFFFFFFF, "Division verification failed");
            }
        }
    }

    #[test]
    fn test_division_by_one() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);
        let one = ZMod2_248::from_u64(1);

        let result = a / one;
        assert_eq!(
            a.equals(&result),
            0xFFFFFFFF,
            "Division by one should return original"
        );
    }

    // ========================================================================
    // Inversion Tests
    // ========================================================================

    #[test]
    fn test_inversion() {
        for i in 0..TEST_ITERATIONS {
            let (_, vb) = generate_test_vectors(i);
            let b = ZMod2_248::decode_reduce(&vb);

            let one = ZMod2_248::from_u64(1);
            let b_inv = b.invert();
            let product = b * b_inv;

            // Ensure b is invertible (odd)
            if (vb[0] & 1) == 0 {
                assert_eq!(
                    b_inv.is_zero(),
                    0xFFFFFFFF,
                    "inversion of even value should be zero"
                );
            } else {
                assert_eq!(
                    product.equals(&one),
                    0xFFFFFFFF,
                    "Inversion failed at iteration"
                );
            }
        }
    }

    // ========================================================================
    // Equality and Zero Tests
    // ========================================================================

    #[test]
    fn test_equals_reflexive() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);

        assert_eq!(a.equals(&a), 0xFFFFFFFF, "Equals should be reflexive");
    }

    #[test]
    fn test_equals_symmetric() {
        let (va, _) = generate_test_vectors(0);
        let a = ZMod2_248::decode_reduce(&va);
        let b = ZMod2_248::decode_reduce(&va);

        assert_eq!(a.equals(&b), 0xFFFFFFFF, "Equals should be symmetric");
        assert_eq!(b.equals(&a), 0xFFFFFFFF, "Equals should be symmetric");
    }

    #[test]
    fn test_is_zero_for_zero() {
        let zero = ZMod2_248::ZERO;
        assert_eq!(
            zero.is_zero(),
            0xFFFFFFFF,
            "Zero should be detected as zero"
        );
    }

    #[test]
    fn test_is_zero_for_nonzero() {
        let (va, _) = generate_test_vectors(0);
        let mut a = ZMod2_248::decode_reduce(&va);

        // Ensure it's not zero
        if a.is_zero() == 0xFFFFFFFF {
            a.0[0] = 1;
        }

        assert_eq!(a.is_zero(), 0, "Non-zero should not be detected as zero");
    }

    // ========================================================================
    // Comprehensive Property Tests
    // ========================================================================

    #[test]
    fn test_comprehensive_operations() {
        for i in 0..TEST_ITERATIONS {
            // Test case with specific known values
            let (va, vb) = generate_test_vectors(i);
            let a = ZMod2_248::decode_reduce(&va);
            let b = ZMod2_248::decode_reduce(&vb);

            // Test that (a + b) - b == a
            let sum = a + b;
            let diff = sum - b;
            assert_eq!(
                a.equals(&diff),
                0xFFFFFFFF,
                "Addition/subtraction identity failed at index {}",
                i
            );

            // Test that (a * b) / b == a (when b is invertible)
            let mut b_inv = b;
            if (vb[0] & 1) == 0 {
                b_inv.0[0] |= 1u64;
            }
            let prod = a * b_inv;
            let quot = prod / b_inv;
            assert_eq!(
                a.equals(&quot),
                0xFFFFFFFF,
                "Multiplication/division identity failed at index {}",
                i
            );
        }
    }
}
