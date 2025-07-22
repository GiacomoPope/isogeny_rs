
/// Given an integer `a` represented with little endian u64 words, return the number
/// of leading zeros of the binary representation.
fn bn_leading_zeros_vartime(a: &[u64]) -> u32 {
    let mut leading_zeros: u32 = 0;
    for word in a.iter().rev() {
        leading_zeros += word.leading_zeros();
        if *word != 0 {
            break;
        }
    }

    leading_zeros
}

/// Return the bit length of an integer `a` represented as little endian bytes.
pub fn bn_bit_length_vartime(a: &[u64]) -> usize {
    (a.len() << 6) - (bn_leading_zeros_vartime(a) as usize)
}


// This is a small wrapper to use a bn as constant
#[macro_export]
macro_rules! create_bn {

    ($typename:ident, $maxlimbs:expr) => {

        pub trait BnTrait {
            fn from_factorisation(factorisation: &[(usize, usize)]) -> $typename;
            fn from_primepower(x: usize, e : usize) -> $typename;
        }

        
        use std::ops::Mul;
        use std::cmp::{PartialEq, PartialOrd, Ordering};
        use std::cmp::Ordering::{Greater, Equal, Less};

        use forward_ref::forward_ref_binop;
        
        use fp2::utils64::{addcarry_u64, umull, umull_add};

        #[derive(Clone, Copy, Debug)]
        pub struct $typename {
            pub limbs: [u64; $maxlimbs],
            pub len: usize,
        }

        impl $typename {
            pub fn new() -> Self {
                $typename {
                    limbs: [0u64; $maxlimbs],
                    len: 0,
                }
            }
        }

        impl $typename {
            pub const fn default() -> Self {
                Self{ limbs: [0; $maxlimbs], len: 1}
            }

            pub const fn One() -> Self {
                let mut one = Self::default();
                one.limbs[0] = 1;

                one
            }
        }

        impl BnTrait for $typename {
            fn from_primepower(x: usize, e : usize) -> Self {
                // If x^e fits inside a word, then we can just finish here.
               let n_bitlength = ((usize::BITS - x.leading_zeros()) as usize) * e;
               if n_bitlength <= 64 {
                   let n_u64 = (x as u64).pow(e as u32);
                   return Self::from(n_u64);
               }

               // If x^e fits into a u128, then we can also handle that as a special-case too
               if n_bitlength <= 128 {
                   let n_u128 = (x as u128).pow(e as u32);
                   let n0 = n_u128 as u64;
                   let n1 = (n_u128 >> 64) as u64;
                   return Self::from((n0, n1));
               }

               // Otherwise we halve the exponent and compute two big numbers which
               // we then need to combine. We use that n = (x^(e / 2))^2 when n is
               // even and n = x * (x^(e / 2))^2 when e is odd.
               let e_half = e / 2;
               let n_lo = Self::from_primepower(x, e_half);
               let n_lo_sqr = &n_lo * &n_lo;

               // Calculate x^e based on whether e is even or odd
               let n = if e % 2 == 0 {
                   // If e is even, n = (x^(e/2))^2
                   n_lo_sqr
               } else {
                   // If e is odd, n = x * (x^(e/2))^2
                   &n_lo_sqr * (x as u64)
               };

               n
           }

            fn from_factorisation(factorisation: &[(usize, usize)]) -> Self {
               let (p0, e0) = factorisation[0];
               let mut n = Self::from_primepower(p0, e0);
               for (pi, ei) in factorisation.iter().skip(1) {
                   let ni = Self::from_primepower(*pi, *ei);
                   n = &n * &ni;
               }

               n
           }
        }

        impl From<u64> for $typename {  
            fn from(x: u64) -> Self {
                let mut new = Self::default();
                new.limbs[0] = x;

                new
            }
        }

        impl From<&[u64]> for $typename {  
            fn from(x: &[u64]) -> Self {
                assert!(x.len() <= $maxlimbs);

                let mut new = Self::default();
                new.len = x.len();
                new.limbs[..x.len()].clone_from_slice(x);

                println!("{:?}", new);
                new
            }
        }

        impl From<(u64, u64)> for $typename {
            fn from(x: (u64, u64)) -> Self {
                let mut new = Self::default();
                new.len = 2;
                new.limbs[0] = x.0;
                new.limbs[1] = x.1;

                new
            }
        }


        impl AsRef<[u64]> for $typename {
            fn as_ref(&self) -> &[u64] {
                &self.limbs[..self.len]
            }
        }

        impl Mul<u64> for $typename {
            type Output = $typename;

            fn mul(self, rhs: u64) -> Self {
                // If a has length 1 then we can do a single double wide multiplication.
                if self.len == 1 {
                    let (n0, n1) = umull(self.limbs[0], rhs);
                    if n1 == 0 {
                        return Self::from(n0);
                    }
                    return Self::from((n0, n1));
                }

                // Compute b * (a0 + a1 * 2^64 + a1 * 2^128 + ...)
                let mut res = $typename::default();
                res.len = self.len;
                let mut carry: u64;
                (res.limbs[0], carry) = umull(self.limbs[0], rhs);
                for i in 1..self.len {
                    (res.limbs[i], carry) = umull_add(self.limbs[i], rhs, carry);
                }

                // If there was a carry at the end of the multiplication, we push this carry
                // to the end of the vector.
                if carry != 0 {
                    res.limbs[res.len] = carry;
                    res.len += 1;
                }

                res
            }
        }

        impl Mul<u64> for &$typename {
            type Output = $typename;

            fn mul(self, rhs : u64) -> $typename {
                *self*rhs
            }
        }


        impl Mul for $typename {
            type Output = $typename;

            fn mul(self, rhs: Self) -> Self {
                assert!(self.len + rhs.len <= $maxlimbs);

                // Assume that the length of b is smaller than the length of a for logic.
                if rhs.len > self.len {
                    return rhs * self;
                }
            
                // If a has length 1 then so does b and we can do simple multiplication.
                if self.len == 1 {
                    let (n0, n1) = umull(self.limbs[0], rhs.limbs[0]);
                    if n1 == 0 {
                        return Self::from(n0);
                    }
                    return Self::from((n0, n1));
                }
            
                // If b has length 1 then we can do a simple scalar multiplication
                if rhs.len == 1 {
                    return self * rhs.limbs[0];
                }
            
                // General case
                let mut res = vec![0u64; self.len + rhs.len];
                let mut carry: u64 = 0;
                let mut cc: u8 = 0;
            
                // TODO: clean up the carry here...
                for i in 0..self.len {
                    for j in 0..rhs.len {
                        let (lo, hi) = umull_add(self.limbs[i], rhs.limbs[j], carry);
                        (res[i + j], cc) = addcarry_u64(res[i + j], lo, cc);
                        carry = hi;
                    }
                    res[i + rhs.len] += carry;
                    carry = 0;
                }
            
                Self::from(&res[..])
            }
        }

        forward_ref_binop!(impl Mul, mul for $typename, $typename);

        impl PartialEq for $typename {
            fn eq(&self, other: &Self) -> bool {
                if self.len != other.len {
                    return false;
                }

                for i in 0..self.len {
                    if self.limbs[i] != other.limbs[i] {
                        return false;
                    }
                }

                true
            }
        }

        impl PartialOrd for $typename {

            fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            
                if self.len > other.len {
                    return Some(Greater);
                } else if self.len < other.len {
                    return Some(Less);
                }
            
                for i in (0..self.len).rev() {
                    if self.limbs[i] > other.limbs[i] {
                        return Some(Greater);
                    } else if self.limbs[i] < other.limbs[i] {
                        return Some(Less);
                    }
                }
            
                Some(Equal)
            }
        }
    };
}

