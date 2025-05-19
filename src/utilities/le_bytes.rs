use super::ct::ct_u32_neq_zero;

/// Given an integer `a` represented as little endian bytes, return the number
/// of leading zeros of the binary representation.
fn le_bytes_leading_zeros(a: &[u8]) -> u32 {
    let mut leading_zeros: u32 = 0;
    let mut mask = u32::MAX;

    for byte in a.iter().rev() {
        leading_zeros += byte.leading_zeros() & mask;
        mask ^= mask & ct_u32_neq_zero(*byte as u32);
    }
    leading_zeros
}

/// Return the bit length of an integer `a` represented as little endian bytes.
pub fn le_bytes_bit_length(a: &[u8]) -> usize {
    (a.len() << 3) - (le_bytes_leading_zeros(a) as usize)
}

/// Compute a - (b + c) where a, b < 256 and c = 0 or 1
#[inline]
fn carrying_sub(x: u8, y: u8, carry: bool) -> (u8, bool) {
    let (a, b) = x.overflowing_sub(y);
    let (c, d) = a.overflowing_sub(carry as u8);
    (c, b | d)
}

/// Given two integers `a` and `b` represented as little endian bytes, compute
/// the value of a = (a - b) modulo 2^n in place.
// TODO: make this work when n is not a perfect multiple of the slice len?
pub fn byte_slice_difference_into(a: &mut [u8], b: &[u8]) {
    let mut borrow = false;
    for (i, val) in b.iter().enumerate() {
        (a[i], borrow) = carrying_sub(a[i], *val, borrow)
    }
}

/// Given an integer `a` of bit length `a_bitlen` represented as
/// little endian bytes, compute the value of `a` when `a` is odd
/// and `a - 1` when `a` is even encoded as binary bits.
pub fn encode_to_odd_binary(a_bits: &mut [u8], a: &[u8], a_bitlen: usize) {
    let mut flip_bit = a[0] & 1;
    for i in 0..a_bitlen {
        // We want to compute the binary of a when a is odd and the
        // binary of (a - 1) when a is even. To do this, we compute
        // each bit of a.
        // When a is odd we set flip_bit to zero and compute a ^ 0 = a
        // When a is even we set flip_bit to 1 and compute ai ^ 1
        // then set flip_bit to flip_bit & (ai ^ 1). This sets flip_bit
        // to zero after encountering ai = 1 terminating the flipping.
        a_bits[i] = ((a[i >> 3] >> (i & 7)) & 1) ^ flip_bit;
        flip_bit &= a_bits[i];
    }
}
