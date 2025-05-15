/// Given two integers `a` and `b` represented as little endian bytes, compute
/// the value of (a - b) modulo 2^n.
// TODO: make this work when n is not a perfect multiple of the slice len?
pub fn byte_slice_difference(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut a_diff = a.to_vec();
    let mut borrow = false;
    for (i, val) in b.iter().enumerate() {
        (a_diff[i], borrow) = a_diff[i].overflowing_sub(val + (borrow as u8))
    }
    a_diff
}

/// Given an integer `a` of bit length `a_bitlen` represented as
/// little endian bytes, compute the value of `a` when `a` is odd
/// and `a - 1` when `a` is even encoded as binary bits.
pub fn encode_to_odd_binary(a_bits: &mut Vec<u8>, a: &[u8], a_bitlen: usize) {
    let mut flip_bit = (a[0] & 1) as u8;
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
