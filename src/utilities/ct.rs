/// Returns `0x00..00` when `x = 0` and `0xFF..FF` otherwise in constant time.
pub fn ct_u32_neq_zero(x: u32) -> u32 {
    // (x | -x) will be zero only when y is zero.
    let q = x | x.wrapping_neg();

    // Similarly, r = (q | - q) >> 31 will be 0 only if x is zero, otherwise it will be 1
    // So then -r is  0xFF..FF when x != 0 and 0x00..00 otherwise.
    ((q | q.wrapping_neg()) >> 31).wrapping_neg()
}

/// Returns `0x00..00` when `x = 0` and `0xFF..FF` otherwise in constant time.
fn ct_u32_eq_zero(x: u32) -> u32 {
    !ct_u32_neq_zero(x)
}

/// Returns `0xFF..FF` when two values are equal and zero otherwise
pub fn ct_u32_eq(x: u32, y: u32) -> u32 {
    ct_u32_eq_zero(x ^ y)
}
