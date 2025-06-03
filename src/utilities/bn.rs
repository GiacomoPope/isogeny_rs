use fp2::utils64::{addcarry_u64, umull, umull_add};

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

fn mul_bn_by_u64_vartime(a: &[u64], b: u64) -> Vec<u64> {
    // If a has length 1 then we can do a single double wide multiplication.
    if a.len() == 1 {
        let (n0, n1) = umull(a[0], b);
        if n1 == 0 {
            return vec![n0];
        }
        return vec![n0, n1];
    }

    // Compute b * (a0 + a1 * 2^64 + a1 * 2^128 + ...)
    let mut res = vec![0u64; a.len()];
    let mut carry: u64;
    (res[0], carry) = umull(a[0], b);
    for i in 1..a.len() {
        (res[i], carry) = umull_add(a[i], b, carry);
    }

    // If there was a carry at the end of the multiplication, we push this carry
    // to the end of the vector.
    if carry != 0 {
        res.push(carry);
    }

    res
}

/// Given two integers represented as u64 words (little endian) compute their
/// product.
fn mul_bn_vartime(a: &[u64], b: &[u64]) -> Vec<u64> {
    // Assume that the length of b is smaller than the length of a for logic.
    if b.len() > a.len() {
        return mul_bn_vartime(b, a);
    }

    // If a has length 1 then so does b and we can do simple multiplication.
    if a.len() == 1 {
        let (n0, n1) = umull(a[0], b[0]);
        if n1 == 0 {
            return vec![n0];
        }
        return vec![n0, n1];
    }

    // If b has length 1 then we can do a simple scalar multiplication
    if b.len() == 1 {
        return mul_bn_by_u64_vartime(a, b[0]);
    }

    // General case
    let mut res = vec![0u64; a.len() + b.len()];
    let mut carry: u64 = 0;
    let mut cc: u8 = 0;

    // TODO: clean up the carry here...
    for i in 0..a.len() {
        for j in 0..b.len() {
            let (lo, hi) = umull_add(a[i], b[j], carry);
            (res[i + j], cc) = addcarry_u64(res[i + j], lo, cc);
            carry = hi;
        }
        res[i + b.len()] += carry;
        carry = 0;
    }

    res
}

/// Represent an integer n = x^e as little endian u64 words. Both x and e are assumed
/// to be public, no constant-time guarentees are made.
pub fn prime_power_to_bn_vartime(x: usize, e: usize) -> Vec<u64> {
    // If x^e fits inside a word, then we can just finish here.
    let n_bitlength = ((usize::BITS - x.leading_zeros()) as usize) * e;
    if n_bitlength <= 64 {
        let n_u64 = (x as u64).pow(e as u32);
        return vec![n_u64];
    }

    // If x^e fits into a u128, then we can also handle that as a special-case too
    if n_bitlength <= 128 {
        let n_u128 = (x as u128).pow(e as u32);
        let n0 = n_u128 as u64;
        let n1 = (n_u128 >> 64) as u64;
        return vec![n0, n1];
    }

    // Otherwise we halve the exponent and compute two big numbers which
    // we then need to combine. We use that n = (x^(e / 2))^2 when n is
    // even and n = x * (x^(e / 2))^2 when e is odd.
    let e_half = e / 2;
    let n_lo = prime_power_to_bn_vartime(x, e_half);
    let n_lo_sqr = mul_bn_vartime(&n_lo, &n_lo);

    // Calculate x^e based on whether e is even or odd
    let mut n = if e % 2 == 0 {
        // If e is even, n = (x^(e/2))^2
        n_lo_sqr
    } else {
        // If e is odd, n = x * (x^(e/2))^2
        mul_bn_by_u64_vartime(&n_lo_sqr, x as u64)
    };

    // Remove trailing zeros: TODO better.
    while *n.last().unwrap() == 0 {
        n.pop();
    }

    n
}

/// Represent an integer n = pi^ei as little endian u64 words given the factorisation
/// [(p0, e0), ..., (pn, en)]
pub fn factorisation_to_bn_vartime(factorisation: &[(usize, usize)]) -> Vec<u64> {
    let (p0, e0) = factorisation[0];
    let mut n: Vec<u64> = prime_power_to_bn_vartime(p0, e0);
    for (pi, ei) in factorisation.iter().skip(1) {
        let ni = prime_power_to_bn_vartime(*pi, *ei);
        n = mul_bn_vartime(&n, &ni);
    }
    n
}
