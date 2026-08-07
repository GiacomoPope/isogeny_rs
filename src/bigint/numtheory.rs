use crate::bigint::BigIntAlg;

/// Computes base^exp mod m
pub fn powmod<T: BigIntAlg>(base: &T, exp: &T, m: &T) -> T {
    let mut res = T::one();
    let mut b = base.clone() % m.clone();
    if b < T::zero() {
        b = b + m.clone();
    }
    let mut e = exp.clone();
    let two = T::from_i32(2);

    while !e.is_zero() {
        if !(e.clone() % two.clone()).is_zero() {
            res = (res * b.clone()) % m.clone();
        }
        b = (b.clone() * b.clone()) % m.clone();
        e = e / two.clone();
    }
    res
}

/// Computes the Legendre symbol (a / p) using Euler's criterion.
pub fn legendre<T: BigIntAlg>(a: &T, p: &T) -> i32 {
    let a_mod = a.clone() % p.clone();
    if a_mod.is_zero() {
        return 0;
    }
    let exp = (p.clone() - T::one()) / T::from_i32(2);
    let res = powmod(&a_mod, &exp, p);

    if res.ct_eq(&T::one()) == u32::MAX {
        1
    } else if res.ct_eq(&(p.clone() - T::one())) == u32::MAX {
        -1
    } else {
        0
    }
}

/// Computes the modular square root of n modulo p using Tonelli-Shanks.
/// Requires p to be prime. Returns None if n is a quadratic non-residue.
pub fn sqrt_mod_p<T: BigIntAlg>(n: &T, p: &T) -> Option<T> {
    let n_mod = n.clone() % p.clone();
    if n_mod.is_zero() {
        return Some(T::zero());
    }
    if p.ct_eq(&T::from_i32(2)) == u32::MAX {
        return Some(n_mod);
    }
    if legendre(n, p) != 1 {
        return None;
    }

    let mut q = p.clone() - T::one();
    let mut s = 0_u32;
    let two = T::from_i32(2);

    while (q.clone() % two.clone()).is_zero() {
        q = q / two.clone();
        s += 1;
    }

    let mut z = T::from_i32(2);
    while legendre(&z, p) != -1 {
        z = z + T::one();
    }

    let mut m = s;
    let mut c = powmod(&z, &q, p);
    let mut t = powmod(n, &q, p);
    let mut r = powmod(n, &((q + T::one()) / two.clone()), p);

    loop {
        if t.ct_eq(&T::zero()) == u32::MAX {
            return Some(T::zero());
        }
        if t.ct_eq(&T::one()) == u32::MAX {
            return Some(r);
        }

        let mut t2 = t.clone();
        let mut i = 0_u32;
        while i < m {
            if t2.ct_eq(&T::one()) == u32::MAX {
                break;
            }
            t2 = powmod(&t2, &two, p);
            i += 1;
        }

        if i == m {
            return None;
        }

        let mut b = c.clone();
        for _ in 0..(m - i - 1) {
            b = powmod(&b, &two, p);
        }

        m = i;
        c = powmod(&b, &two, p);
        t = (t * c.clone()) % p.clone();
        r = (r * b) % p.clone();
    }
}

/// Solves the Diophantine equation x^2 + d*y^2 = m using Cornacchia's Algorithm.
/// Returns Some((x, y)) if a solution exists, otherwise None.
pub fn cornacchia<T: BigIntAlg>(d: &T, m: &T) -> Option<(T, T)> {
    let two = T::from_i32(2);

    let minus_d = (m.clone() - (d.clone() % m.clone())) % m.clone();
    let mut r0 = sqrt_mod_p(&minus_d, m)?;

    if r0 > (m.clone() / two.clone()) {
        r0 = m.clone() - r0;
    }

    let mut a = m.clone();
    let mut b = r0;
    let limit = m.sqrt();

    while b > limit {

        let r = a.clone() % b.clone();
        a = b;
        b = r;
    }

    let b_sq = b.clone() * b.clone();
    let diff = m.clone() - b_sq;

    if !(diff.clone() % d.clone()).is_zero() {
        return None;
    }

    let c = diff / d.clone();
    let c_sqrt = c.sqrt();

    // check if c is a perfect square
    if c_sqrt.clone() * c_sqrt.clone() == c {
        Some((b, c_sqrt))
    } else {
        None
    }
}
