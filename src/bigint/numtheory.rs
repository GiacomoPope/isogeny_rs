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

/// Modular square root of n modulo p using Tonelli-Shanks.
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

/// Miller-Rabin primality test
pub fn is_probable_prime<T: BigIntAlg>(n: &T, reps: u32) -> bool {
    let two = T::from_i32(2);
    let three = T::from_i32(3);

    if n.ct_eq(&two) == u32::MAX || n.ct_eq(&three) == u32::MAX { return true; }
    if n < &two || (n.clone() % two.clone()).is_zero() { return false; }

    let mut d = n.clone() - T::one();
    let mut s = 0;
    while (d.clone() % two.clone()).is_zero() {
        d = d / two.clone();
        s += 1;
    }

    let mut a = T::from_i32(2);
    let n_minus_one = n.clone() - T::one();
    let n_minus_two = n.clone() - two.clone();

    for _ in 0..reps {
        let mut x = powmod(&a, &d, n);

        if x.ct_eq(&T::one()) == u32::MAX || x.ct_eq(&n_minus_one) == u32::MAX {
            a = a + T::one();
            if a > n_minus_two {
                a = T::from_i32(2);
            }
            continue;
        }

    let mut composite = true;
        for _ in 1..s {
            x = powmod(&x, &two, n);
            if x.ct_eq(&n_minus_one) == u32::MAX {
                composite = false;
                break;
            }
        }
        if composite { return false; }

        a = a + T::one();
        if a > n_minus_two {
            a = T::from_i32(2);
        }
    }
    true
}

/// Solves the Diophantine equation x^2 + d*y^2 = m using Cornacchia's Algorithm.
pub fn cornacchia<T: BigIntAlg>(d: &T, m: &T) -> Option<(T, T)> {
    if m <= &T::zero() {
        return None;
    }

    // TODO: not sure how well modular square roots behave when m is not prime
    //       I forgot if we can keep them, I need to check this later
    if !is_probable_prime(m, 40) {
        return None;
    }

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

    if c_sqrt.clone() * c_sqrt.clone() == c {
        Some((b, c_sqrt))
    } else {
        None
    }
}

pub fn invmod<T: BigIntAlg>(a: &T, m: &T) -> Option<T> {
    let (gcd, x, _y) = a.xgcd(m);

    if gcd.ct_eq(&T::one()) == u32::MAX {
        let mut inv = x % m.clone();
        if inv < T::zero() {
            inv = inv + m.clone();
        }
        Some(inv)
    } else {
        None
    }
}
