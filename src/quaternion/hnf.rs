use crate::quaternion::algebra::{IntQuat, QuatConfig};
use crate::quaternion::hnf_ext::HnfModExt;
use crate::bigint::BigIntAlg;

pub fn xgcd_with_u_not_0<T: BigIntAlg>(a: &T, b: &T) -> (T, T, T) {
    let (d, mut u, mut v) = a.xgcd(b);
    if u.is_zero() {
        u = u + (b.clone() / d.clone());
        v = v - (a.clone() / d.clone());
    }
    (d, u, v)
}

/// Hermite Normal Form Modulo algorithm from Henri Cohen.
/// Mutates `generators` in-place, returning the reduced basis as 4 IntQuats.
pub fn quat_hnf_mod_core<T: BigIntAlg, P: QuatConfig<T>>(
    generators: &mut [IntQuat<T, P>],
    modulo: &T,
) -> [IntQuat<T, P>; 4] {
    let n = generators.len();
    debug_assert!(n >= 4, "HNF modulo requires at least 4 generators");

    let mut w = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];
    let mut m = modulo.clone().abs();

    let mut i: isize = 3;
    let mut k: isize = (n as isize) - 1;
    let mut j: isize = (n as isize) - 1;

    // TODO: Since we operate column by column (i from 3 down to 0),
    // coordinates > i are already zero. Full IntQuat scalar multiplications here 
    // compute on zeroes unnecessarily.
    while i >= 0 {
        let ui = i as usize;

        while j > 0 {
            j -= 1;
            let uj = j as usize;
            let uk = k as usize;

            if !generators[uj].coords[ui].is_zero() {
                let (d, u, v) = xgcd_with_u_not_0(&generators[uk].coords[ui], &generators[uj].coords[ui]);

                let c = &generators[uk] * &u + &generators[uj] * &v;
                let coeff_1 = generators[uk].coords[ui].clone() / d.clone();
                let coeff_2 = -(generators[uj].coords[ui].clone() / d.clone());

                generators[uj] = &generators[uj] * &coeff_1 + &generators[uk] * &coeff_2;
                generators[uj].centered_mod_mut(&m);

                generators[uk] = c;
                generators[uk].centered_mod_mut(&m);
            }
        }

        let uk = k as usize;
        let (d, u, _v) = xgcd_with_u_not_0(&generators[uk].coords[ui], &m);

        for x in 0..4 {
            w[ui].coords[x] = (u.clone() * generators[uk].coords[x].clone()).positive_mod(&m);
        }

        if w[ui].coords[ui].is_zero() {
            w[ui].coords[ui] = m.clone();
        }

        for h in (ui + 1)..4 {
            let mut q = w[h].coords[ui].clone() / w[ui].coords[ui].clone();
            let rem = w[h].coords[ui].clone() % w[ui].coords[ui].clone();
            if rem < T::zero() {
                q = q - T::from_i32(1); // floor adjustment
            }
            q = -q;
            w[h] = &w[h] + &w[ui] * &q;
        }

        m = m / d;

        if i != 0 {
            k -= 1;
            i -= 1;
            j = k;
            if generators[k as usize].coords[i as usize].is_zero() {
                generators[k as usize].coords[i as usize] = m.clone();
            }
        } else {
            k -= 1;
            i -= 1;
            j = k;
        }
    }

    w
}
