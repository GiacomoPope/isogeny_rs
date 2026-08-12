use crate::bigint::numtheory::{cornacchia, invmod};
use crate::bigint::BigIntAlg;
use crate::quaternion::algebra::{IntQuat, QuatConfig, RatQuat};
use crate::quaternion::dim2::{
    dim2_lattice_short_basis, rounded_div, Mat2x2, Vec2,
};
use crate::quaternion::ideal::QuatLeftIdeal;
use crate::quaternion::lll::quat_lideal_reduce_basis;
use rand::Rng;

pub struct QlapotiEnumParams<'a, T: BigIntAlg> {
    pub a_alpha: &'a T,
    pub b_alpha: &'a T,
    pub m: &'a T,
    pub n: &'a T,
}

pub fn quat_lideal_shortest_equivalent<T: BigIntAlg, P: QuatConfig<T>>(
    equiv: &mut QuatLeftIdeal<T, P>,
    elem: &mut RatQuat<T, P>,
    lideal: &QuatLeftIdeal<T, P>,
) {
    let mut gram = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];
    let mut red = crate::quaternion::lattice::IntLattice::zero();

    quat_lideal_reduce_basis(&mut red, &mut gram, lideal);

    let new_alpha_coord = [T::one(), T::zero(), T::zero(), T::zero()];
    let mut eval_coord = [T::zero(), T::zero(), T::zero(), T::zero()];

    for i in 0..4 {
        let mut sum = T::zero();
        for j in 0..4 {
            sum = sum + red.generators[j].coords[i].clone() * new_alpha_coord[j].clone();
        }
        eval_coord[i] = sum;
    }

    let mut new_alpha = RatQuat::new(
        IntQuat::new(eval_coord[0].clone(), eval_coord[1].clone(), eval_coord[2].clone(), eval_coord[3].clone()),
        lideal.lattice.denom.clone(),
    );

    debug_assert!(
        lideal.lattice.contains(&new_alpha).is_some(),
        "Shortest equivalent generator must belong to the ideal lattice"
    );

    *elem = new_alpha.clone();
    new_alpha.denom = new_alpha.denom * lideal.norm.clone();
    new_alpha.num.coords[0] = -new_alpha.num.coords[0].clone();

    *equiv = QuatLeftIdeal::mul(lideal, &new_alpha);
}

pub fn quat_lideal_generator_small_coprime<T: BigIntAlg, P: QuatConfig<T>>(
    gene: &mut RatQuat<T, P>,
    lideal: &QuatLeftIdeal<T, P>,
    sampling_bound_bits: u32,
) {
    let mut found = false;
    gene.denom = lideal.lattice.denom.clone();
    let n2 = lideal.norm.clone() * lideal.norm.clone();
    let mut rng = rand::rng(); 

    while !found {
        let mut coeffs = [T::zero(), T::zero(), T::zero(), T::zero()];
        for i in 0..4 {
            let limit = 1_i64 << sampling_bound_bits;
            let val = rng.random_range(-limit..=limit); 
            coeffs[i] = T::from_i32(val as i32); 
        }

        for i in 0..4 {
            let mut sum = T::zero();
            for j in 0..4 {
                sum = sum + lideal.lattice.basis.generators[j].coords[i].clone() * coeffs[j].clone();
            }
            gene.num.coords[i] = sum;
        }

        if gene.denom.ct_eq(&T::one()) == u32::MAX {
            let gcd_val = (gene.num.coords[0].clone() * T::from_i32(2)).gcd(&lideal.norm);
            found = gcd_val.ct_eq(&T::one()) == u32::MAX;
        } else {
            found = true;
            if gene.denom.ct_eq(&T::from_i32(2)) == u32::MAX {
                let gcd_val = lideal.norm.gcd(&gene.num.coords[0]);
                found = gcd_val.ct_eq(&T::one()) == u32::MAX;
            }
        }

        if found {
            let (n_val, d_val) = gene.norm();
            let actual_norm = n_val / d_val;
            let gcd_val = n2.gcd(&actual_norm);
            found = gcd_val.ct_eq(&lideal.norm) == u32::MAX;
        }
    }
}

pub fn quat_qlapoti_check_mod_condition<T: BigIntAlg>(m: &T, a: &T, b: &T) -> bool {
    let mut ok = true;

    let a_mod2 = (a.clone() % T::from_i32(2)).abs().ct_eq(&T::one()) == u32::MAX;
    let b_mod2 = (b.clone() % T::from_i32(2)).abs().ct_eq(&T::one()) == u32::MAX;

    let mut m_mod8 = m.clone() % T::from_i32(8);
    if m_mod8 < T::zero() { m_mod8 = m_mod8 + T::from_i32(8); }

    let mut m_mod4 = m_mod8.clone() % T::from_i32(4);
    if m_mod4 < T::zero() { m_mod4 = m_mod4 + T::from_i32(4); }

    if m_mod8.is_zero() {
        ok = false;
    }

    if ok {
        if (a_mod2 == b_mod2) && (!a_mod2) {
            if !m_mod4.is_zero() { ok = false; }
        } else {
            if (a_mod2 == b_mod2) && (a_mod2) {
                if m_mod4.ct_eq(&T::from_i32(2)) != u32::MAX { ok = false; }
            } else {
                if m_mod4.ct_eq(&T::one()) != u32::MAX { ok = false; }
            }
        }
    }

    ok
}

pub fn quat_dim2_lattice_qlapoti_cvp_condition<T: BigIntAlg, P: QuatConfig<T>>(
    elem: &mut RatQuat<T, P>,
    vec: &Vec2<T>,
    q_params: &QlapotiEnumParams<T>,
) -> bool {
    let mut found = true;
    let a_val = -vec.coords[0].clone();
    let b_val = -vec.coords[1].clone();

    let tmp = a_val.clone() * q_params.a_alpha.clone() + b_val.clone() * q_params.b_alpha.clone();
    let mut m2 = q_params.m.clone() - tmp.clone();
    m2 = m2 / q_params.n.clone();

    let a_sq = a_val.clone() * a_val.clone();
    let b_sq = b_val.clone() * b_val.clone();
    m2 = m2 - a_sq.clone() - b_sq.clone();

    m2 = m2.clone() * T::from_i32(2) + a_sq.clone() + b_sq.clone();

    if found {
        found = quat_qlapoti_check_mod_condition(&m2, &a_val, &b_val);
    }

    if found {
        if let Some((x_res, y_res)) = cornacchia(&T::one(), &m2) {
            let mut a0 = x_res;
            let mut b0 = y_res;

            let a_val_parity = (a_val.clone() % T::from_i32(2)).abs().ct_eq(&T::one()) == u32::MAX;
            let a0_parity = (a0.clone() % T::from_i32(2)).abs().ct_eq(&T::one()) == u32::MAX;

            if a_val_parity != a0_parity {
                std::mem::swap(&mut a0, &mut b0);
            }

            let a1 = (a0.clone() + a_val.clone()) / T::from_i32(2);
            let b1 = (b0.clone() + b_val.clone()) / T::from_i32(2);

            let a2 = a_val.clone() - a1.clone();
            let b2 = b_val.clone() - b1.clone();

            elem.num.coords[0] = a1;
            elem.num.coords[1] = a2;
            elem.num.coords[2] = b1;
            elem.num.coords[3] = b2;
        } else {
            found = false;
        }
    }
    found
}

pub fn quat_elem_is_odd_norm<T: BigIntAlg, P: QuatConfig<T>>(elem: &RatQuat<T, P>) -> bool {
    let mut found = false;
    let c0_mod2 = (elem.num.coords[0].clone() % T::from_i32(2)).abs().ct_eq(&T::one()) == u32::MAX;
    let c2_mod2 = (elem.num.coords[2].clone() % T::from_i32(2)).abs().ct_eq(&T::one()) == u32::MAX;

    if c0_mod2 == c2_mod2 {
        return false;
    }

    for i in 0..4 {
        let mut val_mod4 = elem.num.coords[i].clone() % T::from_i32(4);
        if val_mod4 < T::zero() { val_mod4 = val_mod4 + T::from_i32(4); }

        if found {
            if val_mod4.ct_eq(&T::from_i32(2)) == u32::MAX { return false; }
        } else {
            if val_mod4.ct_eq(&T::from_i32(2)) == u32::MAX { found = true; }
        }
    }
    found
}

pub fn get_endtype<T: BigIntAlg, P: QuatConfig<T>>(elem: &RatQuat<T, P>) -> i32 {
    let mut t = [0_i32; 4];
    for i in 0..4 {
        let mut val_mod4 = elem.num.coords[i].clone() % T::from_i32(4);
        if val_mod4 < T::zero() { val_mod4 = val_mod4 + T::from_i32(4); }

        if val_mod4.is_zero() { t[i] = 0; }
        else if val_mod4.ct_eq(&T::one()) == u32::MAX { t[i] = 1; }
        else if val_mod4.ct_eq(&T::from_i32(2)) == u32::MAX { t[i] = 2; }
        else { t[i] = 3; }
    }

    let (t1, t2, t3, t4) = (t[0], t[1], t[3], t[2]);

    if t1 == 2 {
        if t2 == 2 {
            if (t3 == 1 && t4 == 0) || (t3 == 3 && t4 == 2) { return 1; }
            if t3 == 0 && (t4 == 1 || t4 == 3) { return 2; }
            return 0;
        }
        if t2 == 0 && ((t3 == 1 && t4 == 2) || (t3 == 3 && t4 == 2)) { return 1; }
        return 0;
    }
    if t1 == 0 && t2 == 2 && t3 == 2 && (t4 == 1 || t4 == 3) { return 2; }
    0
}

pub fn quat_qlapoti<T: BigIntAlg, P: QuatConfig<T>>(
    mu1: &mut RatQuat<T, P>,
    mu2: &mut RatQuat<T, P>,
    theta: &mut RatQuat<T, P>,
    smallest: &mut RatQuat<T, P>,
    lideal: &QuatLeftIdeal<T, P>,
    max_counter_alpha: i32,
    gen_sampling_bound_bits: u32,
    two_power: u32,
) -> i32 {
    let mut found = false;
    let mut transformed = false;
    let mut keep_alpha = false;

    let mut small = QuatLeftIdeal {
        lattice: crate::quaternion::lattice::RatLattice::zero(),
        norm: T::zero(),
        parent_order: lideal.parent_order.clone(),
    };
    quat_lideal_shortest_equivalent(&mut small, smallest, lideal);

    let mut gram = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];
    let mut red_basis = crate::quaternion::lattice::IntLattice::zero();

    quat_lideal_reduce_basis(&mut red_basis, &mut gram, &small);
    small.lattice.basis = red_basis;

    let n = small.norm.clone();
    let two_e = T::one() << two_power; 
    let psqrt = P::p().sqrt();

    let mut alpha = RatQuat::new(IntQuat::zero(), T::one());
    let mut alpha_0 = RatQuat::new(IntQuat::zero(), T::one());
    let mut lam = T::one();
    let mut alpha0norm = T::zero();

    let mut l_inv = Mat2x2::zero();
    let mut l_det = T::one();
    let mut l_red = Mat2x2::zero();

    for counter in 0..(max_counter_alpha * 20) {
        found = true;

        if counter % max_counter_alpha == 0 {
            keep_alpha = false;
            lam = T::one();
        }

        if !keep_alpha {
            quat_lideal_generator_small_coprime(&mut alpha_0, &small, gen_sampling_bound_bits);
            alpha = alpha_0.clone();

            let (n_val, d_val) = alpha_0.norm();
            alpha0norm = n_val / d_val;

            lam = T::one();
        } else {
            lam = lam + T::one();
            alpha = alpha + alpha_0.clone();
            let mut gcd_check = small.norm.gcd(&lam);
            while gcd_check.ct_eq(&T::one()) != u32::MAX {
                lam = lam + T::one();
                alpha = alpha + alpha_0.clone();
                gcd_check = small.norm.gcd(&lam);
            }
        }

        let norm_n_val = alpha0norm.clone() * T::from_i32(2);
        let norm_d_val = lam.clone() * lam.clone();
        let sub_term = (norm_n_val * norm_d_val) / n.clone();
        let m = two_e.clone() - sub_term.clone();

        if m <= T::zero() {
            found = false;
            continue;
        }

        let mut a_alpha = alpha.num.coords[0].clone();
        let mut b_alpha = alpha.num.coords[1].clone();
        if alpha.denom.ct_eq(&T::one()) == u32::MAX {
            a_alpha = a_alpha.clone() * T::from_i32(2);
            b_alpha = b_alpha.clone() * T::from_i32(2);
        }

        let tmp = match invmod(&a_alpha, &n) {
            Some(inv) => inv,
            None => {
                continue;
            }
        };

        let mut t_val = (m.clone() * tmp.clone()) % n.clone();
        if t_val < T::zero() {
            t_val = t_val + n.clone();
        }

        let v_target = Vec2::new(-t_val.clone(), T::zero());

        if !keep_alpha {
            let mut x_val = (b_alpha.clone() * tmp.clone()) % n.clone();
            if x_val < T::zero() {
                x_val = x_val + n.clone();
            }

            let mut l_mat = Mat2x2::zero();
            l_mat.m[0][0] = n.clone() - x_val;
            l_mat.m[0][1] = n.clone();
            l_mat.m[1][0] = T::one();
            l_mat.m[1][1] = T::zero();

            l_red = dim2_lattice_short_basis(&l_mat, &T::one());
            let (inv, det) = l_red.inv_with_det_as_denom();
            l_inv = inv;
            l_det = det;

            let gcd_val = a_alpha.gcd(&b_alpha);
            keep_alpha = gcd_val.ct_eq(&T::one()) == u32::MAX;

            let norm_n_check = l_red.m[1][1].clone() * l_red.m[1][1].clone() 
                             + l_red.m[0][1].clone() * l_red.m[0][1].clone();
            keep_alpha = keep_alpha && (norm_n_check < psqrt);
        }

        let v_close_unrounded = l_inv.eval(&v_target);
        let mut v_close = v_close_unrounded.clone();
        v_close.coords[0] = rounded_div(&v_close.coords[0], &l_det);
        v_close.coords[1] = rounded_div(&v_close.coords[1], &l_det);
        v_close = l_red.eval(&v_close);

        let v_diff = Vec2::new(
            v_target.coords[0].clone() - v_close.coords[0].clone(),
            v_target.coords[1].clone() - v_close.coords[1].clone()
        );

        let q_params = QlapotiEnumParams {
            a_alpha: &a_alpha,
            b_alpha: &b_alpha,
            m: &m,
            n: &n,
        };

        let mut encoding = RatQuat::<T, P>::new(IntQuat::zero(), T::one());
        found = quat_dim2_lattice_qlapoti_cvp_condition(&mut encoding, &v_diff, &q_params);

        if !found { continue; }

        let mut gamma1 = RatQuat::new(
            IntQuat::new(
                encoding.num.coords[0].clone() * n.clone(),
                encoding.num.coords[2].clone() * n.clone(),
                T::zero(), T::zero()
            ), T::one()
        );

        let mut gamma2 = RatQuat::new(
            IntQuat::new(
                encoding.num.coords[1].clone() * n.clone(),
                encoding.num.coords[3].clone() * n.clone(),
                T::zero(), T::zero()
            ), T::one()
        );

        gamma1 = gamma1 + alpha.clone();
        gamma2 = gamma2 + alpha.clone();

        *mu1 = gamma1.clone();
        *mu2 = gamma2.clone();

        let mut gamma1_conj = gamma1.clone();
        gamma1_conj.num = gamma1_conj.num.conj();
        gamma1_conj.denom = gamma1_conj.denom * n.clone();
        *theta = mu2.clone() * gamma1_conj;

        if theta.denom.ct_eq(&T::one()) == u32::MAX {
            let endtype = get_endtype(theta);
            if endtype == 1 {
                let temp = mu1.clone() + mu2.clone();
                *mu1 = mu1.clone() - mu2.clone();
                *mu2 = temp;
            } else if endtype == 2 {
                let temp_alg = RatQuat::new(IntQuat::new(T::one(), T::zero(), T::one(), T::zero()), T::one());
                *mu2 = temp_alg * mu2.clone();
                let temp = mu1.clone() + mu2.clone();
                *mu1 = mu1.clone() - mu2.clone();
                *mu2 = temp;
            } else {
                keep_alpha = false;
                lam = T::one();
                found = false;
                transformed = false;
                continue;
            }

            gamma1_conj = mu1.clone();
            gamma1_conj.num = gamma1_conj.num.conj();
            gamma1_conj.denom = gamma1_conj.denom * n.clone();
            *theta = mu2.clone() * gamma1_conj;
            transformed = true;
        } else {
            if !quat_elem_is_odd_norm(theta) {
                keep_alpha = false;
                lam = T::one();
                found = false;
                transformed = false;
                continue;
            }
        }

        found = true;
        break;
    }

    if found {
        1 + if transformed { 2 } else { 0 }
    } else {
        0
    }
}
