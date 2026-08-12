use isogeny::quaternion::dim2::dim2_lattice_short_basis;
use isogeny::quaternion::dim2::Mat2x2;
use isogeny::quaternion::lattice::IntLattice;

use isogeny::bigint::BigIntAlg;
use isogeny::quaternion::algebra::{IntQuat, QuatConfig, RatQuat};
use isogeny::quaternion::dim2::Vec2;
use isogeny::quaternion::ideal::QuatLeftIdeal;
use isogeny::quaternion::lattice::RatLattice;
use isogeny::quaternion::qlapoti::{
    quat_dim2_lattice_qlapoti_cvp_condition, quat_lideal_generator_small_coprime,
    quat_lideal_shortest_equivalent, quat_qlapoti, QlapotiEnumParams,
};
use num_bigint::BigInt;
use std::str::FromStr;
use std::sync::LazyLock;

// ========================================================================
// Test Configurations
// ========================================================================

#[derive(Debug, Clone)]
pub struct P103;
static P103_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(103));
impl QuatConfig<BigInt> for P103 {
    fn p() -> &'static BigInt {
        &P103_VAL
    }
}

#[derive(Debug, Clone)]
pub struct PSQISign;
static PSQISIGN_VAL: LazyLock<BigInt> = LazyLock::new(|| {
    let mut p = BigInt::from(1);
    p = p << 248;
    p = p * BigInt::from(5);
    p - BigInt::from(1)
});
impl QuatConfig<BigInt> for PSQISign {
    fn p() -> &'static BigInt {
        &PSQISIGN_VAL
    }
}

#[inline]
fn b_zero() -> BigInt {
    <BigInt as BigIntAlg>::zero()
}
#[inline]
fn b_one() -> BigInt {
    <BigInt as BigIntAlg>::one()
}
#[inline]
fn b(val: i32) -> BigInt {
    <BigInt as BigIntAlg>::from_i32(val)
}

fn o0_lattice<P: QuatConfig<BigInt>>() -> RatLattice<BigInt, P> {
    let mut lat = RatLattice::zero();
    lat.basis.generators = [
        IntQuat::new_i32(2, 0, 0, 0),
        IntQuat::new_i32(0, 2, 0, 0),
        IntQuat::new_i32(0, 1, 1, 0),
        IntQuat::new_i32(1, 0, 0, 1),
    ];
    lat.denom = b(2);
    lat.hnf();
    lat
}

// ========================================================================
// Qlapoti Tests
// ========================================================================

#[test]
fn test_qlapoti_lideal_short_equivalent() {
    let alg_order = o0_lattice::<P103>();

    let elem_int = IntQuat::new_i32(1, 0, -2, 1);
    let elem = RatQuat::new(elem_int, b(1));
    let n = b(3);

    let mut lideal = QuatLeftIdeal::create(&alg_order, &elem, &n);

    let mut equiv = QuatLeftIdeal {
        lattice: RatLattice::zero(),
        norm: b_zero(),
        parent_order: alg_order.clone().into(),
    };

    let mut elem_out = RatQuat::new(IntQuat::zero(), b_one());

    let elem2_int = IntQuat::new_i32(200, 5192, 2397665, -12233);
    let elem2 = RatQuat::new(elem2_int, b(1));
    lideal = QuatLeftIdeal::mul(&lideal, &elem2);

    quat_lideal_shortest_equivalent(&mut equiv, &mut elem_out, &lideal);

    let (_, norm_d) = elem_out.norm();
    assert_eq!(norm_d, b_one(), "Denominator of norm should be 1");
}

#[test]
fn test_qlapoti_lideal_generator_small_coprime() {
    let alg_order = o0_lattice::<P103>();

    let gen_int = IntQuat::new_i32(3, 5, 7, 11);
    let mut gene = RatQuat::new(gen_int, b(1));
    let n = b(51); 

    let lideal = QuatLeftIdeal::create(&alg_order, &gene, &n);

    quat_lideal_generator_small_coprime(&mut gene, &lideal, 16);

    let (_norm_n, norm_d) = gene.norm();
    assert_eq!(norm_d, b_one(), "Generator norm denominator must be 1");

    let lideal2 = QuatLeftIdeal::create(&alg_order, &gene, &n);
    assert_eq!(
        RatLattice::equal(&lideal.lattice, &lideal2.lattice), 
        u32::MAX, 
        "Lattices are not equal modulo N"
    );

    let gen2_int = IntQuat::new_i32(2, 4, 2, 6);
    let mut gene2 = RatQuat::new(gen2_int, b(1));
    let n2 = b(23);

    let lideal3 = QuatLeftIdeal::create(&alg_order, &gene2, &n2);

    quat_lideal_generator_small_coprime(&mut gene2, &lideal3, 16);

    let (_, norm_d2) = gene2.norm();
    assert_eq!(norm_d2, b_one());

    let lideal4 = QuatLeftIdeal::create(&alg_order, &gene2, &n2);
    assert_eq!(
        RatLattice::equal(&lideal3.lattice, &lideal4.lattice), 
        u32::MAX, 
        "Second ideal equality check failed"
    );
}

#[test]
fn test_lattice_qlapoti_cvp_condition() {
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();
    let m = BigInt::from_str("113078212145816597093331034492692987310554024929624970589366518801641159342").unwrap();

    let a_alpha = BigInt::from_str("-55802309567322429328887817479381622545473").unwrap();
    let b_alpha = BigInt::from_str("164832325035431386247229290667677491931168").unwrap();

    let q_params = QlapotiEnumParams {
        a_alpha: &a_alpha,
        b_alpha: &b_alpha,
        m: &m,
        n: &n,
    };

    let ab = Vec2::new(
        BigInt::from_str("-3324930064024875210").unwrap(),
        BigInt::from_str("-2166920436366097869").unwrap(),
    );

    let mut out = RatQuat::<BigInt, P103>::new(IntQuat::zero(), b_one());
    let found = quat_dim2_lattice_qlapoti_cvp_condition(&mut out, &ab, &q_params);

    assert!(found, "CVP Condition failed to find solution");

    let cmp0 = BigInt::from_str("5098956283653089832").unwrap();
    let cmp1 = BigInt::from_str("-1774026219628214622").unwrap();
    let cmp2 = BigInt::from_str("1185358405318603462").unwrap();
    let cmp3 = BigInt::from_str("981562031047494407").unwrap();

    assert_eq!(out.num.coords[0], cmp0, "Coord 0 mismatch");
    assert_eq!(out.num.coords[1], cmp1, "Coord 1 mismatch");
    assert_eq!(out.num.coords[2], cmp2, "Coord 2 mismatch");
    assert_eq!(out.num.coords[3], cmp3, "Coord 3 mismatch");

    let ab2 = Vec2::new(
        BigInt::from_str("6588358983801282790").unwrap(),
        BigInt::from_str("4139758051024177521").unwrap(),
    );

    let mut out2 = RatQuat::<BigInt, P103>::new(IntQuat::zero(), b_one());
    let found2 = quat_dim2_lattice_qlapoti_cvp_condition(&mut out2, &ab2, &q_params);

    assert!(found2, "Second CVP Condition failed to find solution");

    let cmp2_0 = BigInt::from_str("-2743025951076238849").unwrap();
    let cmp2_1 = BigInt::from_str("-3845333032725043941").unwrap();
    let cmp2_2 = BigInt::from_str("-1506191648044069223").unwrap();
    let cmp2_3 = BigInt::from_str("-2633566402980108298").unwrap();

    assert_eq!(out2.num.coords[0], cmp2_0, "Second CVP Coord 0 mismatch");
    assert_eq!(out2.num.coords[1], cmp2_1, "Second CVP Coord 1 mismatch");
    assert_eq!(out2.num.coords[2], cmp2_2, "Second CVP Coord 2 mismatch");
    assert_eq!(out2.num.coords[3], cmp2_3, "Second CVP Coord 3 mismatch");
}

#[test]
fn test_qlapoti_qlapoti() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let c2 = b_zero();
    let c3 = b_one();

    let elem = RatQuat::new(IntQuat::new(c0, c1, c2, c3), b(2));
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let lideal = QuatLeftIdeal::create(&alg_order, &elem, &n);
    assert_eq!(lideal.norm, n, "Initial ideal norm is incorrect");

    let mut mu1 = RatQuat::<BigInt, PSQISign>::new(IntQuat::zero(), b_one());
    let mut mu2 = RatQuat::<BigInt, PSQISign>::new(IntQuat::zero(), b_one());
    let mut theta = RatQuat::<BigInt, PSQISign>::new(IntQuat::zero(), b_one());
    let mut smallest = RatQuat::<BigInt, PSQISign>::new(IntQuat::zero(), b_one());

    let output = quat_qlapoti(
        &mut mu1,
        &mut mu2,
        &mut theta,
        &mut smallest,
        &lideal,
        20000,
        17,
        246,
    );

    assert!(output != 0, "Qlapoti algorithm failed to find a valid solution");

    // `smallest` must belong to the original ideal lattice
    assert!(
        lideal.lattice.contains(&smallest).is_some(),
        "Smallest equivalent generator is not contained in the original ideal"
    );

    // small_equiv = lideal * (conj(smallest) / N(lideal))
    let mut smallest_conj = smallest.clone();
    smallest_conj.num = smallest_conj.num.conj();
    smallest_conj.denom = smallest_conj.denom * lideal.norm.clone();
    let small_equiv = QuatLeftIdeal::mul(&lideal, &smallest_conj);

    // mu1 and mu2 belong to small_equiv
    assert!(
        small_equiv.lattice.contains(&mu1).is_some(),
        "mu1 is not contained in the equivalent ideal"
    );
    assert!(
        small_equiv.lattice.contains(&mu2).is_some(),
        "mu2 is not contained in the equivalent ideal"
    );

    // N(small_equiv) < sqrt(p)
    let p_sqrt = PSQISign::p().sqrt();
    assert!(
        small_equiv.norm < p_sqrt,
        "Norm of small_equiv exceeds sqrt(p)"
    );

    // I1 = small_equiv * (conj(mu1) / N(small_equiv))
    let mut mu1_conj = mu1.clone();
    mu1_conj.num = mu1_conj.num.conj();
    mu1_conj.denom = mu1_conj.denom * small_equiv.norm.clone();
    let i1 = QuatLeftIdeal::mul(&small_equiv, &mu1_conj);

    // theta == mu2 * (conj(mu1) / N(small_equiv))
    let cmp = mu2.clone() * mu1_conj.clone();
    for i in 0..4 {
        let lhs = cmp.num.coords[i].clone() * theta.denom.clone();
        let rhs = theta.num.coords[i].clone() * cmp.denom.clone();
        assert_eq!(
            lhs, rhs,
            "Theta mathematical validation failed at coordinate {}", i
        );
    }

    // I2 = small_equiv * (conj(mu2) / N(small_equiv))
    let mut mu2_conj = mu2.clone();
    mu2_conj.num = mu2_conj.num.conj();
    mu2_conj.denom = mu2_conj.denom * small_equiv.norm.clone();
    let i2 = QuatLeftIdeal::mul(&small_equiv, &mu2_conj);

    // N(I1) + N(I2) == 2^(246 - transformed_flag)
    let sum_norms = i1.norm.clone() + i2.norm.clone();
    let expected_power = 246 - if (output & 2) != 0 { 1 } else { 0 };
    let expected_sum = b_one() << expected_power;

    assert_eq!(
        sum_norms, expected_sum,
        "KLPT equation failed! N(I1) + N(I2) != 2^e"
    );

}


#[test]
fn test_qlapoti_intermediate_ideal_norm() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0, c1, b_zero(), b_one()), b(2));

    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let lideal = QuatLeftIdeal::create(&alg_order, &elem, &n);

    assert_eq!(lideal.norm, n, "The initial left ideal norm MUST equal n.");
}

#[test]
fn test_qlapoti_intermediate_shortest_equivalent() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0, c1, b_zero(), b_one()), b(2));
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let lideal = QuatLeftIdeal::create(&alg_order, &elem, &n);

    let mut smallest = RatQuat::<BigInt, PSQISign>::new(IntQuat::zero(), b_one());
    let mut small = QuatLeftIdeal {
        lattice: RatLattice::zero(),
        norm: b_zero(),
        parent_order: lideal.parent_order.clone(),
    };
    quat_lideal_shortest_equivalent(&mut small, &mut smallest, &lideal);

    assert_eq!(small.norm, n, "The equivalent ideal mathematically bounds to N(I).");
}

#[test]
fn test_qlapoti_intermediate_alpha0_norm_bound() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0, c1, b_zero(), b_one()), b(2));
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let lideal = QuatLeftIdeal::create(&alg_order, &elem, &n);

    let mut smallest = RatQuat::<BigInt, PSQISign>::new(IntQuat::zero(), b_one());
    let mut small = QuatLeftIdeal {
        lattice: RatLattice::zero(),
        norm: b_zero(),
        parent_order: lideal.parent_order.clone(),
    };
    quat_lideal_shortest_equivalent(&mut small, &mut smallest, &lideal);

    let mut gram = [
        [b_zero(), b_zero(), b_zero(), b_zero()],
        [b_zero(), b_zero(), b_zero(), b_zero()],
        [b_zero(), b_zero(), b_zero(), b_zero()],
        [b_zero(), b_zero(), b_zero(), b_zero()]
    ];
    let mut red_basis = IntLattice::zero();
    isogeny::quaternion::lll::quat_lideal_reduce_basis(&mut red_basis, &mut gram, &small);
    small.lattice.basis = red_basis;

    let mut alpha_0 = RatQuat::new(IntQuat::zero(), b_one());
    quat_lideal_generator_small_coprime(&mut alpha_0, &small, 17);

    let (norm_n, norm_d) = alpha_0.norm();
    let actual_norm = norm_n / norm_d;

    let sub_term = (actual_norm * b(2)) / n;
    let target_two_e = b_one() << 246;

    assert!(sub_term < target_two_e, "KLPT sub_term exponentially exceeds the bounds, pushing m into negatives.");
}

#[test]
fn test_lattice_determinant_and_scaling() {
    let lat = o0_lattice::<PSQISign>();

    let det_o0 = lat.basis.inv_with_det(None);
    assert_eq!(det_o0, b(4), "Determinant of O0 basis should be 4");

    let n = b(10);
    let mut lat_scaled = lat.clone();
    for i in 0..4 {
        lat_scaled.basis.generators[i] = &lat_scaled.basis.generators[i] * &n;
    }

    let det_scaled = lat_scaled.basis.inv_with_det(None);
    assert_eq!(det_scaled, b(40000), "Scalar multiplication is bugged in IntQuat * T");
}

#[test]
fn test_lattice_principal_ideal_norm() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0.clone(), c1.clone(), b_zero(), b_one()), b(2));

    let principal = QuatLeftIdeal::<BigInt, PSQISign>::create_principal(&alg_order, &elem);

    let (norm_n, norm_d) = elem.norm();
    assert_eq!(norm_d, b_one(), "Element norm must be integer");

    let computed_index = RatLattice::index(&principal.lattice, &alg_order);
    let computed_norm = computed_index.sqrt();

    assert_eq!(computed_norm, norm_n, "RatLattice::index or alg_elem_mul is corrupting the principal ideal geometry");
}

#[test]
fn test_lattice_add_lazy_gcd_modulo() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0, c1, b_zero(), b_one()), b(2));
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let principal = QuatLeftIdeal::create_principal(&alg_order, &elem);

    let mut order_n = alg_order.clone();
    for i in 0..4 {
        order_n.basis.generators[i] = &order_n.basis.generators[i] * &n;
    }

    let mut scaled1 = principal.lattice.basis.clone();
    scaled1.scale(&order_n.denom);
    let mut scaled2 = order_n.basis.clone();
    scaled2.scale(&principal.lattice.denom);

    let det1 = scaled1.inv_with_det(None);
    let det2 = scaled2.inv_with_det(None);
    let mod_val = det1.gcd(&det2);

    let n_sq = n.clone() * n.clone();
    assert!(mod_val > n_sq, "The HNF modulo evaluation is not correct");
}

#[test]
fn test_lattice_hnf_modulo_reduction_integrity() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0, c1, b_zero(), b_one()), b(2));
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let principal = QuatLeftIdeal::create_principal(&alg_order, &elem);

    let mut order_n = alg_order.clone();
    for i in 0..4 {
        order_n.basis.generators[i] = &order_n.basis.generators[i] * &n;
    }

    let sum_lattice = RatLattice::add_lazy(&principal.lattice, &order_n);

    let expected_det = b(4) * n.clone() * n.clone();

    let mut inv_mat = isogeny::quaternion::lattice::IntLattice::zero();
    let computed_det = sum_lattice.basis.inv_with_det(Some(&mut inv_mat)).abs();

    assert_eq!(
        computed_det, expected_det,
        "HNF reduction failed, The resulting determinant is incorrect."
    );
}

#[test]
fn test_lattice_add_lazy_commutativity() {
    let alg_order = o0_lattice::<PSQISign>();

    let c0 = BigInt::from_str("11679558057699548966295664199498600969").unwrap();
    let c1 = BigInt::from_str("5995913694671076569732436307214390300").unwrap();
    let elem = RatQuat::new(IntQuat::new(c0, c1, b_zero(), b_one()), b(2));
    let n = BigInt::from_str("3588068757273373623184041115224326680").unwrap();

    let principal = QuatLeftIdeal::create_principal(&alg_order, &elem);

    let mut order_n = alg_order.clone();
    for i in 0..4 {
        order_n.basis.generators[i] = &order_n.basis.generators[i] * &n;
    }

    let sum_ab = RatLattice::add_lazy(&principal.lattice, &order_n);
    let sum_ba = RatLattice::add_lazy(&order_n, &principal.lattice);

    assert_eq!(
        RatLattice::equal(&sum_ab, &sum_ba), u32::MAX,
        "Lattice addition is not commutative?"
    );
}



#[test]
fn test_qlapoti_dim2_mat2x2_eval_transposition() {
    let mut mat = Mat2x2::zero();
    mat.m[0][0] = b(1); mat.m[0][1] = b(2);
    mat.m[1][0] = b(3); mat.m[1][1] = b(4);

    let v = Vec2::new(b(10), b(20));

    let res = mat.eval(&v);

    assert_eq!(res.coords[0], b(50), "Matrix evaluation failed on coord 0.");
    assert_eq!(res.coords[1], b(110), "Matrix evaluation failed on coord 1.");
}

#[test]
fn test_qlapoti_dim2_inv_with_det() {
    let mut mat = Mat2x2::zero();
    mat.m[0][0] = b(1); mat.m[0][1] = b(2);
    mat.m[1][0] = b(3); mat.m[1][1] = b(4);

    let (inv, det) = mat.inv_with_det_as_denom();

    assert_eq!(det, b(-2), "Determinant computation is incorrect.");

    assert_eq!(inv.m[0][0], b(4), "Adjugate [0][0] is wrong");
    assert_eq!(inv.m[0][1], b(-2), "Adjugate [0][1] is wrong");
    assert_eq!(inv.m[1][0], b(-3), "Adjugate [1][0] is wrong");
    assert_eq!(inv.m[1][1], b(1), "Adjugate [1][1] is wrong");
}

#[test]
fn test_qlapoti_dim2_short_basis() {
    let mut mat = Mat2x2::zero();
    let n = BigInt::from_str("1000000").unwrap();
    let x = BigInt::from_str("12345").unwrap();

    mat.m[0][0] = n.clone() - x;
    mat.m[0][1] = n.clone();
    mat.m[1][0] = b(1);
    mat.m[1][1] = b(0);

    let red = dim2_lattice_short_basis(&mat, &b(1));

    let (_, det_orig) = mat.inv_with_det_as_denom();
    let (_, det_red) = red.inv_with_det_as_denom();
    assert_eq!(det_orig.abs(), det_red.abs(), "Determinant wrong");

    // ||v1|| * ||v2|| <= 2 * det^2
    let v0_norm = red.m[0][0].clone() * red.m[0][0].clone() + red.m[1][0].clone() * red.m[1][0].clone();
    let v1_norm = red.m[0][1].clone() * red.m[0][1].clone() + red.m[1][1].clone() * red.m[1][1].clone();

    let prod = v0_norm * v1_norm;
    let bound = det_orig.clone() * det_orig.clone() * b(2);

    assert!(prod < bound, "basis is not reduced");
}
