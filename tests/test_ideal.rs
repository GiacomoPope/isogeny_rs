use isogeny::quaternion::algebra::{IntQuat, QuatConfig, RatQuat};
use isogeny::bigint::BigIntAlg;
use isogeny::quaternion::ideal::{quat_order_discriminant, quat_order_is_maximal, QuatLeftIdeal};
use isogeny::quaternion::lattice::RatLattice;
use num_bigint::BigInt;
use std::sync::LazyLock;

// ========================================================================
// Test Configurations
// ========================================================================

pub struct P23;
static P23_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(23));
impl QuatConfig<BigInt> for P23 {
    fn p() -> &'static BigInt {
        &P23_VAL
    }
}

pub struct P43;
static P43_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(43));
impl QuatConfig<BigInt> for P43 {
    fn p() -> &'static BigInt {
        &P43_VAL
    }
}

pub struct P103;
static P103_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(103));
impl QuatConfig<BigInt> for P103 {
    fn p() -> &'static BigInt {
        &P103_VAL
    }
}

pub struct P367;
static P367_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(367));
impl QuatConfig<BigInt> for P367 {
    fn p() -> &'static BigInt {
        &P367_VAL
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

/// Helper function to create the standard O0 order lattice
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
// Tests
// ========================================================================

#[test]
fn test_lideal_norm() {
    let o0 = o0_lattice::<P23>();
    let gene = RatQuat::<BigInt, P23>::new_i32(1, 3, 2, 1, 4);
    let norm_val = b(17);

    let mut lideal = QuatLeftIdeal::create(&o0, &gene, &norm_val);
    let computed_norm = lideal.norm.clone();

    // Verify norm re-computation via index
    lideal.norm = b_zero();
    lideal.update_norm();
    assert_eq!(lideal.norm, computed_norm);
}

#[test]
fn test_lideal_create_principal() {
    let o0 = o0_lattice::<P367>();

    let mut gamma = RatQuat::<BigInt, P367>::new_i32(1, 219, 200, 78, -1);
    let ideal = QuatLeftIdeal::create_principal(&o0, &gamma);

    assert_eq!(ideal.norm, b(2321156));
    assert_eq!(ideal.lattice.denom, b_one());

    let gens = &ideal.lattice.basis.generators;
    assert_eq!(gens[0].coords[0], b(1160578));
    assert_eq!(gens[0].coords[1], b_zero());
    assert_eq!(gens[0].coords[2], b_zero());
    assert_eq!(gens[0].coords[3], b_zero());

    assert_eq!(gens[1].coords[0], b_zero());
    assert_eq!(gens[1].coords[1], b(1160578));
    assert_eq!(gens[1].coords[2], b_zero());
    assert_eq!(gens[1].coords[3], b_zero());

    assert_eq!(gens[2].coords[0], b(310126));
    assert_eq!(gens[2].coords[1], b(182529));
    assert_eq!(gens[2].coords[2], b_one());
    assert_eq!(gens[2].coords[3], b_zero());

    assert_eq!(gens[3].coords[0], b(978049));
    assert_eq!(gens[3].coords[1], b(310126));
    assert_eq!(gens[3].coords[2], b_zero());
    assert_eq!(gens[3].coords[3], b_one());

    // Same test, but with gamma not reduced
    gamma = RatQuat::<BigInt, P367>::new_i32(2, 438, 400, 156, -2);
    let ideal2 = QuatLeftIdeal::create_principal(&o0, &gamma);
    assert_eq!(ideal2.ct_eq(&ideal), u32::MAX);
}

#[test]
fn test_lideal_create_from_primitive() {
    let o0 = o0_lattice::<P367>();
    let gamma = RatQuat::<BigInt, P367>::new_i32(1, 219, 200, 78, -1);
    let n = b(31);
    let ideal = QuatLeftIdeal::create(&o0, &gamma, &n);

    assert_eq!(ideal.norm, n);
    assert_eq!(ideal.lattice.denom, b(2));

    let gens = &ideal.lattice.basis.generators;
    assert_eq!(gens[0].coords[0], b(62));
    assert_eq!(gens[0].coords[1], b_zero());
    assert_eq!(gens[0].coords[2], b_zero());
    assert_eq!(gens[0].coords[3], b_zero());

    assert_eq!(gens[1].coords[0], b_zero());
    assert_eq!(gens[1].coords[1], b(62));
    assert_eq!(gens[1].coords[2], b_zero());
    assert_eq!(gens[1].coords[3], b_zero());

    assert_eq!(gens[2].coords[0], b(2));
    assert_eq!(gens[2].coords[1], b(1));
    assert_eq!(gens[2].coords[2], b(1));
    assert_eq!(gens[2].coords[3], b_zero());

    assert_eq!(gens[3].coords[0], b(61));
    assert_eq!(gens[3].coords[1], b(2));
    assert_eq!(gens[3].coords[2], b_zero());
    assert_eq!(gens[3].coords[3], b(1));
}

#[test]
fn test_lideal_generator() {
    let o0 = o0_lattice::<P103>();
    let gene = RatQuat::<BigInt, P103>::new_i32(1, 3, 5, 7, 11);
    let n = b(17);

    let lideal = QuatLeftIdeal::create(&o0, &gene, &n);
    let gen_out = lideal.generator().expect("Failed to find a valid generator");

    let lideal2 = QuatLeftIdeal::create(&o0, &gen_out, &n);
    assert_eq!(lideal.ct_eq(&lideal2), u32::MAX);
}

#[test]
fn test_lideal_mul() {
    let o0 = o0_lattice::<P103>();
    let gen1 = RatQuat::<BigInt, P103>::new_i32(1, 3, 5, 7, 11);
    let mut gen2 = RatQuat::<BigInt, P103>::new_i32(1, -2, 13, -17, 19);

    let lideal1 = QuatLeftIdeal::create_principal(&o0, &gen1);
    let prod_ideal = QuatLeftIdeal::mul(&lideal1, &gen2);

    let gen_prod = &gen1 * &gen2;
    let lideal2 = QuatLeftIdeal::create_principal(&o0, &gen_prod);
    assert_eq!(prod_ideal.ct_eq(&lideal2), u32::MAX);

    gen2 = RatQuat::<BigInt, P103>::new_i32(2, -2, 13, -17, 19);
    let prod_ideal_frac = QuatLeftIdeal::mul(&lideal1, &gen2);
    let gen_prod_frac = &gen1 * &gen2;
    let lideal_frac2 = QuatLeftIdeal::create_principal(&o0, &gen_prod_frac);
    assert_eq!(prod_ideal_frac.ct_eq(&lideal_frac2), u32::MAX);
}

#[test]
fn test_lideal_add_intersect_equals() {
    let o0 = o0_lattice::<P103>();

    let gen1 = RatQuat::<BigInt, P103>::new_i32(1, 3, 5, 7, 11);
    let n1 = b(17);
    let lideal1 = QuatLeftIdeal::create(&o0, &gen1, &n1);

    let gen2 = RatQuat::<BigInt, P103>::new_i32(1, -2, 13, -17, 19);
    let n2 = b(43);
    let lideal2 = QuatLeftIdeal::create(&o0, &gen2, &n2);

    let gen3 = &gen2 * &gen1;
    let lideal3 = QuatLeftIdeal::create_principal(&o0, &gen3);

    let mut lideal4 = QuatLeftIdeal::add(&lideal1, &lideal2);
    assert_eq!(lideal4.norm, b_one());
    assert_eq!(RatLattice::equal(&lideal4.lattice, &o0), u32::MAX);

    lideal4 = QuatLeftIdeal::intersect(&lideal1, &lideal1);
    assert_eq!(lideal4.ct_eq(&lideal1), u32::MAX);

    lideal4 = QuatLeftIdeal::add(&lideal1, &lideal1);
    assert_eq!(lideal4.ct_eq(&lideal1), u32::MAX);

    lideal4 = QuatLeftIdeal::add(&lideal1, &lideal3);
    assert_eq!(lideal4.ct_eq(&lideal1), u32::MAX);

    lideal4 = QuatLeftIdeal::intersect(&lideal1, &lideal3);
    assert_eq!(lideal4.ct_eq(&lideal3), u32::MAX);

    lideal4 = QuatLeftIdeal::intersect(&lideal1, &lideal2);
    lideal4 = QuatLeftIdeal::add(&lideal4, &lideal2);
    assert_eq!(lideal4.ct_eq(&lideal2), u32::MAX);

    lideal4 = QuatLeftIdeal::intersect(&lideal1, &lideal2);
    let mut lideal5 = QuatLeftIdeal::intersect(&lideal1, &lideal3);
    lideal4 = QuatLeftIdeal::add(&lideal4, &lideal5);

    lideal5 = QuatLeftIdeal::add(&lideal2, &lideal3);
    lideal5 = QuatLeftIdeal::intersect(&lideal1, &lideal5);

    assert_eq!(lideal4.ct_eq(&lideal5), u32::MAX);
    assert_eq!(lideal4.norm, b(17));
}

#[test]
fn test_lideal_order_discriminant() {
    let o0 = o0_lattice::<P43>();
    let disc = quat_order_discriminant(&o0);
    assert_eq!(disc, b(43));
}

#[test]
fn test_lideal_order_is_maximal() {
    let o0 = o0_lattice::<P43>();
    assert_eq!(quat_order_is_maximal(&o0), u32::MAX);

    let mut id = RatLattice::<BigInt, P43>::zero();
    id.basis.generators = [
        IntQuat::new_i32(1, 0, 0, 0),
        IntQuat::new_i32(0, 1, 0, 0),
        IntQuat::new_i32(0, 0, 1, 0),
        IntQuat::new_i32(0, 0, 0, 1),
    ];
    id.denom = b_one();
    assert_eq!(quat_order_is_maximal(&id), 0);
}
