use isogeny::quaternion::algebra::{QuatConfig, IntQuat, RatQuat};
use isogeny::bigint::BigIntAlg;
use isogeny::quaternion::lattice::QuatLattice;
use num_bigint::BigInt;
use std::sync::LazyLock;

// ========================================================================
// Test Configurations
// ========================================================================

pub struct P19;
static P19_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(19));
impl QuatConfig<BigInt> for P19 { fn p() -> &'static BigInt { &P19_VAL } }

pub struct P23;
static P23_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(23));
impl QuatConfig<BigInt> for P23 { fn p() -> &'static BigInt { &P23_VAL } }

#[derive(Debug)]
pub struct P103;
static P103_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(103));
impl QuatConfig<BigInt> for P103 { fn p() -> &'static BigInt { &P103_VAL } }

type TestLattice = QuatLattice<BigInt, P103>;
type TestQuat = IntQuat<BigInt, P103>;
type TestRatQuat = RatQuat<BigInt, P103>;

type TestLattice19 = QuatLattice<BigInt, P19>;
type TestQuat19 = IntQuat<BigInt, P19>;

type TestLattice23 = QuatLattice<BigInt, P23>;
type TestQuat23 = IntQuat<BigInt, P23>;
type TestRatQuat23 = RatQuat<BigInt, P23>;

#[inline]
fn b_zero() -> BigInt { <BigInt as BigIntAlg>::zero() }
#[inline]
fn b(val: i32) -> BigInt { <BigInt as BigIntAlg>::from_i32(val) }

fn identity_gens<P: QuatConfig<BigInt>>() -> [IntQuat<BigInt, P>; 4] {
    [
        IntQuat::new_i32(1, 0, 0, 0),
        IntQuat::new_i32(0, 1, 0, 0),
        IntQuat::new_i32(0, 0, 1, 0),
        IntQuat::new_i32(0, 0, 0, 1),
    ]
}

// ========================================================================
// Tests
// ========================================================================

#[test]
fn test_lattice_equal() {
    let mut lat = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat.generators = identity_gens();
    cmp.generators = identity_gens();
    assert_eq!(TestLattice::equal(&lat, &cmp), u32::MAX);

    lat.denom = b(5);
    cmp.denom = b(4);
    assert_eq!(TestLattice::equal(&lat, &cmp), 0);

    lat.denom = b(1);
    cmp.denom = b(-1);
    assert_eq!(TestLattice::equal(&lat, &cmp), u32::MAX);

    lat.denom = b(3);
    cmp.denom = b(3);
    assert_eq!(TestLattice::equal(&lat, &cmp), u32::MAX);

    // Transposed from original matrix
    lat.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    lat.denom = b(6);
    lat.hnf();

    cmp.generators = lat.generators.clone();
    cmp.denom = b(6);
    assert_eq!(TestLattice::equal(&lat, &cmp), u32::MAX);

    cmp.denom = b(-7);
    assert_eq!(TestLattice::equal(&lat, &cmp), 0);

    cmp.denom = b(6);
    cmp.generators[3].coords[3] = b(165);
    assert_eq!(TestLattice::equal(&lat, &cmp), 0);
}

#[test]
fn test_lattice_inclusion() {
    let mut lat = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat.generators = identity_gens();
    cmp.generators = identity_gens();
    assert_eq!(TestLattice::inclusion(&lat, &cmp), u32::MAX);

    lat.denom = b(5);
    cmp.denom = b(4);
    assert_eq!(TestLattice::inclusion(&lat, &cmp), 0);

    lat.denom = b(1);
    cmp.denom = b(3);
    assert_eq!(TestLattice::inclusion(&lat, &cmp), u32::MAX);

    lat.denom = b(3);
    cmp.denom = b(3);
    assert_eq!(TestLattice::inclusion(&lat, &cmp), u32::MAX);

    lat.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    lat.denom = b(6);
    lat.hnf();

    cmp.generators = lat.generators.clone();
    cmp.denom = b(6);
    assert_eq!(TestLattice::inclusion(&lat, &cmp), u32::MAX);

    cmp.denom = b(12);
    assert_eq!(TestLattice::inclusion(&lat, &cmp), u32::MAX);

    cmp.denom = b(6);
    cmp.generators[3].coords[3] = b(165);
    assert_eq!(TestLattice::inclusion(&lat, &cmp), 0);
}

#[test]
fn test_lattice_reduce_denom() {
    let mut lat = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    let s = 15i32;
    for idx in 0..4 {
        let i = idx as i32;
        lat.generators[idx] = TestQuat::new_i32(i * s, (i + 1) * s, (i + 2) * s, (i + 3) * s);
        cmp.generators[idx] = TestQuat::new_i32(i, i + 1, i + 2, i + 3);
    }
    lat.denom = b(4 * s);
    cmp.denom = b(4);

    let mut red = lat.clone();
    red.reduce_denom();

    assert_eq!(red.ct_eq(&cmp), u32::MAX);

    lat.reduce_denom();
    assert_eq!(lat.ct_eq(&cmp), u32::MAX);
}

#[test]
fn test_lattice_conjugate_without_hnf() {
    let mut lat = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat.generators = [
        TestQuat::new_i32(4, 0, 0, 0),
        TestQuat::new_i32(0, -2, -1, 0),
        TestQuat::new_i32(0, 0, -1, 0),
        TestQuat::new_i32(1, 0, 0, -3),
    ];
    lat.denom = b(6);

    cmp.generators = [
        TestQuat::new_i32(4, 0, 0, 0),
        TestQuat::new_i32(0, 2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(1, 0, 0, 3),
    ];
    cmp.denom = b(6);

    lat.hnf();

    let mut conj = lat.conjugate_without_hnf();
    conj.hnf();
    cmp.hnf();

    assert_eq!(TestLattice::equal(&conj, &cmp), u32::MAX);

    let mut conj_conj = conj.conjugate_without_hnf();
    conj_conj.hnf();
    assert_eq!(TestLattice::equal(&conj_conj, &lat), u32::MAX);
}

#[test]
fn test_lattice_dual_without_hnf() {
    let mut lat = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    lat.denom = b(6);

    cmp.generators = [
        TestQuat::new_i32(6, 0, 0, 0),
        TestQuat::new_i32(0, 3, 0, 0),
        TestQuat::new_i32(0, 0, 6, 0),
        TestQuat::new_i32(0, 0, 0, 2),
    ];
    cmp.denom = b(1);

    lat.hnf();

    let mut dual = lat.dual_without_hnf();
    dual.hnf();
    cmp.hnf();

    assert_eq!(TestLattice::equal(&dual, &cmp), u32::MAX);
    assert_eq!(!TestLattice::equal(&dual, &lat), u32::MAX);

    let mut dual_dual = dual.dual_without_hnf();
    dual_dual.hnf();
    assert_eq!(TestLattice::equal(&dual_dual, &lat), u32::MAX);
}

#[test]
fn test_lattice_add() {
    let mut lat1 = TestLattice::zero();
    let mut lat2 = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat1.generators = [
        TestQuat::new_i32(44, 0, 0, 0),
        TestQuat::new_i32(0, 5, 0, 0),
        TestQuat::new_i32(3, 0, 3, 0),
        TestQuat::new_i32(32, 0, 0, 1),
    ];
    lat2.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, 2, 0, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(0, 0, 0, 3),
    ];
    cmp.generators = [
        TestQuat::new_i32(2, 0, 0, 0),
        TestQuat::new_i32(0, 1, 0, 0),
        TestQuat::new_i32(1, 0, 1, 0),
        TestQuat::new_i32(0, 0, 0, 3),
    ];
    lat1.denom = b(4);
    lat2.denom = b(6);
    cmp.denom = b(12);

    let sum = &lat1 + &lat2;
    assert_eq!(sum.ct_eq(&cmp), u32::MAX);

    lat1.generators = [
        TestQuat::new_i32(4, 0, 0, 0),
        TestQuat::new_i32(0, 5, 0, 0),
        TestQuat::new_i32(3, 0, 3, 0),
        TestQuat::new_i32(0, 0, 0, 7),
    ];
    lat2.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    lat1.denom = b(4);
    lat2.denom = b(6);

    let sum2 = &lat1 + &lat2;
    assert_eq!(sum2.ct_eq(&cmp), u32::MAX);

    cmp.generators = lat2.generators.clone();
    cmp.denom = lat2.denom.clone();
    cmp.hnf();

    let sum_self = &lat2 + &lat2;
    assert_eq!(sum_self.ct_eq(&cmp), u32::MAX);
}

#[test]
fn test_lattice_intersect() {
    let mut lat1 = TestLattice::zero();
    let mut lat2 = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat1.generators = [
        TestQuat::new_i32(4, 0, 0, 0),
        TestQuat::new_i32(0, 5, 0, 0),
        TestQuat::new_i32(3, 0, 3, 0),
        TestQuat::new_i32(0, 0, 0, 7),
    ];
    lat2.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    lat1.denom = b(4);
    lat2.denom = b(6);
    lat1.hnf();
    lat2.hnf();

    cmp.generators = [
        TestQuat::new_i32(2, 0, 0, 0),
        TestQuat::new_i32(0, 10, 0, 0),
        TestQuat::new_i32(1, 0, 3, 0),
        TestQuat::new_i32(0, 0, 0, 7),
    ];
    cmp.denom = b(2);

    let inter = TestLattice::intersect(&lat1, &lat2);
    assert_eq!(TestLattice::equal(&inter, &cmp), u32::MAX);

    let inter2 = TestLattice::intersect(&lat2, &lat1);
    assert_eq!(TestLattice::equal(&inter2, &cmp), u32::MAX);

    cmp.generators = lat1.generators.clone();
    cmp.denom = lat1.denom.clone();
    let inter_self = TestLattice::intersect(&lat1, &lat1);
    assert_eq!(TestLattice::equal(&inter_self, &cmp), u32::MAX);
}

#[test]
fn test_lattice_alg_elem_mul() {
    let mut lat = TestLattice23::zero();
    let mut cmp = TestLattice23::zero();

    let elem = TestRatQuat23::new_i32(2, 3, 4, -1, 0);

    lat.generators = [
        TestQuat23::new_i32(11, 0, 0, 0),
        TestQuat23::new_i32(2, -13, 0, 0),
        TestQuat23::new_i32(0, 0, 15, 0),
        TestQuat23::new_i32(0, -1, 0, -4),
    ];
    lat.denom = b(5);
    lat.hnf();

    let prod = &lat * &elem;

    cmp.generators = [
        TestQuat23::new_i32(33, 44, -11, 0),
        TestQuat23::new_i32(27 - 4 * 13, 36 + 3 * 13, -9, -13),
        TestQuat23::new_i32(15 * 23, 0, 45, -60),
        TestQuat23::new_i32(-4, 3 + 23 * 4, 16, -1 + 4 * 3),
    ];
    cmp.denom = b(10);
    cmp.hnf();

    assert_eq!(TestLattice23::equal(&cmp, &prod), u32::MAX);

    let prod2 = &lat * &elem;
    assert_eq!(TestLattice23::equal(&cmp, &prod2), u32::MAX);
}

#[test]
fn test_lattice_mul() {
    let mut lat1 = TestLattice19::zero();
    let mut lat2 = TestLattice19::zero();
    let mut cmp = TestLattice19::zero();

    lat1.generators = [
        TestQuat19::new_i32(44, 0, 0, 0),
        TestQuat19::new_i32(0, 5, 0, 0),
        TestQuat19::new_i32(3, 0, 3, 0),
        TestQuat19::new_i32(32, 0, 0, 1),
    ];
    lat2.generators = [
        TestQuat19::new_i32(1, 0, 0, 0),
        TestQuat19::new_i32(0, 2, 0, 0),
        TestQuat19::new_i32(0, 0, 1, 0),
        TestQuat19::new_i32(0, 0, 0, 3),
    ];
    cmp.generators = [
        TestQuat19::new_i32(1, 0, 0, 0),
        TestQuat19::new_i32(0, 1, 0, 0),
        TestQuat19::new_i32(0, 0, 1, 0),
        TestQuat19::new_i32(0, 0, 0, 1),
    ];
    lat1.denom = b(4);
    lat2.denom = b(6);
    cmp.denom = b(24);

    let prod = &lat1 * &lat2;
    assert_eq!(prod.ct_eq(&cmp), u32::MAX);

    lat1.generators = [
        TestQuat19::new_i32(4, 0, 0, 0),
        TestQuat19::new_i32(0, 5, 0, 0),
        TestQuat19::new_i32(3, 0, 3, 0),
        TestQuat19::new_i32(0, 0, 0, 7),
    ];
    lat2.generators = [
        TestQuat19::new_i32(1, 0, 0, 0),
        TestQuat19::new_i32(0, -2, 1, 0),
        TestQuat19::new_i32(0, 0, 1, 0),
        TestQuat19::new_i32(-1, 0, 0, -3),
    ];
    lat1.denom = b(4);
    lat2.denom = b(6);

    let prod2 = &lat1 * &lat2;
    assert_eq!(prod2.ct_eq(&cmp), u32::MAX);

    cmp.generators = [
        TestQuat19::new_i32(1, 0, 0, 0),
        TestQuat19::new_i32(0, 1, 0, 0),
        TestQuat19::new_i32(0, 0, 1, 0),
        TestQuat19::new_i32(0, 0, 0, 1),
    ];
    cmp.denom = b(36);

    let prod_self = &lat2 * &lat2;
    assert_eq!(prod_self.ct_eq(&cmp), u32::MAX);
}

#[test]
fn test_lattice_contains() {
    let mut lat = TestLattice::zero();

    lat.generators = [
        TestQuat::new_i32(4, 0, 0, 0),
        TestQuat::new_i32(0, 5, 0, 0),
        TestQuat::new_i32(3, 0, 3, 0),
        TestQuat::new_i32(0, 0, 0, 7),
    ];
    lat.denom = b(4);

    let x = TestRatQuat::new_i32(3, 1, -2, 26, 9);
    assert!(lat.contains(&x).is_none());

    lat.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    lat.denom = b(6);
    lat.hnf();

    let coord_res = lat.contains(&x);
    assert!(coord_res.is_some());
    let coord = coord_res.unwrap();

    assert_eq!(coord[0], b(2));
    assert_eq!(coord[1], b(-2));
    assert_eq!(coord[2], b(52));
    assert_eq!(coord[3], b(6));
}

#[test]
fn test_lattice_index() {
    let mut sublat = TestLattice::zero();
    let mut overlat = TestLattice::zero();

    overlat.generators = identity_gens();
    overlat.denom = b(2);

    sublat.generators = [
        TestQuat::new_i32(2, 0, 0, 0),
        TestQuat::new_i32(0, 4, 0, 0),
        TestQuat::new_i32(1, 2, 1, 0),
        TestQuat::new_i32(0, 3, 0, 1),
    ];
    sublat.denom = b(2);

    let index = TestLattice::index(&sublat, &overlat);
    assert_eq!(index, b(8));
}

#[test]
fn test_lattice_hnf() {
    let mut lat = TestLattice::zero();
    let mut cmp = TestLattice::zero();

    lat.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, -2, 1, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(-1, 0, 0, -3),
    ];
    cmp.generators = [
        TestQuat::new_i32(1, 0, 0, 0),
        TestQuat::new_i32(0, 2, 0, 0),
        TestQuat::new_i32(0, 0, 1, 0),
        TestQuat::new_i32(0, 0, 0, 3),
    ];

    cmp.denom = b(6);
    lat.denom = b(6);

    lat.hnf();

    assert_eq!(lat.ct_eq(&cmp), u32::MAX);
}

#[test]
fn test_lattice_gram() {
    let mut lattice = TestLattice::zero();

    lattice.generators = [
        TestQuat::new_i32(202, 0, 0, 0),
        TestQuat::new_i32(0, 202, 0, 0),
        TestQuat::new_i32(158, 149, 1, 0),
        TestQuat::new_i32(53, 158, 0, 1),
    ];
    lattice.denom = b(2);

    let gram = lattice.gram();

    let elem1 = TestRatQuat::new_i32(2, 360, 149, 1, 0);
    let elem2 = TestRatQuat::new_i32(2, 53, 360, 0, 1);

    let vec1_opt = lattice.contains(&elem1);
    let vec2_opt = lattice.contains(&elem2);
    assert!(vec1_opt.is_some());
    assert!(vec2_opt.is_some());

    let vec1 = vec1_opt.unwrap();
    let vec2 = vec2_opt.unwrap();

    let prod = &elem1 * &elem2.conj();

    let norm1 = (prod.num.coords[0].clone() * b(2)) / prod.denom.clone();

    let mut vec1_gram = [b_zero(), b_zero(), b_zero(), b_zero()];
    for i in 0..4 {
        let mut sum = b_zero();
        for j in 0..4 {
            sum = sum + gram[i][j].clone() * vec1[j].clone();
        }
        vec1_gram[i] = sum;
    }

    let mut norm2 = b_zero();
    for i in 0..4 {
        norm2 = norm2 + vec1_gram[i].clone() * vec2[i].clone();
    }

    norm2 = norm2 / lattice.denom.clone();
    norm2 = norm2 / lattice.denom.clone();

    assert_eq!(norm1, norm2);
}
