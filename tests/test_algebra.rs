use isogeny::quaternion::algebra::{IntQuat, RatQuat, QuatConfig};
use isogeny::bigint::{BigIntAlg, backends::CryptoInt};

use num_bigint::BigInt;
use std::sync::LazyLock;

// ========================================================================
// Test Configurations
// ========================================================================

#[derive(PartialEq)]
#[derive(Clone)]
pub struct P7;
static P7_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(7));
impl QuatConfig<BigInt> for P7 {
    fn p() -> &'static BigInt {
        &P7_VAL
    }
}

type IntQ7 = IntQuat<BigInt, P7>;
type RatQ7 = RatQuat<BigInt, P7>;

const LIMBS: usize = 4;
pub struct P11;
static P11_VAL: LazyLock<CryptoInt<LIMBS>> = LazyLock::new(|| CryptoInt::from_i32(11));

impl QuatConfig<CryptoInt<LIMBS>> for P11 {
    fn p() -> &'static CryptoInt<LIMBS> {
        &P11_VAL
    }
}

type IntQ11 = IntQuat<CryptoInt<LIMBS>, P11>;
type RatQ11 = RatQuat<CryptoInt<LIMBS>, P11>;

// ========================================================================
// Translated SQISign Tests (from algebra.c)
// ========================================================================

#[test]
fn quat_test_init_set_ui() {
    // Tests retrieving the algebra parameter
    assert_eq!(P7::p().clone(), BigInt::from(7));
    assert_eq!(P11::p().clone(), CryptoInt::from_i32(11));
}

#[test]
fn quat_test_alg_coord_mul() {
    // Test for p = 7
    let a7 = IntQ7::new_i32(152, 57, 190, 28);
    let b7 = IntQ7::new_i32(165, 35, 231, 770);
    let c7 = &a7 * &b7;

    assert_eq!(c7.coords[0], BigInt::from(-435065));
    assert_eq!(c7.coords[1], BigInt::from(993549));
    assert_eq!(c7.coords[2], BigInt::from(23552));
    assert_eq!(c7.coords[3], BigInt::from(128177));

    // Test for p = 11 (Same coords, different algebra)
    let a11 = IntQ11::new_i32(152, 57, 190, 28);
    let b11 = IntQ11::new_i32(165, 35, 231, 770);
    let c11 = &a11 * &b11;

    assert_eq!(c11.coords[0], CryptoInt::from_i32(-696865));
    assert_eq!(c11.coords[1], CryptoInt::from_i32(1552877));
    assert_eq!(c11.coords[2], CryptoInt::from_i32(23552));
    assert_eq!(c11.coords[3], CryptoInt::from_i32(128177));

    // Test for ones, p = 7
    let ones = IntQ7::new_i32(1, 1, 1, 1);
    let c_ones = &ones * &ones;

    assert_eq!(c_ones.coords[0], BigInt::from(-14));
    assert_eq!(c_ones.coords[1], BigInt::from(2));
    assert_eq!(c_ones.coords[2], BigInt::from(2));
    assert_eq!(c_ones.coords[3], BigInt::from(2));
}

#[test]
fn quat_test_alg_add() {
    let a = RatQ7::new_i32(9, -12, 0, -7, 19);
    let b = RatQ7::new_i32(3, -6, 2, 7, -19);

    // Using add_lazy to match C's unreduced outputs
    let c = a.add_lazy(&b);
    assert_eq!(c.denom, BigInt::from(9));
    assert_eq!(c.num.coords[0], BigInt::from(-30));
    assert_eq!(c.num.coords[1], BigInt::from(6));
    assert_eq!(c.num.coords[2], BigInt::from(14));
    assert_eq!(c.num.coords[3], BigInt::from(-38));

    let b2 = RatQ7::new_i32(6, -6, 2, 7, -19);
    let c2 = a.add_lazy(&b2);
    assert_eq!(c2.denom, BigInt::from(18));
    assert_eq!(c2.num.coords[0], BigInt::from(-42));
    assert_eq!(c2.num.coords[1], BigInt::from(6));
    assert_eq!(c2.num.coords[2], BigInt::from(7));
    assert_eq!(c2.num.coords[3], BigInt::from(-19));

    let c3 = a.add_lazy(&a).add_lazy(&a);
    assert_eq!(c3.denom, BigInt::from(9));
    assert_eq!(c3.num.coords[0], BigInt::from(-36));
    assert_eq!(c3.num.coords[2], BigInt::from(-21));
}

#[test]
fn quat_test_alg_sub() {
    let a = RatQ7::new_i32(9, -12, 0, -7, 19);
    let b = RatQ7::new_i32(3, -6, 2, 7, -19);

    let c = a.sub_lazy(&b);
    assert_eq!(c.denom, BigInt::from(9));
    assert_eq!(c.num.coords[0], BigInt::from(6));   // -12 - 3*(-6)
    assert_eq!(c.num.coords[1], BigInt::from(-6));  // -3*2
    assert_eq!(c.num.coords[2], BigInt::from(-28)); // -7 - 3*7
    assert_eq!(c.num.coords[3], BigInt::from(76));  // 19 - 3*(-19)

    let b2 = RatQ7::new_i32(6, -6, 2, 7, -19);
    let c2 = a.sub_lazy(&a).sub_lazy(&b2);
    assert_eq!(c2.denom, BigInt::from(18));
    assert_eq!(c2.num.coords[0], BigInt::from(18));
}

#[test]
fn quat_test_alg_mul() {
    let a7 = RatQ7::new_i32(76, 152, 57, 190, 28);
    let b7 = RatQ7::new_i32(385, 165, 35, 231, 770);

    let c7 = a7.mul_lazy(&b7);
    assert_eq!(c7.denom, BigInt::from(29260));
    assert_eq!(c7.num.coords[0], BigInt::from(-435065));
    assert_eq!(c7.num.coords[1], BigInt::from(993549));
    assert_eq!(c7.num.coords[2], BigInt::from(23552));
    assert_eq!(c7.num.coords[3], BigInt::from(128177));

    let a11 = RatQ11::new_i32(76, 152, 57, 190, 28);
    let b11 = RatQ11::new_i32(385, 165, 35, 231, 770);

    let c11 = a11.mul_lazy(&b11);
    assert_eq!(c11.denom, CryptoInt::from_i32(29260));
    assert_eq!(c11.num.coords[0], CryptoInt::from_i32(-696865));
    assert_eq!(c11.num.coords[1], CryptoInt::from_i32(1552877));
}

#[test]
fn quat_test_alg_norm() {
    let a11 = RatQ11::new_i32(2, 1, 5, 7, 2);
    let (num, denom) = a11.norm();
    assert_eq!(num, CryptoInt::from_i32(609));
    assert_eq!(denom, CryptoInt::from_i32(4));

    // Unreduced input
    let a11_unreduced = RatQ11::new_i32(4, 2, 10, 14, 4);
    let (num_u, denom_u) = a11_unreduced.norm();
    assert_eq!(num_u, CryptoInt::from_i32(609));
    assert_eq!(denom_u, CryptoInt::from_i32(4));

    let b11 = RatQ11::new_i32(76, 152, 57, 190, 28);
    let (num_b, denom_b) = b11.norm();
    assert_eq!(num_b, CryptoInt::from_i32(432077));
    assert_eq!(denom_b, CryptoInt::from_i32(5776));
}

#[test]
fn quat_test_alg_scalar() {
    let elem1 = RatQ7::new_i32(1, 1, 0, 0, 0);
    assert_eq!(elem1.denom, BigInt::from(1));
    assert_eq!(elem1.num.coords[0], BigInt::from(1));

    let elem2 = RatQ7::new_i32(9, 5, 0, 0, 0);
    assert_eq!(elem2.denom, BigInt::from(9));
    assert_eq!(elem2.num.coords[0], BigInt::from(5));
}

#[test]
fn quat_test_alg_conj() {
    let a = RatQ7::new_i32(25, 0, 0, 0, 7);
    let conj_a = a.conj();
    assert_eq!(conj_a.denom, BigInt::from(25));
    assert_eq!(conj_a.num.coords[3], BigInt::from(-7));

    let b = RatQ7::new_i32(25, -125, 2, 0, -30);
    let conj_b = b.conj();
    assert_eq!(conj_b.num.coords[0], BigInt::from(-125));
    assert_eq!(conj_b.num.coords[1], BigInt::from(-2));
    assert_eq!(conj_b.num.coords[3], BigInt::from(30));
}

#[test]
fn quat_test_alg_normalize() {
    let mut x = RatQ7::new_i32(-25, -125, 2, 0, -30);
    x.normalize();
    assert_eq!(x.denom, BigInt::from(25));
    assert_eq!(x.num.coords[0], BigInt::from(125));
    assert_eq!(x.num.coords[1], BigInt::from(-2));
    assert_eq!(x.num.coords[3], BigInt::from(30));

    let mut y = RatQ7::new_i32(48, -36, 18, 0, -300);
    y.normalize();
    assert_eq!(y.denom, BigInt::from(8));
    assert_eq!(y.num.coords[0], BigInt::from(-6));
    assert_eq!(y.num.coords[1], BigInt::from(3));
    assert_eq!(y.num.coords[3], BigInt::from(-50));
}

#[test]
fn quat_test_alg_elem_equal() {
    fn check_semantic_eq(a: &RatQ7, b: &RatQ7) -> bool {
        let mut a_norm = a.clone();
        let mut b_norm = b.clone();
        a_norm.normalize();
        b_norm.normalize();
        a_norm.ct_eq(&b_norm) == u32::MAX
    }

    let a = RatQ7::new_i32(5, 1, -3, -2, 2);
    let b = RatQ7::new_i32(15, 3, -9, -6, 6);
    assert!(check_semantic_eq(&a, &b));

    let c = RatQ7::new_i32(15, 3, -9, -6, 3);
    assert!(!check_semantic_eq(&a, &c));
}

#[test]
fn quat_test_alg_elem_is_zero() {
    let a = RatQ7::new_i32(1, 0, 0, 0, 0);
    assert!(a.is_zero());

    let b = RatQ7::new_i32(56865, 0, 0, 0, 0);
    assert!(b.is_zero());

    let c = RatQ7::new_i32(56865, 0, 0, 0, 1);
    assert!(!c.is_zero());

    let d = RatQ7::new_i32(56865, 19, -500, 20, -2);
    assert!(!d.is_zero());
}

#[test]
fn quat_test_alg_elem_mul_by_scalar() {
    let elem = RatQ7::new_i32(2, 2, -4, 5, 25);
    let scalar = BigInt::from(6);

    let prod_num = &elem.num * &scalar;
    let prod = RatQ7::new(prod_num, elem.denom.clone());

    assert_eq!(prod.denom, BigInt::from(2));
    assert_eq!(prod.num.coords[0], BigInt::from(12));
    assert_eq!(prod.num.coords[1], BigInt::from(-24));
    assert_eq!(prod.num.coords[2], BigInt::from(30));
    assert_eq!(prod.num.coords[3], BigInt::from(150));
}
