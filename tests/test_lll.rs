use isogeny::bigint::rat::BigRat;
use isogeny::bigint::BigIntAlg;
use isogeny::quaternion::algebra::{IntQuat, QuatConfig, RatQuat};
use isogeny::quaternion::ideal::QuatLeftIdeal;
use isogeny::quaternion::lattice::{IntLattice, RatLattice};
use isogeny::quaternion::lll::{quat_lideal_reduce_basis, quat_lll_core};
use num_bigint::BigInt;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::sync::LazyLock;

// ========================================================================
// Test Configurations
// ========================================================================

pub struct P3;
static P3_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(3));
impl QuatConfig<BigInt> for P3 {
    fn p() -> &'static BigInt {
        &P3_VAL
    }
}

pub struct P19;
static P19_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(19));
impl QuatConfig<BigInt> for P19 {
    fn p() -> &'static BigInt {
        &P19_VAL
    }
}

pub struct P103;
static P103_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(103));
impl QuatConfig<BigInt> for P103 {
    fn p() -> &'static BigInt {
        &P103_VAL
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
// LLL Verification Helpers (Translated from lll_verification.c)
// ========================================================================

fn lll_bilinear<T: BigIntAlg, P: QuatConfig<T>>(
    vec0: &[BigRat<T>; 4],
    vec1: &[BigRat<T>; 4],
) -> BigRat<T> {
    let p = BigRat::from_integer(P::p().clone());
    let mut sum = vec0[0].clone() * vec1[0].clone() + vec0[1].clone() * vec1[1].clone();
    sum = sum + p.clone() * vec0[2].clone() * vec1[2].clone();
    sum = sum + p * vec0[3].clone() * vec1[3].clone();
    sum
}

fn lll_gram_schmidt<T: BigIntAlg, P: QuatConfig<T>>(
    mat: &IntLattice<T, P>,
) -> [[BigRat<T>; 4]; 4] {
    let mut work = [
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
    ];

    for i in 0..4 {
        for j in 0..4 {
            work[i][j] = BigRat::from_integer(mat.generators[i].coords[j].clone());
        }
    }

    for i in 0..4 {
        let norm = lll_bilinear::<T, P>(&work[i], &work[i]);
        let norm_inv = BigRat::new(norm.den, norm.num); // Invert norm

        for j in (i + 1)..4 {
            let b_val = lll_bilinear::<T, P>(&work[i], &work[j]);
            let coeff = b_val * norm_inv.clone();

            for k in 0..4 {
                let prod = coeff.clone() * work[i][k].clone();
                work[j][k] = work[j][k].clone() - prod;
            }
        }
    }

    work
}

fn lll_verify<T: BigIntAlg, P: QuatConfig<T>>(mat: &IntLattice<T, P>) -> bool {
    let delta = BigRat::new(T::from_i32(99), T::from_i32(100)); // 0.99
    let eta = BigRat::new(T::from_i32(501), T::from_i32(1000)); // 0.501

    let b_star = lll_gram_schmidt(mat);
    let mut res = true;

    // Check size reduction
    for i in 0..4 {
        for j in 0..i {
            let mut b_i = [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()];
            for k in 0..4 {
                b_i[k] = BigRat::from_integer(mat.generators[i].coords[k].clone());
            }

            let b_val = lll_bilinear::<T, P>(&b_star[j], &b_i);
            let norm = lll_bilinear::<T, P>(&b_star[j], &b_star[j]);

            let mut mu = b_val * BigRat::new(norm.den, norm.num); // Division
            if mu < BigRat::zero() {
                mu = BigRat::zero() - mu;
            }

            if mu > eta {
                res = false;
            }
        }
    }

    // Check Lovász condition
    for i in 1..4 {
        let mut b_i = [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()];
        for k in 0..4 {
            b_i[k] = BigRat::from_integer(mat.generators[i].coords[k].clone());
        }

        let b_val = lll_bilinear::<T, P>(&b_star[i - 1], &b_i);
        let norm_prev = lll_bilinear::<T, P>(&b_star[i - 1], &b_star[i - 1]);

        let mu = b_val * BigRat::new(norm_prev.den.clone(), norm_prev.num.clone());
        let mu_sq = mu.clone() * mu;

        let factor = delta.clone() - mu_sq;
        let lhs = lll_bilinear::<T, P>(&b_star[i], &b_star[i]);
        let rhs = factor * norm_prev;

        if lhs < rhs {
            res = false;
        }
    }

    res
}

// ========================================================================
// Tests
// ========================================================================

#[test]
fn test_lll_bigrat_consts() {
    let mut t = BigRat::<BigInt>::new(b(123), b(-123));

    assert!(!t.num.is_zero());
    assert_eq!(t.num, b(-1));
    assert_eq!(t.den, b(1));

    t = BigRat::new(b(123), b(123));
    assert_eq!(t.num, b(1));
    assert_eq!(t.den, b(1));

    t = BigRat::new(b(0), b(123));
    assert!(t.num.is_zero());
    assert_eq!(t.den, b(1));
}

#[test]
fn test_lll_verify_fail_conditions() {
    // Tests that unreduced matrices correctly fail the LLL verification[cite: 10]
    let mut mat = IntLattice::<BigInt, P3>::zero();
    mat.generators[0] = IntQuat::new_i32(0, 2, 3, -14);
    mat.generators[1] = IntQuat::new_i32(2, -1, -4, -8);
    mat.generators[2] = IntQuat::new_i32(1, -2, 1, 0);
    mat.generators[3] = IntQuat::new_i32(1, 1, 0, 7);

    assert!(!lll_verify(&mat), "Unreduced matrix 1 should fail verification");

    let mut mat2 = IntLattice::<BigInt, P103>::zero();
    mat2.generators[0] = IntQuat::new_i32(3, 0, 90, -86);
    mat2.generators[1] = IntQuat::new_i32(11, 15, 12, 50);
    mat2.generators[2] = IntQuat::new_i32(1, -2, 0, 3);
    mat2.generators[3] = IntQuat::new_i32(-1, 0, 5, 5);

    assert!(!lll_verify(&mat2), "Unreduced matrix 2 should fail verification");
}

#[test]
fn test_lll_lattice_lll() {
    let mut lat = RatLattice::<BigInt, P103>::zero();

    lat.denom = b(60);
    lat.basis.generators[0] = IntQuat::new_i32(3, 1, 0, -19);
    lat.basis.generators[1] = IntQuat::new_i32(7, 0, 12, 0);
    lat.basis.generators[2] = IntQuat::new_i32(0, 0, 5, 0);
    lat.basis.generators[3] = IntQuat::new_i32(0, -6, 0, 3);
    lat.hnf();

    let mut red = lat.basis.clone();
    let mut gram = red.gram();
    quat_lll_core(&mut gram, &mut red);

    assert!(lll_verify(&red), "Reduced basis failed LLL verification parameters");

    let mut test_lat = RatLattice::<BigInt, P103>::zero();
    test_lat.denom = lat.denom.clone();
    test_lat.basis = red;
    test_lat.hnf();

    assert_eq!(RatLattice::equal(&test_lat, &lat), u32::MAX, "Reduced lattice is not equal to original");
}

#[test]
fn test_lll_randomized_lattice_lll() {
    let mut rng = StdRng::seed_from_u64(987654321);

    for _ in 0..20 {
    let mut lat = RatLattice::<BigInt, P103>::zero();

        // Random invertible matrix generator
        let mut det = b_zero();
        while det.is_zero() {
            for i in 0..4 {
                for j in 0..4 {
                    lat.basis.generators[i].coords[j] = b(rng.random_range(-20..20));
                }
            }
            det = lat.basis.inv_with_det(None);
        }

        lat.denom = b(rng.random_range(1..100));
        lat.hnf();

        let mut red = lat.basis.clone();
        let mut gram = red.gram();
        quat_lll_core(&mut gram, &mut red);

        assert!(lll_verify(&red), "Randomized reduced basis failed verification");

        let mut test_lat = RatLattice::<BigInt, P103>::zero();
        test_lat.denom = lat.denom.clone();
        test_lat.basis = red;
        test_lat.hnf();

        assert_eq!(RatLattice::equal(&test_lat, &lat), u32::MAX, "Reduced lattice altered lattice space");
    }
}

#[test]
fn test_lideal_reduce_basis() {
    let alg_order = o0_lattice::<P19>();

    // Test generator defined as 1, 1, 2, 8, 8 mapping to denom=1, c0=1, c1=2, c2=8, c3=8
    let init_helper_int = IntQuat::new_i32(1, 2, 8, 8);
    let init_helper = RatQuat::new(init_helper_int, b_one());

    let lideal = QuatLeftIdeal::create_principal(&alg_order, &init_helper);

    let mut red = IntLattice::zero();
    let mut gram = [
        [b_zero(), b_zero(), b_zero(), b_zero()],
        [b_zero(), b_zero(), b_zero(), b_zero()],
        [b_zero(), b_zero(), b_zero(), b_zero()],
        [b_zero(), b_zero(), b_zero(), b_zero()]
    ];

    quat_lideal_reduce_basis(&mut red, &mut gram, &lideal);

    assert!(lll_verify(&red), "Ideal reduced basis failed LLL verification");

    let mut test_lat = RatLattice::<BigInt, P19>::zero();
    test_lat.basis = red;
    test_lat.denom = lideal.lattice.denom.clone();
    test_lat.hnf();

    assert_eq!(RatLattice::equal(&lideal.lattice, &test_lat), u32::MAX, "Reduced ideal lattice not equal to original");
}
