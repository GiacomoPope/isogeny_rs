use isogeny::quaternion::algebra::{IntQuat, QuatConfig};
use isogeny::quaternion::hnf_ext::HnfModExt;
use isogeny::bigint::BigIntAlg;
use isogeny::quaternion::hnf::*;
use isogeny::quaternion::lattice::IntLattice;
use num_bigint::BigInt;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::sync::LazyLock;

// ========================================================================
// Test Configuration
// ========================================================================

pub struct P103;
static P103_VAL: LazyLock<BigInt> = LazyLock::new(|| BigInt::from(103));
impl QuatConfig<BigInt> for P103 {
    fn p() -> &'static BigInt {
        &P103_VAL
    }
}

type TestQuat = IntQuat<BigInt, P103>;

#[inline]
fn b_zero() -> BigInt { <BigInt as BigIntAlg>::zero() }
#[inline]
fn b_one() -> BigInt { <BigInt as BigIntAlg>::one() }

fn assert_is_hnf(mat: &[TestQuat; 4]) {
    let zero = b_zero();

    for i in 0..4 {
        for j in (i + 1)..4 {
            assert_eq!(
                mat[i].coords[j], zero,
                "HNF Failed: Not upper triangular at col {}, row {}. Value: {:?}", i, j, mat[i].coords[j]
            );
        }

        assert!(
            mat[i].coords[i] > zero,
            "HNF Failed: Diagonal pivot not positive at col {}, row {}. Value: {:?}", i, i, mat[i].coords[i]
        );

        for j in (i + 1)..4 {
            assert!(
                mat[j].coords[i] >= zero,
                "HNF Failed: Element to right is negative at col {}, row {}. Value: {:?}", j, i, mat[j].coords[i]
            );
            assert!(
                mat[j].coords[i] < mat[i].coords[i],
                "HNF Failed: Element to right >= row pivot at col {}, row {}. Value: {:?} >= {:?}",
                j, i, mat[j].coords[i], mat[i].coords[i]
            );
        }
    }
}

#[test]
fn test_positive_mod() {
    let modulo = BigInt::from_i32(3);

    assert_eq!(BigInt::from_i32(5).positive_mod(&modulo), BigInt::from_i32(2));
    assert_eq!(BigInt::from_i32(6).positive_mod(&modulo), b_zero());
    assert_eq!(BigInt::from_i32(-1).positive_mod(&modulo), BigInt::from_i32(2));
}

#[test]
fn test_centered_mod() {
    let modulo = BigInt::from_i32(3);

    assert_eq!(BigInt::from_i32(5).centered_mod(&modulo), BigInt::from_i32(-1));
    assert_eq!(BigInt::from_i32(4).centered_mod(&modulo), BigInt::from_i32(1));
    assert_eq!(BigInt::from_i32(6).centered_mod(&modulo), b_zero());
}

#[test]
fn test_xgcd_with_u_not_0() {
    let a = b_zero();
    let b = BigInt::from_i32(5);

    let (d, u, v) = xgcd_with_u_not_0(&a, &b);

    assert_eq!(d, BigInt::from_i32(5));
    assert!(!BigIntAlg::is_zero(&u), "xgcd_with_u_not_0 failed to guarantee u != 0");

    let res = a.clone() * u.clone() + b.clone() * v.clone();
    assert_eq!(res, d, "Bezout identity failed");
}

#[test]
fn test_mat_4x4_inv_with_det_as_denom() {
    let mut mat: IntLattice<BigInt, P103> = IntLattice::zero();
    mat.generators[0].coords[0] = BigInt::from_i32(1);
    mat.generators[1].coords[1] = BigInt::from_i32(4);
    mat.generators[2].coords[2] = BigInt::from_i32(1);
    mat.generators[3].coords[3] = BigInt::from_i32(1);
    mat.generators[1].coords[0] = BigInt::from_i32(2);
    mat.generators[0].coords[1] = BigInt::from_i32(3);

    let mut adj = IntLattice::zero();
    let det = mat.inv_with_det(Some(&mut adj));

    assert_eq!(det, BigInt::from_i32(-2));

    let mut prod = [TestQuat::zero(), TestQuat::zero(), TestQuat::zero(), TestQuat::zero()];
    for i in 0..4 {
        for j in 0..4 {
            let mut sum = b_zero();
            for k in 0..4 {
                sum = sum + adj.generators[j].coords[k].clone() * mat.generators[k].coords[i].clone();
            }
            prod[j].coords[i] = sum;
        }
    }

    let mut expected = [TestQuat::zero(), TestQuat::zero(), TestQuat::zero(), TestQuat::zero()];
    for i in 0..4 {
        expected[i].coords[i] = det.clone();
    }

    for i in 0..4 {
        assert_eq!(prod[i].ct_eq(&expected[i]), u32::MAX, "Adjugate multiplication property failed");
    }
}

#[test]
fn test_quat_hnf_mod_core_exact_vectors() {
    let mut generators = vec![TestQuat::zero(); 8];
    generators[0].coords[0] = BigInt::from(4);
    generators[2].coords[0] = BigInt::from(3);
    generators[4].coords[0] = BigInt::from(1);
    generators[7].coords[0] = BigInt::from(-1);
    generators[1].coords[1] = BigInt::from(5);
    generators[5].coords[1] = BigInt::from(-2);
    generators[2].coords[2] = BigInt::from(3);
    generators[6].coords[2] = BigInt::from(1);
    generators[5].coords[2] = BigInt::from(1);
    generators[3].coords[3] = BigInt::from(7);
    generators[7].coords[3] = BigInt::from(-3);

    let det = BigInt::from(1);
    let hnf = quat_hnf_mod_core(&mut generators, &det);

    let mut identity = [TestQuat::zero(), TestQuat::zero(), TestQuat::zero(), TestQuat::zero()];
    for i in 0..4 {
        identity[i].coords[i] = b_one();
    }
    for i in 0..4 {
        assert_eq!(identity[i].ct_eq(&hnf[i]), u32::MAX, "Identity vector test failed");
    }

    let mut generators = vec![TestQuat::zero(); 8];
    generators[4].coords[0] = BigInt::from(438);
    generators[4].coords[1] = BigInt::from(400);
    generators[4].coords[2] = BigInt::from(156);
    generators[4].coords[3] = BigInt::from(-2);
    generators[5].coords[0] = BigInt::from(-400);
    generators[5].coords[1] = BigInt::from(438);
    generators[5].coords[2] = BigInt::from(2);
    generators[5].coords[3] = BigInt::from(156);
    generators[6].coords[0] = BigInt::from(-28826);
    generators[6].coords[1] = BigInt::from(-148);
    generators[6].coords[2] = BigInt::from(220);
    generators[6].coords[3] = BigInt::from(-122);
    generators[7].coords[0] = BigInt::from(586);
    generators[7].coords[1] = BigInt::from(-28426);
    generators[7].coords[2] = BigInt::from(278);
    generators[7].coords[3] = BigInt::from(218);

    let mut expected_cmp = [TestQuat::zero(), TestQuat::zero(), TestQuat::zero(), TestQuat::zero()];
    expected_cmp[0].coords[0] = BigInt::from(2321156);
    expected_cmp[1].coords[1] = BigInt::from(2321156);
    expected_cmp[2].coords[0] = BigInt::from(620252);
    expected_cmp[2].coords[1] = BigInt::from(365058);
    expected_cmp[2].coords[2] = BigInt::from(2);
    expected_cmp[3].coords[0] = BigInt::from(1956098);
    expected_cmp[3].coords[1] = BigInt::from(620252);
    expected_cmp[3].coords[3] = BigInt::from(2);

    let det = BigInt::parse_bytes(b"21551060705344", 10).unwrap();
    let hnf = quat_hnf_mod_core(&mut generators, &det);

    for i in 0..4 {
        assert_eq!(expected_cmp[i].ct_eq(&hnf[i]), u32::MAX, "GitHub #38 exact vector test failed");
    }
}

#[test]
fn test_quat_hnf_mod_core_fuzzing() {
    let mut rng = StdRng::seed_from_u64(123456789);

    for iteration in 0..250 {
        let n_gens = rng.random_range(4..12);

        let mut expected_hnf = [TestQuat::zero(), TestQuat::zero(), TestQuat::zero(), TestQuat::zero()];

        for i in 0..4 {
            expected_hnf[i].coords[i] = BigInt::from_i32(rng.random_range(1..25));
        }

        for i in 0..4 {
            for j in (i + 1)..4 {
                let d_val = i32::try_from(&expected_hnf[i].coords[i]).unwrap_or(1);
                expected_hnf[j].coords[i] = BigInt::from_i32(rng.random_range(0..d_val));
            }
        }

        assert_is_hnf(&expected_hnf);

        let mut det = b_one();
        for i in 0..4 {
            det = det * expected_hnf[i].coords[i].clone();
        }

        let mut gens = vec![TestQuat::zero(); n_gens];
        for i in 0..4 {
            gens[i] = expected_hnf[i].clone();
        }

        let n_ops = 300;
        for _ in 0..n_ops {
            let op = rng.random_range(0..3);
            if op == 0 {
                let r1 = rng.random_range(0..n_gens);
                let r2 = rng.random_range(0..n_gens);
                if r1 != r2 {
                    let scalar = BigInt::from_i32(rng.random_range(-5..5));
                    for c in 0..4 {
                        gens[r2].coords[c] = gens[r2].coords[c].clone() + scalar.clone() * gens[r1].coords[c].clone();
                    }
                }
            } else if op == 1 {
                let r1 = rng.random_range(0..n_gens);
                let r2 = rng.random_range(0..n_gens);
                if r1 != r2 {
                    gens.swap(r1, r2);
                }
            } else {
                let r1 = rng.random_range(0..n_gens);
                for c in 0..4 {
                    gens[r1].coords[c] = -gens[r1].coords[c].clone();
                }
            }
        }

        for i in 0..n_gens {
            for c in 0..4 {
                let scalar = BigInt::from_i32(rng.random_range(-3..4));
                gens[i].coords[c] = gens[i].coords[c].clone() + scalar * det.clone();
            }
        }

        let computed_hnf = quat_hnf_mod_core(&mut gens, &det);

        let mut eq = u32::MAX;
        for i in 0..4 { eq &= expected_hnf[i].ct_eq(&computed_hnf[i]); }
        if eq == 0 {
            println!("\nFUZZING FAILED it {}", iteration);
            println!("Determinant (Modulo): {}", det);
            println!("Expected HNF Matrix:");
            for i in 0..4 { println!("  {:?}", expected_hnf[i]); }
            println!("Computed HNF Matrix:");
            for i in 0..4 { println!("  {:?}", computed_hnf[i]); }
            panic!("quat_hnf_mod_core failed to reconstruct the original lattice.");
        }
    }
}
