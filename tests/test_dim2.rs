use isogeny::bigint::BigIntAlg;
use isogeny::quaternion::dim2::{
    dim2_lattice_bilinear, dim2_lattice_norm, dim2_lattice_short_basis, rounded_div, Mat2x2, Vec2,
};
use num_bigint::BigInt;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

// ========================================================================
// Helpers
// ========================================================================

#[inline]
fn b_zero() -> BigInt {
    <BigInt as BigIntAlg>::zero()
}

#[inline]
fn b(val: i32) -> BigInt {
    <BigInt as BigIntAlg>::from_i32(val)
}

// ========================================================================
// Tests
// ========================================================================

#[test]
fn test_rounded_div() {
    // Exact divisions
    assert_eq!(rounded_div(&b(10), &b(2)), b(5));
    assert_eq!(rounded_div(&b(-10), &b(2)), b(-5));
    assert_eq!(rounded_div(&b(10), &b(-2)), b(-5));
    assert_eq!(rounded_div(&b(-10), &b(-2)), b(5));

    // Rounding down (positive)
    assert_eq!(rounded_div(&b(10), &b(3)), b(3));
    assert_eq!(rounded_div(&b(13), &b(4)), b(3));

    // Rounding up (positive)
    assert_eq!(rounded_div(&b(11), &b(3)), b(4));
    assert_eq!(rounded_div(&b(15), &b(4)), b(4));

    // Rounding towards zero (negative)
    assert_eq!(rounded_div(&b(-10), &b(3)), b(-3));
    assert_eq!(rounded_div(&b(-13), &b(4)), b(-3));

    // Rounding away from zero (negative)
    assert_eq!(rounded_div(&b(-11), &b(3)), b(-4));
    assert_eq!(rounded_div(&b(-15), &b(4)), b(-4));

    // Tie-breaking behavior (should truncates towards zero)
    assert_eq!(rounded_div(&b(5), &b(2)), b(2));
    assert_eq!(rounded_div(&b(-5), &b(2)), b(-2));
}

#[test]
fn test_mat2x2_eval() {
    let mut mat = Mat2x2::zero();
    mat.m[0][0] = b(1);
    mat.m[0][1] = b(2);
    mat.m[1][0] = b(3);
    mat.m[1][1] = b(4);

    let vec = Vec2::new(b(5), b(6));
    let res = mat.eval(&vec);

    assert_eq!(res.coords[0], b(1 * 5 + 2 * 6));
    assert_eq!(res.coords[1], b(3 * 5 + 4 * 6));
}

#[test]
fn test_mat2x2_inv_with_det() {
    let mut mat = Mat2x2::zero();
    mat.m[0][0] = b(3);
    mat.m[0][1] = b(8);
    mat.m[1][0] = b(4);
    mat.m[1][1] = b(6);

    let (inv, det) = mat.inv_with_det_as_denom();

    assert_eq!(det, b(-14));

    // Inverse check: M * adj(M) = det(M) * I
    let mut prod = Mat2x2::zero();
    for i in 0..2 {
        for j in 0..2 {
            let mut sum = b_zero();
            for k in 0..2 {
                sum = sum + mat.m[i][k].clone() * inv.m[k][j].clone();
            }
            prod.m[i][j] = sum;
        }
    }

    assert_eq!(prod.m[0][0], det);
    assert_eq!(prod.m[1][1], det);
    assert_eq!(prod.m[0][1], b_zero());
    assert_eq!(prod.m[1][0], b_zero());
}

#[test]
fn test_dim2_lattice_metrics() {
    let norm_q = b(3);
    let v1 = Vec2::new(b(2), b(5));
    let v2 = Vec2::new(b(-1), b(4));

    // c0^2 + norm_q * c1^2
    let n1 = dim2_lattice_norm(&v1, &norm_q);
    assert_eq!(n1, b(2 * 2 + 3 * 5 * 5)); // 79

    // v1.c0 * v2.c0 + norm_q * v1.c1 * v2.c1
    let b12 = dim2_lattice_bilinear(&v1, &v2, &norm_q);
    assert_eq!(b12, b(2 * -1 + 3 * 5 * 4)); // 58
}

#[test]
fn test_dim2_lattice_short_basis_exact() {
    let norm_q = b(1);
    let mut mat = Mat2x2::zero();

    mat.m[0][0] = b(1);
    mat.m[0][1] = b(100);
    mat.m[1][0] = b(0);
    mat.m[1][1] = b(1);

    let reduced = dim2_lattice_short_basis(&mat, &norm_q);

    let v0 = Vec2::new(reduced.m[0][0].clone(), reduced.m[1][0].clone());
    let v1 = Vec2::new(reduced.m[0][1].clone(), reduced.m[1][1].clone());

    let norm_0 = dim2_lattice_norm(&v0, &norm_q);
    let norm_1 = dim2_lattice_norm(&v1, &norm_q);

    assert!(norm_0 <= b(2), "Basis 0 is not optimally reduced");
    assert!(norm_1 <= b(2), "Basis 1 is not optimally reduced");

    // determinant check
    let (_, det_orig) = mat.inv_with_det_as_denom();
    let (_, det_red) = reduced.inv_with_det_as_denom();

    assert_eq!(det_orig.abs(), det_red.abs(), "Lattice determinant altered during Gauss reduction");
}

#[test]
fn test_dim2_lattice_short_basis_fuzzing() {
    let mut rng = StdRng::seed_from_u64(42);
    let norm_q = b(103);

    for _ in 0..100 {
        let mut mat = Mat2x2::zero();

        // random 2D basis
        let mut det = b_zero();
        while det.is_zero() {
            mat.m[0][0] = b(rng.random_range(-50..50));
            mat.m[0][1] = b(rng.random_range(-50..50));
            mat.m[1][0] = b(rng.random_range(-50..50));
            mat.m[1][1] = b(rng.random_range(-50..50));

            let (_, d) = mat.inv_with_det_as_denom();
            det = d;
        }

        let reduced = dim2_lattice_short_basis(&mat, &norm_q);

        let v0 = Vec2::new(reduced.m[0][0].clone(), reduced.m[1][0].clone());
        let v1 = Vec2::new(reduced.m[0][1].clone(), reduced.m[1][1].clone());

        let norm_0 = dim2_lattice_norm(&v0, &norm_q);
        let norm_1 = dim2_lattice_norm(&v1, &norm_q);

        assert!(norm_0 <= norm_1, "Cohen reduction did not strictly order the basis vectors");
        let (_, det_red) = reduced.inv_with_det_as_denom();
        assert_eq!(det.abs(), det_red.abs(), "Determinant volume changed");
    }
}
