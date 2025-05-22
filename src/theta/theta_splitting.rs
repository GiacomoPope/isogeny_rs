// ========================================================
// Compting the symplectic transform to expose the
// product structure and then compute the correct
// splitting to Montgomery curves.
// ========================================================

use fp2::traits::Fp2 as FqTrait;

use crate::elliptic::curve::Curve;
use crate::elliptic::point::PointX;
use crate::utilities::ct::ct_u32_eq;

use super::elliptic_product::{EllipticProduct, ProductPoint};
use super::theta_point::ThetaPoint;
use super::theta_structure::ThetaStructure;
use super::theta_util::apply_base_change;

fn splitting_maps<Fq: FqTrait>() -> [[Fq; 16]; 10] {
    [
        [
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::MINUS_ONE,
            <Fq>::ZERO,
        ],
        [
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
        ],
        [
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::MINUS_ONE,
        ],
        [
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
        ],
        [
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::MINUS_ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
        ],
        [
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
            <Fq>::ZERO,
            <Fq>::ONE,
            <Fq>::ZERO,
        ],
        [
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::ONE,
        ],
        [
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
        ],
        [
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
            <Fq>::MINUS_ONE,
            <Fq>::ONE,
        ],
        [
            <Fq>::ONE,
            <Fq>::ZETA,
            <Fq>::ONE,
            <Fq>::ZETA,
            <Fq>::ONE,
            <Fq>::MINUS_ZETA,
            <Fq>::ONE,
            <Fq>::ZETA,
            <Fq>::ONE,
            <Fq>::ZETA,
            <Fq>::ONE,
            <Fq>::MINUS_ZETA,
            <Fq>::ONE,
            <Fq>::ZETA,
            <Fq>::ONE,
            <Fq>::ZETA,
        ],
    ]
}

/// Conditionally set the matrix M1 to M2 when ctl = `0xFF..FF`.
fn set_cond_matrix<Fq: FqTrait>(M1: &mut [Fq], M2: &[Fq], ctl: u32) {
    for i in 0..M1.len() {
        M1[i].set_cond(&M2[i], ctl);
    }
}

/// We can precompute 10 different symplectic transforms which
/// correspond to each of the possible 10 even indicies which could be
/// zero. We can select the right change of basis by using the above
/// functions and then selecting the correct map accordingly. Return this
/// matrix, together with a u32 with `0xFF..FF` indicating success and
/// `0x00..00` failure.
fn compute_splitting_matrix<Fq: FqTrait>(null_point: &ThetaPoint<Fq>) -> ([Fq; 16], u32) {
    const EVEN_INDICIES: [[usize; 2]; 10] = [
        [0, 2],
        [3, 3],
        [0, 3],
        [2, 1],
        [0, 1],
        [1, 2],
        [2, 0],
        [3, 0],
        [1, 0],
        [0, 0],
    ];
    const CHI_EVAL: [[i8; 4]; 4] = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]];

    // M is the matrix we return to compute the change of basis to find the split product.
    let mut M = [Fq::ZERO; 16];

    // M are the maps we will iterate over.
    let maps: [[Fq; 16]; 10] = splitting_maps();

    // We must iterate through all 10 even indicies and select the correct splitting matrix.
    let mut t1: Fq;
    let mut t2: Fq;
    let mut U_sqr: Fq;
    let null_coords = null_point.to_list();

    // The number of zeros found should be exactly one.
    let mut count: u32 = 0;

    for i in 0..10 {
        U_sqr = Fq::ZERO;
        for j in 0..4 {
            t1 = null_coords[j];
            t2 = null_coords[j ^ EVEN_INDICIES[i][1]];
            t1 *= t2;

            // If chi(i, t) is +1 we want ctl to be 0x00..00
            // If chi(i, t) is -1 we want ctl to be 0xFF..FF
            let ctl = (CHI_EVAL[EVEN_INDICIES[i][0]][j] >> 1) as u32;
            t1.set_condneg(ctl);
            U_sqr += t1;
        }

        let ctl = U_sqr.is_zero();
        count += ctl & 1;
        set_cond_matrix(&mut M, &maps[i], ctl);
    }
    let ok = ct_u32_eq(count, 1);
    (M, ok)
}

/// Map from a theta point to one which admits a splitting to elliptic
/// products. Essentially requires computing the correct splitting
/// matrix and then applying the isomorphism. Return the correct structure
/// together with a u32 with `0xFF..FF` indicating success and
/// `0x00..00` failure (determined by the ability to find a unique matrix).
pub fn splitting_isomorphism<Fq: FqTrait>(
    Th: ThetaStructure<Fq>,
    image_points: &mut [ThetaPoint<Fq>],
) -> (ThetaStructure<Fq>, u32) {
    // Compute the correct splitting matrix
    let mut O0 = Th.null_point();
    let (M, ok) = compute_splitting_matrix(&O0);

    // Map the Theta Structure through the symplectic transform
    apply_base_change(&mut O0, &M);

    // Map the points through the symplectic transform
    for P in image_points.iter_mut() {
        apply_base_change(P, &M);
    }

    (ThetaStructure::new_from_point(&O0), ok)
}

/// Given a Theta point in the correct representation, compute two
/// dimension 1 theta points.
/// Algorithm from:
/// Models of Kummer lines and Galois representation,
/// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
fn split_theta_point<Fq: FqTrait>(P: &ThetaPoint<Fq>) -> ((Fq, Fq), (Fq, Fq)) {
    let (a, b, _, d) = P.coords();

    let P1 = (a, b);
    let P2 = (b, d);

    (P1, P2)
}

/// Given a dimension one null theta point, compute the corresponding
/// elliptic curve in the Montgomery model by recovering the Montgomery
/// coefficient A
/// Algorithm from:
/// Models of Kummer lines and Galois representation,
/// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
fn null_point_to_montgomery_curve<Fq: FqTrait>(O0: &(Fq, Fq)) -> Curve<Fq> {
    let (a, b) = O0;

    let aa = a.square();
    let bb = b.square();

    let T1 = aa + bb;
    let T2 = aa - bb;

    let A = -(T1.square() + T2.square()) / (T1 * T2);

    Curve::new(&A)
}

/// Given a dimension one theta point, compute the corresponding
/// elliptic curve point on the Kummer line (X : Z)
/// Algorithm from:
/// Models of Kummer lines and Galois representation,
/// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
fn theta_point_to_montgomery_point<Fq: FqTrait>(O0: &(Fq, Fq), P: &(Fq, Fq)) -> PointX<Fq> {
    let (a, b) = *O0;
    let (U, V) = *P;

    let X = a * V + b * U;
    let Z = a * V - b * U;

    // TODO: rather than use this Type we could directly lift here instead...
    // exposing PointX to be public was the easiest way to keep everything from
    // eccore the same, which will help in the future
    PointX::new(&X, &Z)
}

/// Given a ThetaStructure and set of points in the compatible form,
/// compute the product of elliptic curves and affine points on these
/// curves.
pub fn split_to_product<Fq: FqTrait>(
    Th: &ThetaStructure<Fq>,
    eval_points: &[ThetaPoint<Fq>],
    img_points: &mut [ProductPoint<Fq>],
) -> EllipticProduct<Fq> {
    // First we take the domain theta null point and
    // split this to two level-1 theta null points
    let null_point = Th.null_point();
    let (O1, O2) = split_theta_point(&null_point);

    // Compute Montgomery curve from dimension one
    // null points
    let E3 = null_point_to_montgomery_curve(&O1);
    let E4 = null_point_to_montgomery_curve(&O2);
    let E3E4 = EllipticProduct::new(&E3, &E4);

    // Now compute points on E3 x E4
    for (i, P) in eval_points.iter().enumerate() {
        // Split to level 1
        let (P1, P2) = split_theta_point(P);
        // Compute the XPoint (X : Z) from each theta point
        let Q1X = theta_point_to_montgomery_point(&O1, &P1);
        let Q2X = theta_point_to_montgomery_point(&O2, &P2);

        // Lift these points to (X : Y : Z) on the curves
        let (Q1, _) = E3.lift_pointx(&Q1X);
        let (Q2, _) = E4.lift_pointx(&Q2X);

        // Package these points into a CouplePoint on E3 x E4
        img_points[i] = ProductPoint::new(&Q1, &Q2);
    }

    E3E4
}
