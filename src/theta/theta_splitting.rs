// ========================================================
// Compting the symplectic transform to expose the
// product structure and then compute the correct
// splitting to Montgomery curves.
// ========================================================

use fp2::fq::Fq as FqTrait;

use crate::elliptic::curve::Curve;
use crate::elliptic::point::PointX;
use crate::elliptic::product::{CouplePoint, EllipticProduct};

use super::theta_point::ThetaPoint;
use super::theta_structure::ThetaStructure;
use super::theta_util::apply_base_change;

// This function is a bit of a mess. Ultimately, we want to know whether
// given some pair of indices whether we should multiply by minus one.
// We do this by returning either 0: do nothing, or 0xFF...FF: negate
// the value, which concretely is performed with set_negcond() on the
// field element.
//
// Mathematically we have a few things to juggle. Firstly, although the
// index should really be tuples (x, y) for x,y in {0,1} we simply index
// from {0, ..., 3}. So there is first the identification of:
//
// 0 : (0, 0)
// 1 : (1, 0)
// 2 : (0, 1)
// 3 : (1, 1)
//
// The next thing we need is the dot product of these indices
// For example:
// Take i . j is the dot product, so input (x, y) = (1, 3)
// corresponds to computing:
// (1, 0) . (1, 1) = 1*1 + 0*1 = 1
//
// This means evaluation of chi means the sign is dictated by
// => (-1)^(i.j) = (-1)^1 = -1
//
// A similar thing is done for all pairs of indices below.
//
// TODO: there may be a nicer way to organise this function, but
// I couldn't find a nice closed form for ±1 from a pair (i, j)
// which i could compute on the fly without first matching from
// x,y in {0,..,3} to i,j in {(0,0)...(1,1)} (which would mean
// using a match anyway!!).
fn chi_eval(x: &usize, y: &usize) -> u32 {
    match (x, y) {
        (0, 0) => 0,
        (0, 1) => 0,
        (0, 2) => 0,
        (0, 3) => 0,
        (1, 0) => 0,
        (1, 1) => u32::MAX,
        (1, 2) => 0,
        (1, 3) => u32::MAX,
        (2, 0) => 0,
        (2, 1) => 0,
        (2, 2) => u32::MAX,
        (2, 3) => u32::MAX,
        (3, 0) => 0,
        (3, 1) => u32::MAX,
        (3, 2) => u32::MAX,
        (3, 3) => 0,
        _ => 1,
    }
}

/// For a given index (chi, i) compute the level 2,2 constants (square).
/// The purpose of this is to identify for which (chi, i) this constant
/// is zero.
fn level_22_constants_sqr<Fq: FqTrait>(null_point: &ThetaPoint<Fq>, chi: &usize, i: &usize) -> Fq {
    let mut U_constant = Fq::ZERO;
    let null_coords = null_point.list();

    for t in 0..4 {
        let mut U_it = null_coords[t] * null_coords[i ^ t];
        U_it.set_condneg(chi_eval(chi, &t));
        U_constant += U_it;
    }
    U_constant
}

/// For each possible even index compute the level 2,2 constant. Return
/// the even index for which this constant is zero. This only fails for
/// bad input in which case the whole chain would fail. Evaluates all
/// positions, and so should run in constant time.
fn identify_even_index<Fq: FqTrait>(null_point: &ThetaPoint<Fq>) -> (usize, usize) {
    const EVEN_INDICIES: [(usize, usize); 10] = [
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 2),
        (2, 0),
        (2, 1),
        (3, 0),
        (3, 3),
    ];
    // Initialise the return tuple
    let mut chi_zero = 0;
    let mut i_zero = 0;

    for (chi, i) in EVEN_INDICIES.iter() {
        let U_sqr = level_22_constants_sqr(null_point, chi, i);

        // When U_sqr is zero, U_sqr_is_zero = 0xFF...FF
        // and 0 otherwise, so we can use this as a mask
        // to select the non-zero index through the loop
        let U_sqr_is_zero = U_sqr.is_zero();
        chi_zero |= *chi as u32 & U_sqr_is_zero;
        i_zero |= *i as u32 & U_sqr_is_zero;
    }
    (chi_zero as usize, i_zero as usize)
}

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

/// We can precompute 10 different symplectic transforms which
/// correspond to each of the possible 10 even indicies which could be
/// zero. We can select the right change of basis by using the above
/// functions and then selecting the correct map accordingly.
fn compute_splitting_matrix<Fq: FqTrait>(null_point: &ThetaPoint<Fq>) -> [Fq; 16] {
    // Identity the current location of the zero
    let zero_location = identify_even_index(null_point);

    // Compute the corresponding matrix to map the zero to
    // the desired place
    // TODO: is a match like this the best thing to do in Rust??
    let M: [Fq; 16];
    let maps = splitting_maps();
    match zero_location {
        (0, 2) => M = maps[0],
        (3, 3) => M = maps[1],
        (0, 3) => M = maps[2],
        (2, 1) => M = maps[3],
        (0, 1) => M = maps[4],
        (1, 2) => M = maps[5],
        (2, 0) => M = maps[6],
        (3, 0) => M = maps[7],
        (1, 0) => M = maps[8],
        (0, 0) => M = maps[9],
        // The above locations are an exhaustive list of possible inputs, not sure how to tell rust this...
        _ => panic!("Unreachable"),
    }

    M
}

/// Map from a theta point to one which admits a splitting to elliptic
/// products. Essentially requires computing the correct splitting
/// matrix and then applying the isomorphism
pub fn splitting_isomorphism<Fq: FqTrait>(
    Th: ThetaStructure<Fq>,
    image_points: &mut [ThetaPoint<Fq>],
) -> ThetaStructure<Fq> {
    // Compute the correct splitting matrix
    let mut O0 = Th.null_point();
    let M = compute_splitting_matrix(&O0);

    // Map the Theta Structure through the symplectic transform
    apply_base_change(&mut O0, M);

    // Map the points through the symplectic transform
    for P in image_points.iter_mut() {
        apply_base_change(P, M);
    }

    ThetaStructure::new_from_point(&mut O0)
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
    image_points: &[ThetaPoint<Fq>],
    num_image_points: usize,
) -> (EllipticProduct<Fq>, Vec<CouplePoint<Fq>>) {
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
    let mut C: CouplePoint<Fq>;
    let mut couple_points: Vec<CouplePoint<Fq>> = vec![];

    for P in image_points.iter().take(num_image_points) {
        // Split to level 1
        let (P1, P2) = split_theta_point(P);
        // Compute the XPoint (X : Z) from each theta point
        let Q1X = theta_point_to_montgomery_point(&O1, &P1);
        let Q2X = theta_point_to_montgomery_point(&O2, &P2);

        // Lift these points to (X : Y : Z) on the curves
        let (Q1, _) = E3.complete_pointX(&Q1X);
        let (Q2, _) = E4.complete_pointX(&Q2X);

        // Package these points into a CouplePoint on
        // E3 x E4
        C = CouplePoint::new(&Q1, &Q2);

        // Push this into the output
        couple_points.push(C);
    }

    (E3E4, couple_points)
}
