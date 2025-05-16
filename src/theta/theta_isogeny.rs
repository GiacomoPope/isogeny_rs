// ===================================================================
// Compting general (2,2)-isogenies between theta structures
//
// NOTE: For the two steps before a product structure is reached, we
// need additional symplectic transforms which is controlled by the
// `hadamard` array of `bool`s. The purpose of these is to avoid null
// points (or dual null points) which have zero elements, which are
// incompatible with the doubling formula.
// ===================================================================

use fp2::fq::Fq as FqTrait;

use super::theta_point::ThetaPoint;
use super::theta_structure::ThetaStructure;
use crate::theta::theta_util::{to_hadamard, to_squared_theta};

/// Returns `0xFF..FF` when <T1, T2> are isotropic and `0x00..00` otherwise.
fn check_isotropic<Fq: FqTrait>(
    T1: &ThetaPoint<Fq>,
    T2: &ThetaPoint<Fq>,
    I: &ThetaPoint<Fq>,
) -> u32 {
    let mut ok = u32::MAX;

    let (x1, y1, z1, t1) = T1.coords();
    let (x2, y2, z2, t2) = T2.coords();
    let (A_inv, B_inv, C_inv, D_inv) = I.coords();

    ok &= (x1 * A_inv).equals(&(y1 * B_inv));
    ok &= (z1 * C_inv).equals(&(t1 * D_inv));
    ok &= (x2 * A_inv).equals(&(z2 * C_inv));
    ok &= (y2 * B_inv).equals(&(t2 * D_inv));

    ok
}

/// Given the 8-torsion above the kernel, compute the codomain of the
/// (2,2)-isogeny and the image of all points in `image_points`
/// Cost:
/// Codomain: 8S + 9M
/// Image: 4S + 4M
pub fn two_isogeny<Fq: FqTrait>(
    T1: &ThetaPoint<Fq>,
    T2: &ThetaPoint<Fq>,
    image_points: &mut [ThetaPoint<Fq>],
    hadamard: [bool; 2],
) -> (ThetaStructure<Fq>, u32) {
    // Compute the squared theta transform of both elements
    // of the kernel
    let (xA, xB, yC, yD) = T1.squared_theta();
    let (zA, tB, zC, tD) = T2.squared_theta();

    // The kernel is faulty if any coordinate of T1 or T2 is zero.
    // TODO: we don't actually have to do this, as the chain will
    // just produce (0 : 0 : 0 : 0) at the end and this is caught
    // by splitting.
    let mut ok = u32::MAX;
    ok &= !T1.has_zero_coordinate();
    ok &= !T2.has_zero_coordinate();

    // Compute the codomain coordinates
    let xAtB = xA * tB;
    let zAxB = zA * xB;
    let zCtD = zC * tD;

    let mut A = zA * xAtB;
    let mut B = tB * zAxB;
    let mut C = zC * xAtB;
    let mut D = tD * zAxB;

    // Inverses are precomputed for evaluation below
    let A_inv = xB * zCtD;
    let B_inv = xA * zCtD;
    let C_inv = D;
    let D_inv = C;

    // With the above inverses, we can now check if the kernel
    // is maximally isotropic.
    ok &= check_isotropic(
        &ThetaPoint::from_coords(&xA, &xB, &yC, &yD),
        &ThetaPoint::from_coords(&zA, &tB, &zC, &tD),
        &ThetaPoint::from_coords(&A_inv, &B_inv, &C_inv, &D_inv),
    );

    // Finish computing the codomain coordinates
    // For the penultimate case, we skip the hadamard transformation
    if hadamard[1] {
        (A, B, C, D) = to_hadamard(&A, &B, &C, &D);
    }
    let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

    // Now push through each point through the isogeny
    for P in image_points.iter_mut() {
        let (mut XX, mut YY, mut ZZ, mut TT) = P.coords();
        if hadamard[0] {
            (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
            (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);
        } else {
            (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);
        }

        XX *= A_inv;
        YY *= B_inv;
        ZZ *= C_inv;
        TT *= D_inv;

        if hadamard[1] {
            (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
        }

        P.X = XX;
        P.Y = YY;
        P.Z = ZZ;
        P.T = TT;
    }

    (codomain, ok)
}

/// Special function for the case when we are (2,2)-isogenous to a
/// product of elliptic curves. Essentially the same as above, but with
/// some small changes to deal with that a dual coordinate is now zero.
/// Computes the codomain of the (2,2)-isogeny and the image of all
/// points in `image_points`
/// Cost:
/// Codomain: 8S + 13M
/// Image: 4S + 3M
pub fn two_isogeny_to_product<Fq: FqTrait>(
    T1: &ThetaPoint<Fq>,
    T2: &ThetaPoint<Fq>,
    image_points: &mut [ThetaPoint<Fq>],
) -> (ThetaStructure<Fq>, u32) {
    // Compute the squared theta transform of both elements
    // of the kernel
    let (mut xA, mut xB, yC, yD) = T1.hadamard();
    (xA, xB, _, _) = to_squared_theta(&xA, &xB, &yC, &yD);

    let (mut zA, mut tB, mut zC, mut tD) = T2.hadamard();
    (zA, tB, zC, tD) = to_squared_theta(&zA, &tB, &zC, &tD);

    // Compute the codomain coordinates
    let zAtB = zA * tB;
    let A = xA * zAtB;
    let B = xB * zAtB;
    let C = zC * xA * tB;
    let D = tD * xB * zA;

    // Inverses are precomputed for evaluation below
    let AB = A * B;
    let CD = C * D;
    let A_inv = CD * B;
    let B_inv = CD * A;
    let C_inv = AB * D;
    let D_inv = AB * C;

    let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

    for P in image_points.iter_mut() {
        let (mut XX, mut YY, mut ZZ, mut TT) = P.coords();

        (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
        (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);

        XX *= A_inv;
        YY *= B_inv;
        ZZ *= C_inv;
        TT *= D_inv;

        P.X = XX;
        P.Y = YY;
        P.Z = ZZ;
        P.T = TT;
    }

    // The last step of the chain is not verified, as the verification is done during splitting.
    (codomain, u32::MAX)
}
