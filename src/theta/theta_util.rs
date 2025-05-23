use fp2::traits::Fp2 as FqTrait;

use super::theta_point::ThetaPoint;

/// Given four elements of Fq, compute the hadamard transform using
/// recursive addition.
/// Cost: 8a
#[inline(always)]
pub fn to_hadamard<Fq: FqTrait>(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
    let t1 = (*X) + (*Y);
    let t2 = (*X) - (*Y);
    let t3 = (*Z) + (*T);
    let t4 = (*Z) - (*T);

    let A = t1 + t3;
    let B = t2 + t4;
    let C = t1 - t3;
    let D = t2 - t4;
    (A, B, C, D)
}

/// Given four elements of Fq, first square each coordinate
/// Cost: 4S
#[inline(always)]
pub fn to_squared_coords<Fq: FqTrait>(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
    let XX = X.square();
    let YY = Y.square();
    let ZZ = Z.square();
    let TT = T.square();

    (XX, YY, ZZ, TT)
}

/// Given four elements of Fq, first square each coordinate and
/// then compute the hadamard transform
/// Cost: 4S, 8a
#[inline(always)]
pub fn to_squared_theta<Fq: FqTrait>(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
    let (XX, YY, ZZ, TT) = to_squared_coords(X, Y, Z, T);
    to_hadamard(&XX, &YY, &ZZ, &TT)
}

/// Apply the base change described by M on a ThetaPoint in-place
/// Cost: 16M
#[inline]
pub fn apply_base_change<Fq: FqTrait>(P: &mut ThetaPoint<Fq>, M: &[Fq; 16]) {
    let (x, y, z, t) = P.coords();
    P.X = M[0] * x + M[1] * y + M[2] * z + M[3] * t;
    P.Y = M[4] * x + M[5] * y + M[6] * z + M[7] * t;
    P.Z = M[8] * x + M[9] * y + M[10] * z + M[11] * t;
    P.T = M[12] * x + M[13] * y + M[14] * z + M[15] * t;
}
