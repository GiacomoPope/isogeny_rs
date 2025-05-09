use crate::theta::theta_point::ThetaPoint;
use crate::theta::theta_util::{to_hadamard, to_squared_coords, to_squared_theta};
use fp2::fq::Fq as FqTrait;
// ========================================================
// Functions for working with ThetaStructures
// ========================================================

/// Theta Structure
#[derive(Clone, Copy, Debug)]
pub struct ThetaStructure<Fq: FqTrait> {
    null_point: ThetaPoint<Fq>,
    arithmetic_precom: [Fq; 8],
}

impl<Fq: FqTrait> ThetaStructure<Fq> {
    /// Given the coordinates of a null point, create a null point and
    /// precompute 8 Fp2 elements which are used for doubling and isogeny
    /// computations.
    pub fn new_from_coords(X: &Fq, Z: &Fq, U: &Fq, V: &Fq) -> Self {
        let null_point = ThetaPoint::new(X, Z, U, V);
        Self {
            null_point,
            arithmetic_precom: ThetaStructure::precomputation(&null_point),
        }
    }

    /// Given a null point, store the null point and precompute 8 Fp2
    /// elements which are used for doubling and isogeny computations.
    pub fn new_from_point(null_point: &ThetaPoint<Fq>) -> Self {
        Self {
            null_point: *null_point,
            arithmetic_precom: ThetaStructure::precomputation(null_point),
        }
    }

    /// Return the null point of the ThetaStructure
    pub fn null_point(self) -> ThetaPoint<Fq> {
        self.null_point
    }

    /// For doubling and also computing isogenies, we need the following
    /// constants, which we can precompute once for each ThetaStructure.
    /// Cost: 14M + 5S
    #[inline]
    pub fn precomputation(O0: &ThetaPoint<Fq>) -> [Fq; 8] {
        let (a, b, c, d) = O0.coords();
        let (aa, bb, cc, dd) = to_squared_coords(&a, &b, &c, &d);

        // Compute projectively a/b = a^2*c*d, etc.
        let ab = a * b;
        let cd = c * d;
        let x0 = ab * cd;
        let y0 = aa * cd;
        let aab = aa * b;
        let z0 = aab * d;
        let t0 = aab * c;

        // Compute projectively A^2/B^2 = A^4*C^2*D^2, etc.
        let (AA, BB, CC, DD) = to_hadamard(&aa, &bb, &cc, &dd);
        let A4 = AA.square();
        let AABB = AA * BB;
        let CCDD = CC * DD;
        let X0 = AABB * CCDD;
        let Y0 = A4 * CCDD;
        let A4BB = A4 * BB;
        let Z0 = A4BB * DD;
        let T0 = A4BB * CC;

        [x0, y0, z0, t0, X0, Y0, Z0, T0]
    }

    /// Given a point P, compute it's double [2]P in place.
    /// Cost 8S + 8M
    #[inline(always)]
    pub fn set_double_self(self, P: &mut ThetaPoint<Fq>) {
        let (mut xp, mut yp, mut zp, mut tp) = P.squared_theta();

        // Compute temp. coordinates, 8S and 3M
        xp = self.arithmetic_precom[4] * xp.square();
        yp = self.arithmetic_precom[5] * yp.square();
        zp = self.arithmetic_precom[6] * zp.square();
        tp = self.arithmetic_precom[7] * tp.square();

        // Compute the final coordinates, 3M
        let (mut X, mut Y, mut Z, mut T) = to_hadamard(&xp, &yp, &zp, &tp);
        X *= self.arithmetic_precom[0];
        Y *= self.arithmetic_precom[1];
        Z *= self.arithmetic_precom[2];
        T *= self.arithmetic_precom[3];

        *P = ThetaPoint::new(&X, &Y, &Z, &T);
    }

    /// Compute [2] * self
    #[inline]
    pub fn double_point(self, P: &ThetaPoint<Fq>) -> ThetaPoint<Fq> {
        let mut P2 = *P;
        self.set_double_self(&mut P2);
        P2
    }

    /// Compute [2^n] * self
    #[inline]
    pub fn double_iter(self, P: &ThetaPoint<Fq>, n: usize) -> ThetaPoint<Fq> {
        let mut R = *P;
        for _ in 0..n {
            self.set_double_self(&mut R)
        }
        R
    }
}
