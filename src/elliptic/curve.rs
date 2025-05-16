use super::{point::PointX, projective_point::Point};
use fp2::fq::Fq as FqTrait;
use rand_core::{CryptoRng, RngCore};

/// Curve y^2 = x^3 + A*x^2 + x, for a given constant A
/// (special case of a Montgomery curve).
#[derive(Clone, Copy, Debug)]
pub struct Curve<Fq: FqTrait> {
    pub A: Fq,   // A
    pub A24: Fq, // (A+2)/4
}

impl<Fq: FqTrait> Curve<Fq> {
    /// Create a new curve instance, with the provided constant.
    pub fn new(A: &Fq) -> Self {
        // We check that the curve is not singular, i.e. A^2 != 4.
        let a = *A;
        assert!(a.equals(&<Fq>::TWO) == 0);
        assert!((a + <Fq>::TWO).is_zero() == 0);

        Self {
            A: a,
            A24: (a + <Fq>::TWO).half().half(),
        }
    }

    /// Compute the j-invariant of the curve.
    pub fn j_invariant(self) -> Fq {
        let mut j = self.A.square();
        let mut t1 = Fq::ONE; // This should be C^2, but C = 1
        let mut t0 = Fq::TWO; // This should be 2C^2
        t0 = j - t0;
        t0 -= t1;
        j = t0 - t1;
        t1.set_square();
        j *= t1;
        t0.set_mul4();
        t1 = t0.square();
        t0 *= t1;
        t0.set_mul4();
        j.set_invert();
        j *= t0;

        j
    }

    /// Given a potential x-coordinate, determine if it's a valid point
    /// on the curve
    pub fn is_on_curve(self, x: &Fq) -> u32 {
        let x = *x;
        let mut y = x + self.A; // y = x + A
        y *= x; // y = x^2 + A*x
        y += <Fq>::ONE; // x^2 + A*x + 1
        y *= x; // y = x^3 + A*x^2 + x
        y.is_square()
    }

    /// Given the x-coordinate of a point, lift it to a projective point
    // TODO: should this return an error for the sqrt?
    pub fn lift_point(self, x: &Fq) -> Point<Fq> {
        let x = *x;
        let mut y = x + self.A; // y = x + A
        y *= x; // y = x^2 + A*x
        y += <Fq>::ONE; // x^2 + A*x + 1
        y *= x; // y = x^3 + A*x^2 + x
        y.set_sqrt();
        Point::new(&x, &y, &<Fq>::ONE)
    }

    /// Given the x-coordinates of x(P), x(Q) and x(P - Q) lift the points
    /// onto the curve <P, Q>.
    pub fn lift_basis(self, xP: &Fq, xQ: &Fq, xPQ: &Fq) -> (Point<Fq>, Point<Fq>) {
        let P = self.lift_point(xP);

        // Okeya-Sakurai algorithm to recover Q.Y without a sqrt
        let mut v2 = (*xP) + (*xQ);
        let mut v3 = (*xP) - (*xQ);
        v3.set_square();
        v3 *= *xPQ;
        let mut v1 = self.A.mul2();
        v2 += v1;
        let mut v4 = (*xP) * (*xQ);
        v4 += <Fq>::ONE;
        v2 *= v4;
        v2 -= v1;
        let y = v3 - v2;
        v1 = P.Y + P.Y;
        let x = (*xQ) * v1;
        let Q = Point::new(&x, &y, &v1);

        (P, Q)
    }

    /// Complete an X-only point into a full point;
    /// On error, P3 is set to the point-at-infinity.
    fn lift_pointx_into(self, P3: &mut Point<Fq>, P: &PointX<Fq>) -> u32 {
        let XZ = P.X * P.Z;
        let V = (P.X + P.Z).square() + ((self.A - Fq::TWO) * XZ);
        P3.X = XZ;
        P3.Y = V * XZ;
        let ok = P3.Y.set_sqrt();
        P3.Z = P.Z.square();

        // Set to inf on error.
        P3.Z.set_cond(&Fq::ZERO, !ok);

        ok
    }

    /// Complete an X-only point into a full point;
    /// On error, the output point is set to the point-at-infinity.
    pub fn lift_pointx(self, P: &PointX<Fq>) -> (Point<Fq>, u32) {
        let mut P3 = Point::INFINITY;
        let ok = self.lift_pointx_into(&mut P3, P);
        (P3, ok)
    }

    /// Set P to a random curve point.
    pub fn set_rand_point<R: CryptoRng + RngCore>(self, rng: &mut R, P: &mut Point<Fq>) {
        // This function cannot actually return the point-at-infinity;
        // this is not a problem as long as the curve order is larger
        // than 2^128.
        P.Z = Fq::ONE;
        loop {
            P.X.set_rand(rng);
            P.Y = (((P.X + self.A) * P.X) * P.X) + P.X;
            if P.Y.legendre() >= 0 {
                P.Y.set_sqrt();

                // Randomly chooses the square root to use.
                let mut tmp = [0u8; 1];
                rng.fill_bytes(&mut tmp);
                let ctl = 0u32.wrapping_sub((tmp[0] as u32) & 1);
                P.Y.set_condneg(ctl);
                return;
            }
        }
    }

    /// Return a new random curve point.
    pub fn rand_point<R: CryptoRng + RngCore>(self, rng: &mut R) -> Point<Fq> {
        let mut P = Point::INFINITY;
        self.set_rand_point(rng, &mut P);
        P
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for Curve<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Elliptic Curve: y^2 = x^3 + ({})*x^2 + x", self.A)
    }
}
