use fp2::fq::Fq as FqTrait;

use crate::elliptic::{curve::Curve, projective_point::Point};

#[derive(Clone, Copy, Debug)]
pub struct CouplePoint<Fq: FqTrait> {
    P1: Point<Fq>,
    P2: Point<Fq>,
}

impl<Fq: FqTrait> CouplePoint<Fq> {
    /// Return the pair of points at infinity: O1, O2 on E1 x E2
    pub const INFINITY: Self = Self {
        P1: Point::INFINITY,
        P2: Point::INFINITY,
    };

    /// Create a CouplePoint given a pair of points P1, P2 on E1 x E2
    pub fn new(P1: &Point<Fq>, P2: &Point<Fq>) -> Self {
        Self { P1: *P1, P2: *P2 }
    }

    /// Return the points P1, P2
    pub fn points(self) -> (Point<Fq>, Point<Fq>) {
        (self.P1, self.P2)
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for CouplePoint<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "CouplePoint with Points:\n{}\n{}", self.P1, self.P2)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct EllipticProduct<Fq: FqTrait> {
    E1: Curve<Fq>,
    E2: Curve<Fq>,
}

impl<Fq: FqTrait> EllipticProduct<Fq> {
    /// Create an EllipticProduct given a pair of elliptic curves of
    /// type Curve
    pub fn new(E1: &Curve<Fq>, E2: &Curve<Fq>) -> Self {
        Self { E1: *E1, E2: *E2 }
    }

    /// Return the pair of curves as a tuple
    pub fn curves(self) -> (Curve<Fq>, Curve<Fq>) {
        (self.E1, self.E2)
    }

    /// Addition of elements (P1, P2) and (Q1, Q2) on E1 x E2 is defined
    /// as (P1 + Q1, P2 + Q2). This function calls the add function for
    /// the pair of curves on the EllipticProduct
    pub fn add(self, C1: &CouplePoint<Fq>, C2: &CouplePoint<Fq>) -> CouplePoint<Fq> {
        let mut C3 = CouplePoint::INFINITY;
        C3.P1 = self.E1.add(&C1.P1, &C2.P1);
        C3.P2 = self.E2.add(&C1.P2, &C2.P2);
        C3
    }

    /// Doubles the pair of points (P1, P2) on E1 x E2 as ([2]P1, [2]P2)
    pub fn double(self, C: &CouplePoint<Fq>) -> CouplePoint<Fq> {
        let mut C3 = *C;
        C3.P1 = self.E1.double(&C3.P1);
        C3.P2 = self.E2.double(&C3.P2);
        C3
    }

    /// Repeatedly doubles the pair of points (P1, P2) on E1 x E2 to get
    /// ([2^n]P1, [2^n]P2)
    pub fn double_iter(self, C: &CouplePoint<Fq>, n: usize) -> CouplePoint<Fq> {
        let mut C3 = *C;
        C3.P1 = self.E1.double_iter(&C3.P1, n);
        C3.P2 = self.E2.double_iter(&C3.P2, n);
        C3
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for EllipticProduct<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "EllipticProduct with Curves:\n{}\n{}", self.E1, self.E2)
    }
}
