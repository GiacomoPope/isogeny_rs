use super::point::PointX;
use fp2::fq::Fq as FqTrait;

/// A x-only basis of x(P), x(Q) and x(P - Q)
#[derive(Clone, Copy, Debug)]
pub struct BasisX<Fq: FqTrait> {
    basis: [PointX<Fq>; 3],
}

impl<Fq: FqTrait> BasisX<Fq> {
    /// Create a basis given the x-coordinates of x(P), x(Q) and x(P - Q)
    pub const fn from_x_coords(xP: &Fq, xQ: &Fq, xPQ: &Fq) -> Self {
        let P = PointX::from_x_coord(xP);
        let Q = PointX::from_x_coord(xQ);
        let PQ = PointX::from_x_coord(xPQ);
        Self { basis: [P, Q, PQ] }
    }

    /// Set the basis given an array [P, Q, PQ]
    pub fn from_array(basis: [PointX<Fq>; 3]) -> Self {
        Self { basis: basis }
    }

    /// Return the array of points [P, Q, PQ]
    pub fn to_array(&self) -> [PointX<Fq>; 3] {
        self.basis
    }

    /// Return the point P from the basis.
    #[inline(always)]
    pub fn P(self) -> PointX<Fq> {
        self.basis[0]
    }

    /// Return the point Q from the basis.
    #[inline(always)]
    pub fn Q(self) -> PointX<Fq> {
        self.basis[1]
    }

    /// Return the point P - Q from the basis.
    #[inline(always)]
    pub fn PQ(self) -> PointX<Fq> {
        self.basis[2]
    }
}
