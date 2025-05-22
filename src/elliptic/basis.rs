use super::point::PointX;
use fp2::traits::Fp as FpTrait;

/// A x-only basis of x(P), x(Q) and x(P - Q)
#[derive(Clone, Copy, Debug)]
pub struct BasisX<Fq: FpTrait> {
    pub P: PointX<Fq>,
    pub Q: PointX<Fq>,
    pub PQ: PointX<Fq>,
}

impl<Fq: FpTrait> BasisX<Fq> {
    /// Create a basis given the x-coordinates of x(P), x(Q) and x(P - Q)
    pub const fn from_x_coords(xP: &Fq, xQ: &Fq, xPQ: &Fq) -> Self {
        let P = PointX::from_x_coord(xP);
        let Q = PointX::from_x_coord(xQ);
        let PQ = PointX::from_x_coord(xPQ);
        Self { P, Q, PQ }
    }

    /// Set the basis given the points P, Q and PQ
    pub fn from_points(P: &PointX<Fq>, Q: &PointX<Fq>, PQ: &PointX<Fq>) -> Self {
        Self {
            P: *P,
            Q: *Q,
            PQ: *PQ,
        }
    }

    pub fn from_slice(basis: &[PointX<Fq>]) -> Self {
        // TODO: should we return errors instead?
        if basis.len() != 3 {
            return Self {
                P: PointX::INFINITY,
                Q: PointX::INFINITY,
                PQ: PointX::INFINITY,
            };
        }
        Self {
            P: basis[0],
            Q: basis[1],
            PQ: basis[2],
        }
    }

    /// Return the array of points [P, Q, PQ]
    pub fn to_array(&self) -> [PointX<Fq>; 3] {
        [self.P, self.Q, self.PQ]
    }
}
