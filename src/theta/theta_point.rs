use crate::theta::theta_util::{to_hadamard, to_hadamard_square};
use fp2::traits::Fp2 as FqTrait;

// ========================================================
// Functions for working with ThetaPoints
// ========================================================

/// Theta Point Struct
#[derive(Clone, Copy, Debug)]
pub struct ThetaPoint<Fq: FqTrait> {
    pub X: Fq,
    pub Y: Fq,
    pub Z: Fq,
    pub T: Fq,
}

impl<Fq: FqTrait> ThetaPoint<Fq> {
    /// Compile time, create a new theta point from Fq elements
    pub const fn new(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
        Self {
            X: *X,
            Y: *Y,
            Z: *Z,
            T: *T,
        }
    }

    /// Recover the coordinates of the element
    pub fn coords(&self) -> (Fq, Fq, Fq, Fq) {
        (self.X, self.Y, self.Z, self.T)
    }

    /// Compute the Hadamard transform of the point's coordinates
    pub fn hadamard(&self) -> (Fq, Fq, Fq, Fq) {
        to_hadamard(&self.X, &self.Y, &self.Z, &self.T)
    }

    /// Square each of the point's coordinates and then
    /// compute the hadamard transform
    pub fn hadamard_square(&self) -> (Fq, Fq, Fq, Fq) {
        to_hadamard_square(&self.X, &self.Y, &self.Z, &self.T)
    }

    /// Square each coordinate then perform a hadamard transform after optionally
    /// having applied a Hadamard transform first.
    pub fn cond_hadamard_square_hadamard(&self, hadamard: bool) -> (Fq, Fq, Fq, Fq) {
        if hadamard {
            let (a, b, c, d) = self.hadamard();
            to_hadamard_square(&a, &b, &c, &d)
        } else {
            self.hadamard_square()
        }
    }

    /// Returns if any of the coordinates of `self` are zero
    pub fn has_zero_coordinate(&self) -> u32 {
        self.X.is_zero() | self.Y.is_zero() | self.Z.is_zero() | self.T.is_zero()
    }
}

/// Default element used for initialisation
impl<Fq: FqTrait> Default for ThetaPoint<Fq> {
    fn default() -> Self {
        Self {
            X: Fq::ZERO,
            Y: Fq::ZERO,
            Z: Fq::ZERO,
            T: Fq::ZERO,
        }
    }
}

/// Create a new theta point from a tuple of Fq elements
impl<Fq: FqTrait> From<(Fq, Fq, Fq, Fq)> for ThetaPoint<Fq> {
    fn from((X, Y, Z, T): (Fq, Fq, Fq, Fq)) -> Self {
        Self { X, Y, Z, T }
    }
}

/// Cast the coordinates of the ThetaPoint into a list
impl<Fq: FqTrait> From<&ThetaPoint<Fq>> for [Fq; 4] {
    fn from(P: &ThetaPoint<Fq>) -> Self {
        [P.X, P.Y, P.Z, P.T]
    }
}

/// For debugging, pretty print out the coordinates of a point
impl<Fq: FqTrait> ::std::fmt::Display for ThetaPoint<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "{}\n{}\n{}\n{}\n", self.X, self.Y, self.Z, self.T)
    }
}
