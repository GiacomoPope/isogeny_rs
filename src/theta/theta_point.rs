use crate::theta::theta_util::{to_hadamard, to_squared_theta};
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
    /// Use for initalisation, probably stupid, or at least should have
    /// a different name!
    pub const ZERO: Self = Self {
        X: Fq::ZERO,
        Y: Fq::ZERO,
        Z: Fq::ZERO,
        T: Fq::ZERO,
    };

    /// Compile time, create a new theta point from Fq elements
    pub const fn new(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
        Self {
            X: *X,
            Y: *Y,
            Z: *Z,
            T: *T,
        }
    }

    /// Create a new theta point from Fq elements
    pub fn from_coords(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
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

    /// Recover the coordinates of the element
    pub fn to_list(&self) -> [Fq; 4] {
        [self.X, self.Y, self.Z, self.T]
    }

    /// Compute the Hadamard transform of the point's coordinates
    pub fn hadamard(&self) -> (Fq, Fq, Fq, Fq) {
        to_hadamard(&self.X, &self.Y, &self.Z, &self.T)
    }

    /// Square each of the point's coordinates and then
    /// compute the hadamard transform
    pub fn squared_theta(&self) -> (Fq, Fq, Fq, Fq) {
        to_squared_theta(&self.X, &self.Y, &self.Z, &self.T)
    }

    /// Returns if any of the coordinates of `self` are zero
    pub fn has_zero_coordinate(&self) -> u32 {
        self.X.is_zero() | self.Y.is_zero() | self.Z.is_zero() | self.T.is_zero()
    }
}

/// For debugging, pretty print out the coordinates of a point
impl<Fq: FqTrait> ::std::fmt::Display for ThetaPoint<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "{}\n{}\n{}\n{}\n", self.X, self.Y, self.Z, self.T)
    }
}
