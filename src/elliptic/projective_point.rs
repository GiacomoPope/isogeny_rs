use fp2::traits::Fp as FpTrait;

use super::point::PointX;

/// Projective representation of a point (X : Y : Z)
#[derive(Clone, Copy, Debug)]
pub struct Point<Fq: FpTrait> {
    pub X: Fq,
    pub Y: Fq,
    pub Z: Fq,
}

impl<Fq: FpTrait> Point<Fq> {
    /// The point-at-infinity (neutral element of the group law).
    pub const INFINITY: Self = Self {
        X: Fq::ZERO,
        Y: Fq::ONE,
        Z: Fq::ZERO,
    };

    /// Create a new point.
    /// WARNING no check is made on the validity of the point.
    pub fn new(X: &Fq, Y: &Fq, Z: &Fq) -> Self {
        Self {
            X: *X,
            Y: *Y,
            Z: *Z,
        }
    }

    /// Create a new point from affine coordinates (x, y).
    /// WARNING no check is made on the validity of the point.
    pub fn new_xy(X: &Fq, Y: &Fq) -> Self {
        Self {
            X: *X,
            Y: *Y,
            Z: Fq::ONE,
        }
    }

    /// Get the (x,y) affine coordinates. For the point-at-infinity,
    /// this returns (0,0).
    pub fn to_xy(&self) -> (Fq, Fq) {
        let t = self.Z.invert();
        (self.X * t, self.Y * t)
    }

    /// Returns the X and Z coordinates of the projective point
    pub fn to_xz(&self) -> (Fq, Fq) {
        (self.X, self.Z)
    }

    /// Returns the X and Z coordinates of the projective point
    pub fn to_pointx(&self) -> PointX<Fq> {
        PointX::new(&self.X, &self.Z)
    }

    /// Negate the point in place.
    pub fn set_neg(&mut self) {
        self.Y.set_neg()
    }

    /// Copy rhs into self if ctl == `0xFFFFFFFF`.
    /// Do nothing is ctl == `0x00000000`.
    /// ctl MUST be either `0xFFFFFFFF` or `0x00000000`.
    pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
        self.X.set_cond(&rhs.X, ctl);
        self.Y.set_cond(&rhs.Y, ctl);
        self.Z.set_cond(&rhs.Z, ctl);
    }

    /// Negate this point if ctl == `0xFFFFFFFF`.
    /// Do nothing is ctl == `0x00000000`.
    /// ctl MUST be either `0xFFFFFFFF` or `0x00000000`.
    pub fn set_condneg(&mut self, ctl: u32) {
        self.Y.set_condneg(ctl);
    }

    /// Return `0xFFFFFFFF` if self is the point-at-infinity, `0x00000000` otherwise.
    pub fn is_zero(&self) -> u32 {
        self.Z.is_zero()
    }

    // Return `0xFFFFFFFF` if self and rhs represent the same point.
    /// Otherwise, return `0x00000000`.
    pub fn equals(&self, rhs: &Self) -> u32 {
        // P1 == P2 if and only if:
        //    P1 == inf AND P2 == inf
        //  OR:
        //    P1 != inf AND P2 != inf AND X1*Z2 = X2*Z1 AND Y1*Z2 = Y2*Z1
        let lz = self.Z.is_zero();
        let rz = rhs.Z.is_zero();
        let vx = (self.X * rhs.Z).equals(&(rhs.X * self.Z));
        let vy = (self.Y * rhs.Z).equals(&(rhs.Y * self.Z));
        (lz & rz) | (!lz & !rz & vx & vy)
    }
}

impl<Fq: FpTrait> core::ops::Neg for Point<Fq> {
    type Output = Point<Fq>;

    #[inline(always)]
    fn neg(self) -> Point<Fq> {
        let mut r = self;
        r.set_neg();
        r
    }
}

impl<Fq: FpTrait> core::ops::Neg for &Point<Fq> {
    type Output = Point<Fq>;

    #[inline(always)]
    fn neg(self) -> Point<Fq> {
        let mut r = *self;
        r.set_neg();
        r
    }
}

impl<Fq: FpTrait> ::std::fmt::Display for Point<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        if self.is_zero() == u32::MAX {
            write!(f, "Point: (0 : 1 : 0)")
        } else {
            let (x, y) = self.to_xy();
            write!(f, "Point: ({} : {} : 1)", x, y)
        }
    }
}
