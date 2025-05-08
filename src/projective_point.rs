use fp2::fq::Fq;

/// Projective representation of a point (X : Y : Z)
#[derive(Clone, Copy, Debug)]
pub struct Point<Fp2: Fq> {
    X: Fp2,
    Y: Fp2,
    Z: Fp2,
}

impl<Fp2: Fq> Point<Fp2> {
    /// The point-at-infinity (neutral element of the group law).
    pub const INFINITY: Self = Self {
        X: Fq::ZERO,
        Y: Fq::ONE,
        Z: Fq::ZERO,
    };

    /// Create a new point: WARNING no check is made on the validity of the point.
    fn new(X: &Fp2, Y: &Fp2, Z: &Fp2) -> Self {
        Self {
            X: *X,
            Y: *Y,
            Z: *Z,
        }
    }

    /// Negate the point
    pub fn set_neg(&mut self) {
        self.Y.set_neg()
    }

    /// Copy rhs into self if ctl == 0xFFFFFFFF.
    /// Do nothing is ctl == 0x00000000.
    /// ctl MUST be either 0xFFFFFFFF or 0x00000000.
    pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
        self.X.set_cond(&rhs.X, ctl);
        self.Y.set_cond(&rhs.Y, ctl);
        self.Z.set_cond(&rhs.Z, ctl);
    }

    // Return 0xFFFFFFFF if self and rhs represent the same point.
    /// Otherwise, return 0x00000000.
    pub fn equals(self, rhs: &Self) -> u32 {
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

impl<Fp2: Fq> core::ops::Neg for Point<Fp2> {
    type Output = Point<Fp2>;

    #[inline(always)]
    fn neg(self) -> Point<Fp2> {
        let mut r = self;
        r.set_neg();
        r
    }
}

impl<Fp2: Fq> core::ops::Neg for &Point<Fp2> {
    type Output = Point<Fp2>;

    #[inline(always)]
    fn neg(self) -> Point<Fp2> {
        let mut r = *self;
        r.set_neg();
        r
    }
}
