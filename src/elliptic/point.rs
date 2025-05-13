use fp2::fq::Fq as FqTrait;

/// Special x-only representation of a point (or a pair of points,
/// since two Y coordinates may match a given X).
#[derive(Clone, Copy, Debug)]
pub struct PointX<Fq: FqTrait> {
    pub X: Fq,
    pub Z: Fq,
}

impl<Fq: FqTrait> PointX<Fq> {
    /// The neutral point of the group.
    pub const INFINITY: Self = Self {
        X: <Fq>::ONE,
        Z: <Fq>::ZERO,
    };

    /// Create a point from coordinates. WARNING: no check is made on the point.
    pub fn new(X: &Fq, Z: &Fq) -> Self {
        Self { X: *X, Z: *Z }
    }

    pub const fn from_x_coord(X: &Fq) -> Self {
        Self { X: *X, Z: Fq::ONE }
    }

    /// Recover the (X : Z) coordinates of a PointX
    pub fn coords(self) -> (Fq, Fq) {
        (self.X, self.Z)
    }

    /// Point at infinity is of the form (X : 0)
    /// Returns 0xFFFFFFFF if Z == 0 and 0 otherwise
    pub fn is_zero(self) -> u32 {
        self.Z.is_zero()
    }

    /// Returns the affine `x` coordinate of a point
    /// x = X / Z
    pub fn x(self) -> Fq {
        self.X / self.Z
    }

    /// Return 0xFFFFFFFF if self and rhs represent the same point.
    /// Otherwise, return 0x00000000.
    pub fn equals(self, rhs: &PointX<Fq>) -> u32 {
        let inf1 = self.is_zero();
        let inf2 = rhs.is_zero();
        let e = (self.X * rhs.Z).equals(&(rhs.X * self.Z));
        (inf1 & inf2) | (!inf1 & !inf2 & e)
    }

    /// Affine translation by a two torsion point needed for even degree
    /// Tate pairings
    pub fn translate(self, T: Self) -> Self {
        let (A, B) = T.coords();
        // When we translates three things can happen.
        // - If T = (X : 0) then the translation of P is P
        // - If T = (0 : Z) with Z != 0, the translation of P = (X : Z)
        //   is (Z : X)
        // - Otherwise T = (A : B) and the translation of P = (X : Z) is
        // (A X - B Z : B X - A Z)

        // Compute generic values
        let mut X_new = (A * self.X) - (B * self.Z);
        let mut Z_new = (B * self.X) - (A * self.Z);

        // If A is zero, we should return (Z : X) instead
        let A_is_zero = A.is_zero();
        X_new.set_cond(&self.Z, A_is_zero);
        Z_new.set_cond(&self.X, A_is_zero);

        // If B is zero, we should return (X : Z) instead
        let B_is_zero = B.is_zero();
        X_new.set_cond(&self.X, B_is_zero);
        Z_new.set_cond(&self.Z, B_is_zero);

        Self::new(&X_new, &Z_new)
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for PointX<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "PointX: ({} : {})", self.X, self.Z)
    }
}
