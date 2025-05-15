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

    #[inline]
    pub fn condswap(P: &mut Self, Q: &mut Self, ctl: u32) {
        Fq::condswap(&mut P.X, &mut Q.X, ctl);
        Fq::condswap(&mut P.Z, &mut Q.Z, ctl);
    }
}

impl<Fq: FqTrait> ::std::fmt::Display for PointX<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "PointX: ({} : {})", self.X, self.Z)
    }
}
