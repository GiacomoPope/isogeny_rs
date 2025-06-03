use fp2::traits::Fp as FpTrait;
use rand_core::{CryptoRng, RngCore};

/// Special x-only representation of a point (or a pair of points,
/// since two Y coordinates may match a given X).
#[derive(Clone, Copy, Debug)]
pub struct PointX<Fq: FpTrait> {
    pub X: Fq,
    pub Z: Fq,
}

impl<Fq: FpTrait> PointX<Fq> {
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
    pub fn coords(&self) -> (Fq, Fq) {
        (self.X, self.Z)
    }

    /// Point at infinity is of the form (X : 0)
    /// Returns 0xFFFFFFFF if Z == 0 and 0 otherwise
    pub fn is_zero(&self) -> u32 {
        self.Z.is_zero()
    }

    /// Returns the affine `x` coordinate of a point
    /// x = X / Z
    pub fn x(&self) -> Fq {
        self.X / self.Z
    }

    /// Return 0xFFFFFFFF if self and rhs represent the same point.
    /// Otherwise, return 0x00000000.
    pub fn equals(&self, rhs: &PointX<Fq>) -> u32 {
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

    pub fn batch_normalise(points: &mut [Self]) {
        let mut zs: Vec<Fq> = points.iter().map(|point| point.Z).collect();
        Fq::batch_invert(&mut zs);
        for (i, P) in points.iter_mut().enumerate() {
            P.X *= zs[i];
            P.Z = Fq::ONE;
        }
    }

    /// Return a new random X only point. 
    pub fn rand_point<R: CryptoRng + RngCore>(rng: &mut R) -> PointX<Fq> {
        let x = Fq::rand(rng);
        PointX::from_x_coord(&x)
    }
}

impl<Fq: FpTrait> ::std::fmt::Display for PointX<Fq> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "PointX: ({} : {})", self.X, self.Z)
    }
}
