use fp2::fq::Fq as FqTrait;

use super::{curve::Curve, projective_point::Point};

impl<Fq: FqTrait> Curve<Fq> {
    /// P3 <- P1 + P2
    pub fn add_into(self, P3: &mut Point<Fq>, P1: &Point<Fq>, P2: &Point<Fq>) {
        // Complete routine, to handle all edge cases:
        //   if Z1 == 0:            # P1 == inf
        //       return P2
        //   if Z2 == 0:            # P2 == inf
        //       return P1
        //   L <- Y2*Z1 - Y1*Z2
        //   T <- X2*Z1 - X1*Z2
        //   if T == 0:             # x1 == x2
        //       if L == 0:         # ... and y1 == y2: doubling case
        //           L <- 3*X1^2 + 2*A*X1*Z1 + Z1^2
        //           T <- 2*Y1*Z1
        //       else:              # ... but y1 != y2, thus P2 = -P1
        //           return inf
        //   U <- Z1*Z2*L^2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
        //   X3 <- U*T
        //   Y3 <- L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
        //   Z3 <- Z1*Z2*T^3
        //
        // Constant-time processing:
        //   Cases P1 == inf and P2 == inf are handled at the end.
        //   (L,T) are always computed for both normal and doubling cases.
        //   If P1 == -P2 then we can let T == 0 and L != 0, this will
        //   properly lead to Z3 == 0.
        //
        // Formulas from https://eprint.iacr.org/2015/1060 are faster
        // but do not cover the case when P1 - P2 is a point of order 2,
        // which can happen in all generality.
        //
        // : current formulas have cost 16M+5S; this can probably
        // be improved. Main issues to tackle:
        //   - Multiplications by A are expensive (since A can be any value)
        //   - There are three points of order 2; this makes finding
        //     complete formulas challenging.

        // T = X2*Z1 - X1*Z2
        // L = Y2*Z1 - Y1*Z2
        let x1z2 = P1.X * P2.Z;
        let x2z1 = P2.X * P1.Z;
        let mut T = x2z1 - x1z2;
        let y1z2 = P1.Y * P2.Z;
        let y2z1 = P2.Y * P1.Z;
        let mut L = y2z1 - y1z2;

        // Alternate (T,L) for doubling:
        //   Td = 2*Y1*Z1
        //   Ld = 3*X1^2 + 2*A*X1*Z1 + Z1^2
        let dbl = T.is_zero() & L.is_zero();
        let Td = (P1.Y * P1.Z).mul2();
        let x1x1 = P1.X.square();
        let z1z1 = P1.Z.square();
        let dx1z1 = (P1.X + P1.Z).square() - x1x1 - z1z1;
        let Ld = x1x1.mul3() + z1z1 + self.A * dx1z1;
        T.set_cond(&Td, dbl);
        L.set_cond(&Ld, dbl);

        // U = L^2*Z1*Z2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
        let T2 = T.square();
        let T3 = T * T2;
        let z1z2 = P1.Z * P2.Z;
        let U = (L.square() * z1z2) - ((x1z2 + x2z1 + (self.A * z1z2)) * T2);

        // X3 = U*T
        // Y3 = L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
        // Z3 = Z1*Z2*T^3
        P3.X = U * T;
        P3.Y = (L * ((x1z2 * T2) - U)) - (y1z2 * T3);
        P3.Z = z1z2 * T3;

        // Corrective action in case one of the inputs was the
        // point-at-infinity.
        let inf1 = P1.Z.is_zero();
        let inf2 = P2.Z.is_zero();
        P3.set_cond(P2, inf1);
        P3.set_cond(P1, inf2);
    }

    /// P1 <- P1 + P2
    pub fn addto(self, P1: &mut Point<Fq>, P2: &Point<Fq>) {
        let mut P3 = Point::INFINITY;
        self.add_into(&mut P3, P1, P2);
        *P1 = P3;
    }

    /// Return P1 + P2 as a new point
    pub fn add(self, P1: &Point<Fq>, P2: &Point<Fq>) -> Point<Fq> {
        let mut P3 = Point::INFINITY;
        self.add_into(&mut P3, P1, P2);
        P3
    }

    /// P3 <- P1 - P2
    pub fn sub_into(self, P3: &mut Point<Fq>, P1: &Point<Fq>, P2: &Point<Fq>) {
        let mut nP2 = *P2;
        nP2.set_neg();
        self.add_into(P3, P1, &nP2);
    }

    /// P1 <- P1 - P2
    pub fn subfrom(self, P1: &mut Point<Fq>, P2: &Point<Fq>) {
        let mut nP2 = *P2;
        nP2.set_neg();
        self.addto(P1, &nP2);
    }

    /// Return P1 - P2 as a new point
    pub fn sub(self, P1: &Point<Fq>, P2: &Point<Fq>) -> Point<Fq> {
        let mut nP2 = *P2;
        nP2.set_neg();
        self.add(P1, &nP2)
    }

    #[inline(always)]
    pub fn double_from_coords(self, X: &Fq, Y: &Fq, Z: &Fq) -> (Fq, Fq, Fq) {
        // Doubling formulas in cost 6M+6S
        // These formulas are complete.
        // Formulas from https://eprint.iacr.org/2015/1060 would be
        // more expensive, because multiplications by A are not cheap
        // in the general case.
        //
        // V <- X^2 - Z^2
        // M <- X^2 + Z^2
        // X' <- 2*Y*Z*V^2
        // Y' <- V*(M*(M + 2*A*X*Z) + 4*X^2*Z^2)
        // Z' <- 8*(Y*Z)^3
        let xx = X.square();
        let zz = Z.square();
        let dxz = ((*X) + (*Z)).square() - xx - zz;
        let dyz = ((*Y) * (*Z)).mul2();
        let v = xx - zz;
        let m = xx + zz;
        let X2 = dyz * v.square();
        let Y2 = v * ((m * (m + (self.A * dxz))) + dxz.square());
        let Z2 = dyz * dyz.square();

        (X2, Y2, Z2)
    }

    /// P3 <- 2*P1
    pub fn double_into(self, P3: &mut Point<Fq>, P1: &Point<Fq>) {
        let (X2, Y2, Z2) = self.double_from_coords(&P1.X, &P1.Y, &P1.Z);
        P3.X = X2;
        P3.Y = Y2;
        P3.Z = Z2;
    }

    // This is essentially reuse of the above, but allowing the point
    // itself to be mutated... Maybe it would be better to redesign the
    // below functions to stop the code duplication...
    pub fn double_self(self, P1: &mut Point<Fq>) {
        let (X2, Y2, Z2) = self.double_from_coords(&P1.X, &P1.Y, &P1.Z);
        P1.X = X2;
        P1.Y = Y2;
        P1.Z = Z2;
    }

    /// Return 2*P as a new point
    pub fn double(self, P: &Point<Fq>) -> Point<Fq> {
        let mut P3 = Point::INFINITY;
        self.double_into(&mut P3, P);
        P3
    }

    /// Return [2^n]*P as a new point
    pub fn double_iter(self, P: &Point<Fq>, n: usize) -> Point<Fq> {
        let mut P3 = *P;
        for _ in 0..n {
            self.double_self(&mut P3);
        }
        P3
    }

    /// P3 <- n*P
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn mul_into(self, P3: &mut Point<Fq>, P: &Point<Fq>, n: &[u8], nbitlen: usize) {
        // Montgomery ladder: see https://eprint.iacr.org/2017/212

        // We will need the complete 2*P at the end, to handle some
        // special cases of the formulas.
        let dP = self.double(P);
        let mut X0 = Fq::ONE;
        let mut Z0 = Fq::ZERO;
        let mut X1 = P.X;
        let mut Z1 = P.Z;
        let mut cc = 0u32;
        if nbitlen > 21 {
            // If n is large enough then it is worthwhile to
            // normalize the source point to affine.
            // We do not care if P = inf, since that is handled at
            // the end in the corrective steps.
            let Xp = P.X / P.Z;
            for i in (0..nbitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                Self::xadd_aff(&Xp, &X0, &Z0, &mut X1, &mut Z1);
                self.xdbl(&mut X0, &mut Z0);
                cc = ctl;
            }
        } else {
            for i in (0..nbitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                Self::xadd(&P.X, &P.Z, &X0, &Z0, &mut X1, &mut Z1);
                self.xdbl(&mut X0, &mut Z0);
                cc = ctl;
            }
        }
        Fq::condswap(&mut X0, &mut X1, cc);
        Fq::condswap(&mut Z0, &mut Z1, cc);

        // Special cases:
        //  - ladder fails if P = (0,0) (a point of order 2)
        //  - y is not reconstructed correctly if P has order 2,
        //    or if (n+1)*P = P, -P or infinity.
        let z0z = Z0.is_zero();
        let z1z = Z1.is_zero();
        let x1ex = (X1 * P.Z).equals(&(P.X * Z1));

        // (X0/Z0) is the X coordinate of P0 = n*P
        // (X1/Z1) is the X coordinate of P1 = (n + 1)*P
        // We recompute the Y coordinate of n*P (formulas from
        // Okeya and Sakurai).
        let xxzz = (P.X * X0) + (P.Z * Z0);
        let xpz0 = P.X * Z0;
        let x0zp = X0 * P.Z;
        let zz = P.Z * Z0;
        let zzdA = self.A.mul2() * zz;
        let u = (xxzz * (xpz0 + x0zp + zzdA)) - (zzdA * zz);
        let v = P.Y.mul2() * zz * Z1;
        P3.X = X0 * v;
        P3.Y = (u * Z1) - ((xpz0 - x0zp).square() * X1);
        P3.Z = Z0 * v;

        // Fix result for the special cases.
        //  P = inf                          -> inf
        //  P != inf, 2*P = inf              -> inf or P (depending on n_0)
        //  2*P != inf, P0 = inf             -> inf
        //  2*P != inf, P0 != inf, P1 = inf  -> -P
        //  2*P != inf, P0 != inf, P1 = -P   -> -2*P
        let order1 = P.Z.is_zero();
        let order2 = !order1 & P.Y.is_zero();
        let z0inf = !order1 & !order2 & z0z;
        let z1inf = !order1 & !order2 & !z0z & z1z;
        let p1mp = !order1 & !order2 & !z0z & !z1z & x1ex;

        let n_odd = ((n[0] as u32) & 1).wrapping_neg();
        P3.Z.set_cond(&Fq::ZERO, order1 | (order2 & !n_odd) | z0inf);
        P3.set_cond(P, z1inf | (order2 & n_odd));
        P3.set_cond(&dP, p1mp);
        P3.set_condneg(z1inf | p1mp);
    }

    /// Return n*P as a new point.
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn mul(self, P: &Point<Fq>, n: &[u8], nbitlen: usize) -> Point<Fq> {
        let mut P3 = Point::INFINITY;
        self.mul_into(&mut P3, P, n, nbitlen);
        P3
    }
}
