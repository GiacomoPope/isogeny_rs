use fp2::fq::Fq as FqTrait;

use super::{basis::BasisX, curve::Curve, point::PointX};

impl<Fq: FqTrait> Curve<Fq> {
    /// Compute the x-only double of a given point and return
    /// the X-coords
    #[inline(always)]
    pub fn xdbl_coords(self, X: &Fq, Z: &Fq) -> (Fq, Fq) {
        let mut V1 = (*X + *Z).square();
        let V2 = (*X - *Z).square();
        let X_new = V1 * V2;
        V1 -= V2;
        let mut Z_new = V1;
        Z_new *= self.A24;
        Z_new += V2;
        Z_new *= V1;

        (X_new, Z_new)
    }

    /// x-only doubling formula
    #[inline(always)]
    fn xdbl(self, X: &mut Fq, Z: &mut Fq) {
        let mut V1 = (*X + *Z).square();
        let V2 = (*X - *Z).square();
        *X = V1 * V2;
        V1 -= V2;
        *Z = V1;
        *Z *= self.A24;
        *Z += V2;
        *Z *= V1;
    }

    /// x-only differential formula Note: order of arguments:
    /// (XPQ : ZPQ), (XP : ZP), (XQ : ZQ) For PQ = P - Q
    /// Sets Q  = P + Q in place
    #[inline(always)]
    fn xadd(XPQ: &Fq, ZPQ: &Fq, XP: &Fq, ZP: &Fq, XQ: &mut Fq, ZQ: &mut Fq) {
        let V1 = (*XP - *ZP) * (*XQ + *ZQ);
        let V2 = (*XP + *ZP) * (*XQ - *ZQ);
        *XQ = *ZPQ * (V1 + V2).square();
        *ZQ = *XPQ * (V1 - V2).square();
    }

    /// x-only differential formula Note: order of arguments:
    /// (XPQ : 1), (XP : ZP), (XQ : ZQ) For PQ = P - Q
    /// Sets Q  = P + Q in place
    #[inline(always)]
    fn xadd_aff(XPQ: &Fq, XP: &Fq, ZP: &Fq, XQ: &mut Fq, ZQ: &mut Fq) {
        let V1 = (*XP - *ZP) * (*XQ + *ZQ);
        let V2 = (*XP + *ZP) * (*XQ - *ZQ);
        *XQ = (V1 + V2).square();
        *ZQ = *XPQ * (V1 - V2).square();
    }

    /// P3 <- n*P, X-only variant.
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn xmul_into(self, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: &[u8], nbitlen: usize) {
        // Montgomery ladder: see https://eprint.iacr.org/2017/212
        if nbitlen == 0 {
            P3.X = Fq::ONE;
            P3.Z = Fq::ZERO;
            return;
        }

        let mut X0 = Fq::ONE;
        let mut Z0 = Fq::ZERO;
        let mut X1 = P.X;
        let mut Z1 = P.Z;
        let mut cc = 0u32;
        if nbitlen > 21 {
            // If n is large enough then it is worthwhile to
            // normalize the source point to affine.
            // If P = inf, then this sets Xp to 0; thus, the
            // output of both xdbl() and xadd_aff() has Z = 0,
            // so we correctly get the point-at-infinity at the end.
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

        // The ladder may fail if P = (0,0) (which is a point of
        // order 2) because in that case xadd() (and xadd_aff())
        // return Z = 0 systematically, so the result is considered
        // to be the point-at-infinity, which is wrong is n is odd.
        // We adjust the result in that case.
        let spec = P.X.is_zero() & !P.Z.is_zero() & ((n[0] as u32) & 1).wrapping_neg();
        P3.X = X0;
        P3.Z = Z0;
        P3.X.set_cond(&Fq::ZERO, spec);
        P3.Z.set_cond(&Fq::ONE, spec);
    }

    /// Return n*P as a new point (X-only variant).
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn xmul(self, P: &PointX<Fq>, n: &[u8], nbitlen: usize) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        self.xmul_into(&mut P3, P, n, nbitlen);
        P3
    }

    /// P3 <- (2^e)*P (X-only variant)
    fn xmul_2e_into(self, P3: &mut PointX<Fq>, P: &PointX<Fq>, e: usize) {
        let mut X = P.X;
        let mut Z = P.Z;
        for _ in 0..e {
            let mut V1 = (X + Z).square();
            let V2 = (X - Z).square();
            X = V1 * V2;
            V1 -= V2;
            Z = V1;
            Z *= self.A24;
            Z += V2;
            Z *= V1;
        }
        P3.X = X;
        P3.Z = Z;
    }

    /// Return (2^e)*P (X-only variant).
    pub fn xmul_2e(self, P: &PointX<Fq>, e: usize) -> PointX<Fq> {
        let mut Q = PointX::INFINITY;
        self.xmul_2e_into(&mut Q, P, e);
        Q
    }

    /// x-only doubling and differential addition formula
    /// Note: order of arguments:
    /// (XP : ZP), (XQ : ZQ), (XPQ: ZPQ) For PQ = P - Q
    /// Sets P = [2]P and Q = P + Q in place
    #[inline(always)]
    fn xdbladd(self, XP: &mut Fq, ZP: &mut Fq, XQ: &mut Fq, ZQ: &mut Fq, XQP: &Fq, ZQP: &Fq) {
        // TODO: I think I could just replace P, Q in-place rather than
        // define new mutable elements of Fp2
        let mut t0 = *XP + *ZP;
        let mut t1 = *XP - *ZP;
        let mut X2P = t0.square();
        let mut t2 = *XQ - *ZQ;
        let mut XPQ = *XQ + *ZQ;
        t0 *= t2;
        let mut Z2P = t1.square();
        t1 *= XPQ;
        t2 = X2P - Z2P;
        X2P *= Z2P;
        XPQ = self.A24 * t2;
        let mut ZPQ = t0 - t1;
        Z2P = XPQ + Z2P;
        XPQ = t0 + t1;
        Z2P *= t2;
        ZPQ = ZPQ.square();
        XPQ = XPQ.square();
        ZPQ *= *XQP;
        XPQ *= *ZQP;

        // Modify in place
        *XP = X2P;
        *ZP = Z2P;
        *XQ = XPQ;
        *ZQ = ZPQ;
    }

    /// Return P + n*Q, X-only variant given the x-only basis x(P), x(Q) and x(P - Q).
    /// Integer `n` is encoded as unsigned little-endian, with length `nbitlen bits`.
    /// Bits beyond that length are ignored.
    pub fn three_point_ladder(self, B: &BasisX<Fq>, n: &[u8], nbitlen: usize) -> PointX<Fq> {
        if nbitlen == 0 {
            return B.P();
        }

        // Extract out the coordinates from the basis
        let (mut X0, mut Z0) = B.Q().coords();
        let (mut X1, mut Z1) = B.P().coords();
        let (mut X2, mut Z2) = B.PQ().coords();

        let mut cc = 0u32;
        for i in 0..nbitlen {
            let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
            Fq::condswap(&mut X1, &mut X2, ctl ^ cc);
            Fq::condswap(&mut Z1, &mut Z2, ctl ^ cc);
            self.xdbladd(&mut X0, &mut Z0, &mut X2, &mut Z2, &X1, &Z1);
            cc = ctl;
        }
        Fq::condswap(&mut X1, &mut X2, cc);
        Fq::condswap(&mut Z1, &mut Z2, cc);

        PointX::new(&X1, &Z1)
    }
}
