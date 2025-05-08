use crate::point::PointX;
use fp2::fq::Fq;

/// Curve y^2 = x^3 + A*x^2 + x, for a given constant A
/// (special case of a Montgomery curve).
#[derive(Clone, Copy, Debug)]
pub struct Curve<T: Fq> {
    A: T,   // A
    A24: T, // (A+2)/4
}

impl<Fp2: Fq> Curve<Fp2> {
    /// Create a new curve instance, with the provided constant.
    pub fn new(A: &Fp2) -> Self {
        // We check that the curve is not singular, i.e. A^2 != 4.
        let a = *A;
        assert!(a.equals(&<Fp2>::TWO) == 0);
        assert!((a + <Fp2>::TWO).is_zero() == 0);

        Self {
            A: a,
            A24: (a + <Fp2>::TWO).half().half(),
        }
    }

    /// x-only differential formula with cubical normalisation
    /// (XP : ZP), (XQ : ZQ), (1 : ixPQ) For PQ = P - Q
    #[inline(always)]
    fn cubical_xadd(xP: &Fp2, zP: &Fp2, xQ: &Fp2, zQ: &Fp2, ixPQ: &Fp2) -> (Fp2, Fp2) {
        let V1 = (*xP - *zP) * (*xQ + *zQ);
        let V2 = (*xP + *zP) * (*xQ - *zQ);
        let X = (*ixPQ) * (V1 + V2).square();
        let Z = (V1 - V2).square();
        (X, Z)
    }

    /// x-only doubling and differential addition formula
    /// Assuming P - Q has been normalised.
    /// Note: order of arguments:
    /// (XP : ZP), (XQ : ZQ), (1: ixPQ) For PQ = P - Q
    /// Sets P = [2]P and Q = P + Q in place
    #[inline(always)]
    fn cubical_xdbladd(self, XP: &mut Fp2, ZP: &mut Fp2, XQ: &mut Fp2, ZQ: &mut Fp2, iXQP: &Fp2) {
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
        XPQ *= *iXQP;

        // Modify in place
        *XP = X2P;
        *ZP = Z2P;
        *XQ = XPQ;
        *ZQ = ZPQ;
    }

    /// Given xP, xQ xPQ with
    /// xP = x(P), xQ = x(Q) and xPQ = x(P - Q) sets:
    /// nP <- [n]P
    /// nPQ <- [n]P + Q
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    fn cubical_ladder_into(
        self,
        nP: &mut PointX<Fp2>,
        nPQ: &mut PointX<Fp2>,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        n: &[u8],
        nbitlen: usize,
        div_by_two: bool,
    ) {
        // Batch invert to obtain ixP, ixQ and ixPQ
        let mut inverses = [*xP, *xQ, *xPQ];
        <Fp2>::batch_invert(&mut inverses);
        let ixP = inverses[0];
        let mut ixQ = inverses[1];
        let mut ixPQ = inverses[2];

        // Set initial points for the ladder computation
        let mut xS0 = <Fp2>::ONE;
        let mut zS0 = <Fp2>::ZERO;

        let mut xS1 = *xP;
        let mut zS1 = <Fp2>::ONE;

        let mut xR;
        let mut zR;

        let mut xT = *xQ;
        let mut zT = <Fp2>::ONE;

        // If we want to divide by two, we skip the bottom bit in the ladder
        // this is used for the even pairings where the last factor of two is
        // computed using a translate.
        let start = div_by_two as usize;

        let mut cc = 0u32;
        for i in (start..nbitlen).rev() {
            // First compute R = S0 + S1
            (xR, zR) = Self::cubical_xadd(&xS0, &zS0, &xS1, &zS1, &ixP);

            // Compute [2]Si and T + Si for i = {0, 1} depending on the bit of n
            let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
            <Fp2>::condswap(&mut xS0, &mut xS1, ctl ^ cc);
            <Fp2>::condswap(&mut zS0, &mut zS1, ctl ^ cc);
            <Fp2>::condswap(&mut ixQ, &mut ixPQ, ctl ^ cc);
            self.cubical_xdbladd(&mut xS0, &mut zS0, &mut xT, &mut zT, &ixQ);

            // Update Sj to be R
            xS1 = xR;
            zS1 = zR;

            // Save swap variable
            cc = ctl;
        }
        // Perform a final swap
        <Fp2>::condswap(&mut xS0, &mut xS1, cc);
        <Fp2>::condswap(&mut zS0, &mut zS1, cc);

        // Set the variables to return S0, T = [n]P, [n]P + Q
        *nP = PointX::new(&xS0, &zS0);
        *nPQ = PointX::new(&xT, &zT);
    }

    /// Given x(P), x(Q) and x(P - Q) computes [n]P and [n]P + Q for Tate pairing.
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn cubical_ladder(
        self,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        n: &[u8],
        nbitlen: usize,
        div_by_two: bool,
    ) -> (PointX<Fp2>, PointX<Fp2>) {
        let mut nP = PointX::INFINITY;
        let mut nPQ = PointX::INFINITY;
        self.cubical_ladder_into(&mut nP, &mut nPQ, xP, xQ, xPQ, n, nbitlen, div_by_two);
        (nP, nPQ)
    }

    /// Given xP, xQ xPQ with
    /// xP = x(P), xQ = x(Q) and xPQ = x(P - Q) sets:
    /// nP <- [2^e]P
    /// nPQ <- [2^e]P + Q
    fn cubical_ladder_2exp_into(
        self,
        nP: &mut PointX<Fp2>,
        nPQ: &mut PointX<Fp2>,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        e: usize,
    ) {
        // Batch invert to obtain ixP, ixQ and ixPQ
        let mut inverses = [*xQ, *xPQ];
        <Fp2>::batch_invert(&mut inverses);
        let ixQ = inverses[0];
        let ixPQ = inverses[1];

        let mut X0 = *xP;
        let mut Z0 = <Fp2>::ONE;
        let mut X1 = *xQ;
        let mut Z1 = <Fp2>::ONE;

        // Compute P + Q
        (X1, Z1) = Self::cubical_xadd(&X0, &Z0, &X1, &Z1, &ixPQ);

        // Compute [2^e]P and [2^e]P + Q
        for _ in 0..e {
            self.cubical_xdbladd(&mut X0, &mut Z0, &mut X1, &mut Z1, &ixQ);
        }

        // Set the variables of [2^e]P, [2^e]P + Q
        *nP = PointX::new(&X0, &Z0);
        *nPQ = PointX::new(&X1, &Z1);
    }

    /// Given x(P), x(Q) and x(P - Q) computes the points
    /// [2^e]P and [2^e]P + Q for the Tate pairing.
    pub fn cubical_ladder_2exp(
        self,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        e: usize,
    ) -> (PointX<Fp2>, PointX<Fp2>) {
        let mut nP = PointX::INFINITY;
        let mut nPQ = PointX::INFINITY;
        self.cubical_ladder_2exp_into(&mut nP, &mut nPQ, xP, xQ, xPQ, e);
        (nP, nPQ)
    }

    /// Given x(P), x(Q) and x(P - Q) computes the non-reduced tate pairing
    fn tate_pairing_not_reduced(
        self,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        n: &[u8],
        nbitlen: usize,
    ) -> (Fp2, Fp2) {
        let is_even: bool = n[0] & 1 == 0;
        let (mut nP, mut nPQ) = self.cubical_ladder(xP, xQ, xPQ, n, nbitlen, is_even);
        if is_even {
            nPQ = nPQ.translate(nP);
            nP = nP.translate(nP);
        }
        (nPQ.Z, nP.X)
    }

    /// Given x(P), x(Q) and x(P - Q) computes the non-reduced tate
    /// pairing using points of order 2^e
    fn tate_pairing_not_reduced_2exp(self, xP: &Fp2, xQ: &Fp2, xPQ: &Fp2, e: usize) -> (Fp2, Fp2) {
        let (mut nP, mut nPQ) = self.cubical_ladder_2exp(xP, xQ, xPQ, e - 1);
        nPQ = nPQ.translate(nP);
        nP = nP.translate(nP);
        (nPQ.Z, nP.X)
    }

    /// Compute x^(p^2 - 1) = x^(p + 1)(p - 1) to reduce the tate pairing for E / Fp^2
    /// where Fp^2 has modulus x^2 + 1
    fn reduce_tate_pairing(self, num: Fp2, den: Fp2, d: &[u8], dbitlen: usize) -> Fp2 {
        // We can compute the pth power with a conjugate \pi(x) = x^p
        let num_p = num.conjugate();
        let den_p = den.conjugate();

        // First we compute ePQ^(p-1) using Frobenius x^p / x = x^(p-1)
        let ePQ = (num_p * den) / (den_p * num);

        // Now handle the power of (p + 1) / d with supplied exponent
        ePQ.pow(d, dbitlen)
    }

    /// Given x(P), x(Q) and x(P + Q) computes the reduced tate
    /// pairing
    /// d encodes in little endian the value (p + 1) // 2^e
    pub fn tate_pairing(
        self,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        n: &[u8],
        nbitlen: usize,
        d: &[u8],
        dbitlen: usize,
    ) -> Fp2 {
        // First compute the non-reduced tate pairing which is returned as num / den
        let (num, den) = self.tate_pairing_not_reduced(xP, xQ, xPQ, n, nbitlen);

        // Efficiently compute ePQ^(p^2 - 1)
        self.reduce_tate_pairing(num, den, d, dbitlen)
    }

    /// Given x(P), x(Q) and x(P - Q) computes the reduced tate
    /// pairing using points of order 2^e
    /// d encodes in little endian the value (p + 1) // 2^e
    pub fn tate_pairing_2exp(
        self,
        xP: &Fp2,
        xQ: &Fp2,
        xPQ: &Fp2,
        e: usize,
        d: &[u8],
        dbitlen: usize,
    ) -> Fp2 {
        // First compute the non-reduced tate pairing
        let (num, den) = self.tate_pairing_not_reduced_2exp(xP, xQ, xPQ, e);

        // Efficiently compute ePQ^(p^2 - 1)
        self.reduce_tate_pairing(num, den, d, dbitlen)
    }

    /// Given x(P), x(Q) and x(P - Q) computes the Weil pairing
    pub fn weil_pairing(self, xP: &Fp2, xQ: &Fp2, xPQ: &Fp2, n: &[u8], nbitlen: usize) -> Fp2 {
        let (e1_num, e1_den) = self.tate_pairing_not_reduced(xP, xQ, xPQ, n, nbitlen);
        let (e2_num, e2_den) = self.tate_pairing_not_reduced(xQ, xP, xPQ, n, nbitlen);
        (e1_num * e2_den) / (e2_num * e1_den)
    }

    /// Given x(P), x(Q) and x(P - Q) computes the weil pairing
    /// using points of order = 2^e
    pub fn weil_pairing_2exp(self, xP: &Fp2, xQ: &Fp2, xPQ: &Fp2, e: usize) -> Fp2 {
        let (e1_num, e1_den) = self.tate_pairing_not_reduced_2exp(xP, xQ, xPQ, e);
        let (e2_num, e2_den) = self.tate_pairing_not_reduced_2exp(xQ, xP, xPQ, e);
        (e1_num * e2_den) / (e2_num * e1_den)
    }
}
