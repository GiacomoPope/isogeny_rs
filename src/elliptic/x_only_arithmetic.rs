use fp2::fq::Fq as FqTrait;

use crate::utilities::le_bytes::encode_to_odd_binary;

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
    pub fn xdbl(self, X: &mut Fq, Z: &mut Fq) {
        let mut V1 = (*X + *Z).square();
        let V2 = (*X - *Z).square();
        *X = V1 * V2;
        V1 -= V2;
        *Z = V1;
        *Z *= self.A24;
        *Z += V2;
        *Z *= V1;
    }

    /// x-only doubling, set `R` to the value of \[2\]P
    #[inline]
    fn set_xdouble(self, xR: &mut PointX<Fq>) {
        self.xdbl(&mut xR.X, &mut xR.Z);
    }

    /// Return the value [2]P
    #[inline]
    fn xdouble(self, xP: &PointX<Fq>) -> PointX<Fq> {
        let mut xR = *xP;
        self.set_xdouble(&mut xR);
        xR
    }

    /// Return the value [2^n]P
    #[inline]
    fn xdouble_iter(self, xP: &PointX<Fq>, n: usize) -> PointX<Fq> {
        let mut xR = *xP;
        for _ in 0..n {
            self.set_xdouble(&mut xR);
        }
        xR
    }

    /// x-only differential formula Note: order of arguments:
    /// (XPQ : ZPQ), (XP : ZP), (XQ : ZQ) For PQ = P - Q
    /// Sets Q  = P + Q in place
    #[inline(always)]
    pub fn xadd(XPQ: &Fq, ZPQ: &Fq, XP: &Fq, ZP: &Fq, XQ: &mut Fq, ZQ: &mut Fq) {
        let V1 = (*XP - *ZP) * (*XQ + *ZQ);
        let V2 = (*XP + *ZP) * (*XQ - *ZQ);
        *XQ = *ZPQ * (V1 + V2).square();
        *ZQ = *XPQ * (V1 - V2).square();
    }

    /// x-only differential addition with PointX type, sets `R` to x(P + Q) given x(P)
    /// x(Q) and x(P - Q) as `PointX<Fq>`.
    #[inline]
    fn xdiff_add_into(R: &mut PointX<Fq>, xP: &PointX<Fq>, xQ: &PointX<Fq>, xPmQ: &PointX<Fq>) {
        R.X = xQ.X;
        R.Z = xQ.Z;
        Self::xadd(&xPmQ.X, &xPmQ.Z, &xP.X, &xP.Z, &mut R.X, &mut R.Z);
    }

    /// Return x(P + Q) given x(P), x(Q) and x(P - Q) as `PointX<Fq>`.
    #[inline]
    fn xdiff_add(xP: &PointX<Fq>, xQ: &PointX<Fq>, xPmQ: &PointX<Fq>) -> PointX<Fq> {
        let mut R = PointX::INFINITY;
        Self::xdiff_add_into(&mut R, xP, xQ, xPmQ);
        R
    }

    /// x-only differential formula Note: order of arguments:
    /// (XPQ : 1), (XP : ZP), (XQ : ZQ) For PQ = P - Q
    /// Sets Q  = P + Q in place
    #[inline(always)]
    pub fn xadd_aff(XPQ: &Fq, XP: &Fq, ZP: &Fq, XQ: &mut Fq, ZQ: &mut Fq) {
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

    /// Return (2^e)*P (x-only variant).
    pub fn xmul_2e(self, P: &PointX<Fq>, e: usize) -> PointX<Fq> {
        let mut Q = PointX::INFINITY;
        self.xmul_2e_into(&mut Q, P, e);
        Q
    }

    /// Return (2^e)*R for R in [P, Q, P - Q] (x-only variant).
    pub fn basis_xmul_2e(self, B: &BasisX<Fq>, e: usize) -> BasisX<Fq> {
        let P = self.xmul_2e(&B.P(), e);
        let Q = self.xmul_2e(&B.Q(), e);
        let PQ = self.xmul_2e(&B.PQ(), e);
        BasisX::from_array([P, Q, PQ])
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
    /// Integer `n` is encoded as unsigned little-endian, with length `nbitlen` bits.
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

    /// Helper function for `ladder_biscalar` which re-encodes the scalars `a` and `b` into
    /// two bit-values `(s0, s1)` and a bit-vector `r`.
    fn ladder_biscalar_reencoding(
        a: &[u8],
        b: &[u8],
        a_bitlen: usize,
        b_bitlen: usize,
    ) -> (usize, usize, Vec<u8>) {
        // Compute the max bit-length of a and b to set the length of r
        let k = usize::max(a_bitlen, b_bitlen);

        // First we derive the bit s0 and s1, which are set from
        // if a is even and b is odd then s0, s1 = (1, 0) otherwise (0, 1)
        let a_is_odd = ((a[0] & 1) as u32).wrapping_neg();
        let b_is_odd = ((b[0] & 1) as u32).wrapping_neg();
        let mask: u32 = !a_is_odd & b_is_odd;
        let mut s0: u32 = mask;
        let mut s1: u32 = !mask;

        // We now want to ensure both scalars are odd for the encoding, which
        // means we want to subtract 1 when the scalar is even, and 0 otherwise.
        // We also need to encode each scalar into binary, so we do these two
        // operations together.
        let mut a_bits = vec![0u8; k + 1];
        let mut b_bits = vec![0u8; k + 1];
        encode_to_odd_binary(&mut a_bits, a, a_bitlen);
        encode_to_odd_binary(&mut b_bits, b, b_bitlen);
        let ab_bits: &[Vec<u8>] = &[a_bits, b_bits];

        // Create a new vector of length 2*k
        let mut r = vec![0; 2 * k];
        for i in 0..k {
            let s0_index = (s0 & 1) as usize;
            let s1_index = (s1 & 1) as usize;

            // Set bits in r from bits from the a and b scalars
            r[2 * i] = ab_bits[s0_index][i] ^ ab_bits[s0_index][i + 1];
            r[2 * i + 1] = ab_bits[s1_index][i] ^ ab_bits[s1_index][i + 1];

            // Swap s0 and s1 when r[2*i + 1] is 1.
            // As we know s0 and s1 are 0 or -1 then when r[2*i + 1] = 0
            // mask = 0 and si ^ mask = si. When r[2*i + 1] = 1 then mask
            // is 0xFF..FF and si ^ mask = sj as it will flip all the bits.
            let mask = (r[2 * i + 1] as u32).wrapping_neg();
            s0 ^= mask;
            s1 ^= mask;
        }

        ((s0 & 1) as usize, (s1 & 1) as usize, r)
    }

    /// Return [a]P + [b]*Q, X-only variant given the x-only basis x(P), x(Q) and x(P - Q).
    /// The integers `a` and `b` are encoded as unsigned little-endian.
    pub fn ladder_biscalar(
        self,
        B: &BasisX<Fq>,
        a: &[u8],
        b: &[u8],
        a_bitlen: usize,
        b_bitlen: usize,
    ) -> PointX<Fq> {
        // Encode the scalars a and b into the appropriate form for the biscalar ladder.
        let (s0, s1, r) = Self::ladder_biscalar_reencoding(a, b, a_bitlen, b_bitlen);
        let k = r.len() >> 1;

        let mut T: [PointX<Fq>; 3] = [PointX::INFINITY; 3];
        let mut R: [PointX<Fq>; 3] = [PointX::INFINITY; 3];
        T[0] = B.P();
        T[1] = B.Q();

        R[1] = T[s0];
        R[2] = T[s1];

        // Compute the difference points for T, R
        let mut D1 = R[1];
        let mut D2 = R[2];
        R[2] = Self::xdiff_add(&R[1], &R[2], &B.PQ());
        let mut F1 = R[2];
        let mut F2 = B.PQ();

        // Main ladder loop, compute [a]P + [b]Q
        for i in (0..k).rev() {
            // Compute indices.
            let r1 = r[2 * i] as usize;
            let r2 = r[2 * i + 1] as usize;
            let h = r1 + r2;

            // Compute new T values and swap differential points conditionally.
            T[0] = R[h & 1];
            T[1] = R[2];
            T[0] = self.xdouble(&T[h >> 1]);
            T[1] = R[r2];
            T[2] = R[r2 + 1];
            PointX::condswap(&mut D1, &mut D2, (r2 as u32).wrapping_neg());
            T[1] = Self::xdiff_add(&T[1], &T[2], &D1);
            T[2] = Self::xdiff_add(&R[0], &R[2], &F1);
            PointX::condswap(&mut F1, &mut F2, ((h & 1) as u32).wrapping_neg());

            // Update R values from T values.
            R = T;
        }

        // When a and b are both even we want R[0], if a and b are both odd
        // we want R[2], otherwise we want R[1]
        let index = ((a[0] & 1) + (b[0] & 1)) as usize;
        R[index]
    }
}
