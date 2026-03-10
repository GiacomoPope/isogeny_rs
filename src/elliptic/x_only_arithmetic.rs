use fp2::traits::Fp as FpTrait;

use crate::utilities::bn::{
    bn_div4_vartime, bn_from_le_bytes, bn_is_zero_vartime, bn_lt_vartime, bn_set_div2_vartime,
    bn_sub_into_vartime,
};
use crate::utilities::le_bytes::encode_to_odd_binary;

use super::{basis::BasisX, curve::Curve, point::PointX};

impl<Fq: FpTrait> Curve<Fq> {
    // ============================================================================
    // x-only doubling methods
    // ============================================================================

    /// P <- [2]*P (x-only variant) in place
    #[inline(always)]
    pub fn set_xdbl(&self, P: &mut PointX<Fq>) {
        let mut V1 = (P.X + P.Z).square();
        let V2 = (P.X - P.Z).square();
        P.X = V1 * V2;
        V1 -= V2;
        P.Z = V1;
        P.Z *= self.A24;
        P.Z += V2;
        P.Z *= V1;
    }

    /// Return [2]*P (x-only variant).
    pub fn xdbl(&self, P: &PointX<Fq>) -> PointX<Fq> {
        let mut Q = *P;
        self.set_xdbl(&mut Q);
        Q
    }

    /// P <- (2^e)*P (x-only variant)
    pub fn set_xdbl_iter(&self, P: &mut PointX<Fq>, n: usize) {
        for _ in 0..n {
            self.set_xdbl(P);
        }
    }

    /// Return (2^e)*P (x-only variant).
    pub fn xdbl_iter(&self, P: &PointX<Fq>, n: usize) -> PointX<Fq> {
        let mut Q = *P;
        self.set_xdbl_iter(&mut Q, n);
        Q
    }

    /// Compute [2]P in place using projective (A + 2) / 4 = (A24 : C24)
    /// Cost: 2S + 4M
    #[inline(always)]
    pub fn set_xdbl_proj(P: &mut PointX<Fq>, A24: &Fq, C24: &Fq) {
        let mut t0 = P.X + P.Z;
        t0.set_square();
        let mut t1 = P.X - P.Z;
        t1.set_square();
        let t2 = t0 - t1;
        t1 *= *C24;
        P.X = t0 * t1;
        t0 = t2 * (*A24);
        t0 += t1;
        P.Z = t0 * t2;
    }

    /// Compute \[2^n\]P in place using projective (A + 2) / 4 = (A24 : C24).
    /// Cost: n * (2S + 4M)
    #[inline]
    pub fn set_xdbl_proj_iter(P: &mut PointX<Fq>, A24: &Fq, C24: &Fq, n: usize) {
        for _ in 0..n {
            Self::set_xdbl_proj(P, A24, C24);
        }
    }

    // ============================================================================
    // x-only differential addition methods
    // ============================================================================

    #[inline(always)]
    pub fn xdiff_add_into(P: &PointX<Fq>, Q: &mut PointX<Fq>, PmQ: &PointX<Fq>) {
        let V1 = (P.X - P.Z) * (Q.X + Q.Z);
        let V2 = (P.X + P.Z) * (Q.X - Q.Z);
        Q.X = PmQ.Z * (V1 + V2).square();
        Q.Z = PmQ.X * (V1 - V2).square();
    }

    pub fn xdiff_add(P: &PointX<Fq>, Q: &PointX<Fq>, PmQ: &PointX<Fq>) -> PointX<Fq> {
        let mut R = *Q;
        Self::xdiff_add_into(P, &mut R, PmQ);
        R
    }

    /// x-only differential formula
    #[inline(always)]
    pub fn xdiff_add_aff_into(P: &PointX<Fq>, Q: &mut PointX<Fq>, xPmQ: &Fq) {
        let V1 = (P.X - P.Z) * (Q.X + Q.Z);
        let V2 = (P.X + P.Z) * (Q.X - Q.Z);
        Q.X = (V1 + V2).square();
        Q.Z = *xPmQ * (V1 - V2).square();
    }

    /// Return x(P + Q) given x(P), x(Q) and x(P - Q) as `PointX<Fq>`.
    #[inline]
    pub fn xdiff_add_aff(P: &PointX<Fq>, Q: &PointX<Fq>, xPmQ: &Fq) -> PointX<Fq> {
        let mut R = *Q;
        Self::xdiff_add_aff_into(P, &mut R, xPmQ);
        R
    }

    // ============================================================================
    // x-only double and add methods
    // ============================================================================

    /// x-only doubling and differential addition formula
    /// Note: order of arguments:
    /// (XP : ZP), (XQ : ZQ), (XPQ: ZPQ) For PQ = P - Q
    /// Sets P = [2]P and Q = P + Q in place
    #[inline(always)]
    pub fn xdbladd_into(&self, P: &mut PointX<Fq>, Q: &mut PointX<Fq>, PmQ: &PointX<Fq>) {
        let t0 = P.X + P.Z;
        let t1 = P.X - P.Z;
        let X2P_sq = t0.square();
        let Z2P_sq = t1.square();

        let t2 = X2P_sq - Z2P_sq;
        P.X = X2P_sq * Z2P_sq;
        P.Z = t2 * (Z2P_sq + self.A24 * t2);

        let t0_Q = t0 * (Q.X - Q.Z);
        let t1_Q = t1 * (Q.X + Q.Z);

        Q.X = PmQ.Z * (t0_Q + t1_Q).square();
        Q.Z = PmQ.X * (t0_Q - t1_Q).square();
    }

    #[inline(always)]
    pub fn xdbladd_aff_into(&self, P: &mut PointX<Fq>, Q: &mut PointX<Fq>, xPmQ: &Fq) {
        let t0 = P.X + P.Z;
        let t1 = P.X - P.Z;
        let X2P_sq = t0.square();
        let Z2P_sq = t1.square();

        let t2 = X2P_sq - Z2P_sq;
        P.X = X2P_sq * Z2P_sq;
        P.Z = t2 * (Z2P_sq + self.A24 * t2);

        let t0_Q = t0 * (Q.X - Q.Z);
        let t1_Q = t1 * (Q.X + Q.Z);

        Q.X = (t0_Q + t1_Q).square();
        Q.Z = *xPmQ * (t0_Q - t1_Q).square();
    }

    /// Return x(P + Q) given x(P), x(Q) and x(P - Q) as `PointX<Fq>`.
    #[inline]
    pub fn xdbladd(
        &self,
        P: &PointX<Fq>,
        Q: &PointX<Fq>,
        PmQ: &PointX<Fq>,
    ) -> (PointX<Fq>, PointX<Fq>) {
        let mut R = *P;
        let mut S = *Q;
        self.xdbladd_into(&mut R, &mut S, PmQ);
        (R, S)
    }

    // ============================================================================
    // x-only Montgomery ladder for scalar multiplication
    // ============================================================================

    /// P3 <- n*P, x-only variant.
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn xmul_into(&self, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: &[u8], nbitlen: usize) {
        // Montgomery ladder: see https://eprint.iacr.org/2017/212
        if nbitlen == 0 {
            P3.X = Fq::ONE;
            P3.Z = Fq::ZERO;
            return;
        }

        let mut X0 = PointX::INFINITY;
        let mut X1 = *P;
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
                PointX::condswap(&mut X0, &mut X1, ctl ^ cc);
                self.xdbladd_aff_into(&mut X0, &mut X1, &Xp);
                cc = ctl;
            }
        } else {
            for i in (0..nbitlen).rev() {
                let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                PointX::condswap(&mut X0, &mut X1, ctl ^ cc);
                self.xdbladd_into(&mut X0, &mut X1, P);
                cc = ctl;
            }
        }
        PointX::condswap(&mut X0, &mut X1, cc);

        // The ladder may fail if P = (0,0) (which is a point of
        // order 2) because in that case xadd() (and xadd_aff())
        // return Z = 0 systematically, so the result is considered
        // to be the point-at-infinity, which is wrong is n is odd.
        // We adjust the result in that case.
        let spec = P.X.is_zero() & !P.Z.is_zero() & ((n[0] as u32) & 1).wrapping_neg();
        *P3 = X0;
        P3.set_cond(&PointX::INFINITY, spec);
    }

    /// Return n*P as a new point (x-only variant).
    /// Integer n is encoded as unsigned little-endian, with length
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn xmul(&self, P: &PointX<Fq>, n: &[u8], nbitlen: usize) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        self.xmul_into(&mut P3, P, n, nbitlen);
        P3
    }

    /// P3 <- n*P, x-only variant.
    /// Integer n is represented as a u64 and the scalar n is assumed to be public
    /// nbitlen bits. Bits beyond that length are ignored.
    pub fn set_xmul_u64_vartime(&self, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: u64) {
        // Handle small cases.
        match n {
            0 => {
                *P3 = PointX::INFINITY;
            }
            1 => {
                *P3 = *P;
            }
            2 => *P3 = self.xdbl(P),
            _ => {
                let nbitlen = (64 - n.leading_zeros()) as usize;

                let mut X0 = PointX::INFINITY;
                let mut X1 = *P;
                let mut cc = 0u32;
                if nbitlen > 21 {
                    // If n is large enough then it is worthwhile to
                    // normalize the source point to affine.
                    // If P = inf, then this sets Xp to 0; thus, the
                    // output of both xdbl() and xadd_aff() has Z = 0,
                    // so we correctly get the point-at-infinity at the end.
                    let Xp = P.X / P.Z;
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n >> i) as u32) & 1).wrapping_neg();
                        PointX::condswap(&mut X0, &mut X1, ctl ^ cc);
                        self.xdbladd_aff_into(&mut X0, &mut X1, &Xp);
                        cc = ctl;
                    }
                } else {
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n >> i) as u32) & 1).wrapping_neg();
                        PointX::condswap(&mut X0, &mut X1, ctl ^ cc);
                        self.xdbladd_into(&mut X0, &mut X1, P);
                        cc = ctl;
                    }
                }
                PointX::condswap(&mut X0, &mut X1, cc);

                // The ladder may fail if P = (0,0) (which is a point of
                // order 2) because in that case xadd() (and xadd_aff())
                // return Z = 0 systematically, so the result is considered
                // to be the point-at-infinity, which is wrong is n is odd.
                // We adjust the result in that case.
                let spec = P.X.is_zero() & !P.Z.is_zero() & (((n & 1) as u32) & 1).wrapping_neg();
                *P3 = X0;
                P3.set_cond(&PointX::INFINITY, spec);
            }
        }
    }

    /// Return n*P as a new point (x-only variant).
    /// Integer n is encoded as a u64 which is assumed to be a public value.
    pub fn xmul_u64_vartime(&self, P: &PointX<Fq>, n: u64) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        self.set_xmul_u64_vartime(&mut P3, P, n);
        P3
    }

    /// Return (2^e)*R for R in [P, Q, P - Q] (x-only variant).
    pub fn basis_double_iter(&self, B: &BasisX<Fq>, e: usize) -> BasisX<Fq> {
        let P = self.xdbl_iter(&B.P, e);
        let Q = self.xdbl_iter(&B.Q, e);
        let PQ = self.xdbl_iter(&B.PQ, e);
        BasisX::from_points(&P, &Q, &PQ)
    }

    // ============================================================================
    // x-only specialised ladders
    // ============================================================================

    /// Return P + n*Q, x-only variant given the x-only basis x(P), x(Q) and x(P - Q).
    /// Integer `n` is encoded as unsigned little-endian, with length `nbitlen` bits.
    /// Bits beyond that length are ignored.
    pub fn three_point_ladder(&self, B: &BasisX<Fq>, n: &[u8], nbitlen: usize) -> PointX<Fq> {
        if nbitlen == 0 {
            return B.P;
        }

        // Extract out the coordinates from the basis
        let mut X0 = B.Q;
        let mut X1 = B.P;
        let mut X2 = B.PQ;

        let mut cc = 0u32;
        for i in 0..nbitlen {
            let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
            PointX::condswap(&mut X1, &mut X2, ctl ^ cc);
            self.xdbladd_into(&mut X0, &mut X2, &X1);
            cc = ctl;
        }
        PointX::condswap(&mut X1, &mut X2, cc);
        X1
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
        let k = a_bitlen.max(b_bitlen);

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

    /// Return [a]P + [b]*Q, x-only variant given the x-only basis x(P), x(Q) and x(P - Q).
    /// The integers `a` and `b` are encoded as unsigned little-endian.
    pub fn ladder_biscalar(
        &self,
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
        T[0] = B.P;
        T[1] = B.Q;

        R[1] = T[s0];
        R[2] = T[s1];

        // Compute the difference points for T, R
        let D1 = R[1];
        let D2 = R[2];
        R[2] = Self::xdiff_add(&R[1], &R[2], &B.PQ);
        let F1 = R[2];
        let F2 = B.PQ;

        // The cost for the main loop is k doubles and 2*k differential adds.
        // If we normalise D1, D2, F1, F2 then we can save one mul per diff.
        // add, saving 2*k multiplications in total. As this function is usually
        // called with scalars of size log(p)/2 > 30, then it's worth normalising
        // the points.
        let mut inverses: [Fq; 4] = [D1.Z, D2.Z, F1.Z, F2.Z];
        Fq::batch_invert(&mut inverses);
        let mut xD1 = D1.X * inverses[0];
        let mut xD2 = D2.X * inverses[1];
        let mut xF1 = F1.X * inverses[2];
        let mut xF2 = F2.X * inverses[3];

        // Main ladder loop, compute [a]P + [b]Q
        for i in (0..k).rev() {
            // Compute indices.
            let r1 = r[2 * i] as usize;
            let r2 = r[2 * i + 1] as usize;
            let h = r1 + r2;

            // Compute new T values and swap differential points conditionally.
            let to_double = if (h >> 1) == 0 { R[h & 1] } else { R[2] };
            T[0] = self.xdbl(&to_double);

            Fq::condswap(&mut xD1, &mut xD2, (r2 as u32).wrapping_neg());
            T[1] = Self::xdiff_add_aff(&R[r2], &R[r2 + 1], &xD1);

            T[2] = Self::xdiff_add_aff(&R[0], &R[2], &xF1);
            Fq::condswap(&mut xF1, &mut xF2, ((h & 1) as u32).wrapping_neg());

            // Update R values from T values.
            R = T;
        }

        // When a and b are both even we want R[0], if a and b are both odd
        // we want R[2], otherwise we want R[1]
        let index = ((a[0] & 1) + (b[0] & 1)) as usize;
        R[index]
    }

    /// Variable time computation of [a]P + [b]Q given the x-only basis x(P), x(Q), x(P - Q)
    /// Based off Euclid2D: Algorithm 9 of https://eprint.iacr.org/2017/212.pdf
    pub fn ladder_biscalar_vartime(
        &self,
        B: &BasisX<Fq>,
        a: &[u8],
        b: &[u8],
        a_bitlen: usize,
        b_bitlen: usize,
    ) -> PointX<Fq> {
        // Convert from le bytes to le u64 big numbers
        let mut s0 = bn_from_le_bytes(a, a_bitlen);
        let mut s1 = bn_from_le_bytes(b, b_bitlen);

        // Ensure s0, s1 have the same length for the arithmetic logic
        let s_len = s0.len().max(s1.len());
        s0.resize(s_len, 0);
        s1.resize(s_len, 0);
        let mut bn_tmp: Vec<u64> = vec![0; s_len];

        // Define points we need for arithmetic
        let (mut x0, mut x1, mut xd) = (B.P, B.Q, B.PQ);

        // Euclidian loop, at the end s0 = 0, s1 = gcd(a, b) and x1 = [a/gcd(a, b) * P + b/gcd(a, b) * Q]
        while !bn_is_zero_vartime(&s0) {
            // Swap to ensure s1 >= s0
            if bn_lt_vartime(&s1, &s0) {
                (s0, s1) = (s1, s0);
                (x0, x1) = (x1, x0);

                // s1 >= s0. If the leading term of s1 is zero, remove it
                if s1[s1.len() - 1] == 0 {
                    s0.pop();
                    s1.pop();
                    bn_tmp.pop();
                }
            }

            // s1 <= 4 * s0
            bn_div4_vartime(&mut bn_tmp, &s1);
            if bn_lt_vartime(&bn_tmp, &s0) {
                // s0, s1 = s0, s1 - s0
                bn_sub_into_vartime(&mut s1, &s0);
                (x0, xd) = (Self::xdiff_add(&x1, &x0, &xd), x0);

            // s0 % 2 == s1 % 2:
            } else if s0[0] & 1 == s1[0] & 1 {
                // s0, s1 = s0, (s1 - s0) // 2
                bn_sub_into_vartime(&mut s1, &s0);
                bn_set_div2_vartime(&mut s1);
                (x1, x0) = self.xdbladd(&x1, &x0, &xd);

            // s1 % 2 == 0:
            } else if s1[0] & 1 == 0 {
                // s0, s1 = s0, s1 // 2
                bn_set_div2_vartime(&mut s1);
                (x1, xd) = self.xdbladd(&x1, &xd, &x0);
            } else {
                // s0, s1 = s0 // 2, s1
                bn_set_div2_vartime(&mut s0);
                (x0, xd) = self.xdbladd(&x0, &xd, &x1);
            }
        }

        // Clean trailing padded zeros
        while s1.len() > 1 && s1.last() == Some(&0) {
            s1.pop();
        }

        if s1.is_empty() {
            return PointX::INFINITY;
        }

        // Finalize by multiplying against the GCD eval
        if s1.len() == 1 {
            let mut n = s1[0];
            while n & 1 == 0 && n > 0 {
                n >>= 1;
                self.set_xdbl(&mut x1);
            }
            if n > 1 {
                x1 = self.xmul_u64_vartime(&x1, n);
            }
        } else {
            let mut gcd_bytes = Vec::with_capacity(s1.len() * 8);
            for w in s1 {
                gcd_bytes.extend_from_slice(&w.to_le_bytes());
            }
            x1 = self.xmul(&x1, &gcd_bytes, gcd_bytes.len() * 8);
        }

        x1
    }
}
