// TODO:
//
// Optimise Sqrt Velu, current implementation is very slow!
//
// Cofactor clearing. At the moment clearing the cofactor means multiplying by each ell_i e_i times,
// which seems silly. I think we should probably have some function which converts the prime factorisation into &[u64; ...]
// and then we use these limbs to do var time point multiplication.

// use std::time::Instant;

use fp2::traits::Fp as FqTrait;

use crate::{
    polynomial_ring::poly::Poly,
    utilities::bn::{
        bn_bit_length_vartime, factorisation_to_bn_vartime, prime_power_to_bn_vartime,
    },
};

use super::{curve::Curve, point::PointX};

const VELU_SQRT_THRESHOLD: usize = 200;

/// A structure which allows iterating over [i]P = (X : Z)
struct PointXMultiples<Fq: FqTrait> {
    P: PointX<Fq>,
    Q: PointX<Fq>,
    R: PointX<Fq>,
    i: usize,
}

impl<Fq: FqTrait> PointXMultiples<Fq> {
    pub fn new(A24: &Fq, C24: &Fq, P: &PointX<Fq>) -> Self {
        // precompute [2]P for the second output of multiplies
        let mut P2 = *P;
        Curve::xdbl_proj(A24, C24, &mut P2.X, &mut P2.Z);

        Self {
            P: *P,
            Q: *P,
            R: P2,
            i: 0,
        }
    }
}

impl<Fq: FqTrait> Iterator for PointXMultiples<Fq> {
    type Item = PointX<Fq>;

    fn next(&mut self) -> Option<Self::Item> {
        // Once R = [i]P = 0, we stop iterating as we have considered all non-zero
        // multiples.
        if self.R.is_zero() == u32::MAX {
            return None;
        }

        // We have to handle [1]P and [2]P differently than the rest
        // so for now I include a counter, would be nice to remove this
        // though...
        self.i += 1;

        if self.i == 1 {
            // For the first call, we just want to return [i]P = P
            return Some(self.P);
        } else if self.i == 2 {
            // For the second call, we want to return [2]P which has been computed
            // on creation of the type
            return Some(self.R);
        }
        // For all other calls, we want to set R = [i]P using differential addition
        let S = Curve::xdiff_add(&self.R, &self.P, &self.Q);
        (self.Q, self.R) = (self.R, S);

        Some(self.R)
    }
}

impl<Fq: FqTrait> Curve<Fq> {
    //============================================================
    // Variable time methods to compute [n]P and [n^ell]P for
    // cofactor clearing during isogeny computations. Should be
    // refactored to improve performance for large cofactors.

    /// P3 <- n*P, x-only variant using (A24 : C24).
    /// Integer n is represented as a u64 and is assumed to be public.
    pub fn set_xmul_proj_u64_vartime(A24: &Fq, C24: &Fq, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: u64) {
        // Handle small cases.
        match n {
            0 => {
                *P3 = PointX::INFINITY;
            }
            1 => {
                *P3 = *P;
            }
            2 => {
                *P3 = *P;
                Self::xdbl_proj(A24, C24, &mut P3.X, &mut P3.Z);
            }
            3 => {
                *P3 = *P;
                Self::xdbl_proj(A24, C24, &mut P3.X, &mut P3.Z);
                Self::xadd(&P.X, &P.Z, &P.X, &P.Z, &mut P3.X, &mut P3.Z);
            }
            4 => {
                *P3 = *P;
                Self::xdbl_proj(A24, C24, &mut P3.X, &mut P3.Z);
                Self::xdbl_proj(A24, C24, &mut P3.X, &mut P3.Z);
            }
            5 => {
                let mut P2 = *P;
                *P3 = *P;
                Self::xdbl_proj(A24, C24, &mut P2.X, &mut P2.Z);
                Self::xadd(&P.X, &P.Z, &P2.X, &P2.Z, &mut P3.X, &mut P3.Z);
                Self::xadd(&P.X, &P.Z, &P2.X, &P2.Z, &mut P3.X, &mut P3.Z);
            }

            _ => {
                let nbitlen = u64::BITS - n.leading_zeros();

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
                        let ctl = (((n >> i) as u32) & 1).wrapping_neg();
                        Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                        Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                        Self::xadd_aff(&Xp, &X0, &Z0, &mut X1, &mut Z1);
                        Self::xdbl_proj(A24, C24, &mut X0, &mut Z0);
                        cc = ctl;
                    }
                } else {
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n >> i) as u32) & 1).wrapping_neg();
                        Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                        Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                        Self::xadd(&P.X, &P.Z, &X0, &Z0, &mut X1, &mut Z1);
                        Self::xdbl_proj(A24, C24, &mut X0, &mut Z0);
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
                let spec = P.X.is_zero() & !P.Z.is_zero() & (((n & 1) as u32) & 1).wrapping_neg();
                P3.X = X0;
                P3.Z = Z0;
                P3.X.set_cond(&Fq::ZERO, spec);
                P3.Z.set_cond(&Fq::ONE, spec);
            }
        }
    }

    /// P3 <- n*P, x-only variant using (A24 : C24).
    /// Integer n is represented as a big integer with u64 words, little endian and is assumed to be public.
    fn set_xmul_proj_bn_vartime(
        A24: &Fq,
        C24: &Fq,
        P3: &mut PointX<Fq>,
        P: &PointX<Fq>,
        n: &[u64],
    ) {
        // When n has length n, call the lower level function.
        if n.len() == 1 {
            Self::set_xmul_proj_u64_vartime(A24, C24, P3, P, n[0]);
            return;
        }

        // Compute the bitlength of the big integer, at this point we know the
        // bit length will be more than 64 as we have a multi-limb bn.
        let nbitlen = bn_bit_length_vartime(n);

        let mut X0 = Fq::ONE;
        let mut Z0 = Fq::ZERO;
        let mut X1 = P.X;
        let mut Z1 = P.Z;
        let mut cc = 0u32;

        // As n is large enough, it is worthwhile to
        // normalize the source point to affine.
        // If P = inf, then this sets Xp to 0; thus, the
        // output of both xdbl() and xadd_aff() has Z = 0,
        // so we correctly get the point-at-infinity at the end.
        let Xp = P.X / P.Z;
        for i in (0..nbitlen).rev() {
            let ctl = (((n[i >> 6] >> (i & 63)) as u32) & 1).wrapping_neg();
            Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
            Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
            Self::xadd_aff(&Xp, &X0, &Z0, &mut X1, &mut Z1);
            Self::xdbl_proj(A24, C24, &mut X0, &mut Z0);
            cc = ctl;
        }
        Fq::condswap(&mut X0, &mut X1, cc);
        Fq::condswap(&mut Z0, &mut Z1, cc);

        // The ladder may fail if P = (0,0) (which is a point of
        // order 2) because in that case xadd() (and xadd_aff())
        // return Z = 0 systematically, so the result is considered
        // to be the point-at-infinity, which is wrong is n is odd.
        // We adjust the result in that case.
        let spec = P.X.is_zero() & !P.Z.is_zero() & (((n[0] & 1) as u32) & 1).wrapping_neg();
        P3.X = X0;
        P3.Z = Z0;
        P3.X.set_cond(&Fq::ZERO, spec);
        P3.Z.set_cond(&Fq::ONE, spec);
    }

    /// Return n*P as a new point (x-only variant) using (A24 : C24).
    /// Integer n is encoded as a u64 which is assumed to be a public value.
    pub fn xmul_proj_u64_vartime(A24: &Fq, C24: &Fq, P: &PointX<Fq>, n: u64) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        Self::set_xmul_proj_u64_vartime(A24, C24, &mut P3, P, n);
        P3
    }

    /// Return n*P as a new point (x-only variant) using (A24 : C24).
    /// Integer n is encoded as big number with little endian u64 words, which is assumed to be public.
    pub fn xmul_proj_bn_vartime(A24: &Fq, C24: &Fq, P: &PointX<Fq>, n: &[u64]) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        Self::set_xmul_proj_bn_vartime(A24, C24, &mut P3, P, n);
        P3
    }

    /// Return [n^e]*P as a new point (x-only variant) using (A24 : C24).
    /// Integer n is encoded as a u64 which is assumed to be a public value.
    fn xmul_proj_u64_iter_vartime(
        A24: &Fq,
        C24: &Fq,
        P: &PointX<Fq>,
        x: usize,
        e: usize,
    ) -> PointX<Fq> {
        // Convert x^e to a big integer represented as u64 words.
        let n = prime_power_to_bn_vartime(x, e);
        Self::xmul_proj_bn_vartime(A24, C24, P, &n)
    }

    // ============================================================
    // Internal functions for Velu isogenies with complexity O(ell)
    // suitable for all prime ell. Expects as input (A24 : C24) for
    // a domain instead of the Curve type to avoid inversions during
    // isogeny curve computations.

    /// Compute (X + Z) and (X - Z) for [i]P for 0 < i <= (ell - 1) / 2
    fn edwards_multiples(A24: &Fq, C24: &Fq, constants: &mut [(Fq, Fq)], P: &PointX<Fq>) {
        let mut iP = PointXMultiples::new(A24, C24, P);
        for c in constants.iter_mut() {
            let (X, Z) = iP.next().unwrap().coords();
            *c = (X - Z, X + Z)
        }
    }

    fn velu_two_isogeny_proj(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        img_points: &mut [PointX<Fq>],
    ) {
        debug_assert!(kernel.X.is_zero() != u32::MAX); // TODO: Handle this edge case

        let mut A_codomain = kernel.X.square();
        let C_codomain = kernel.Z.square();

        A_codomain.set_mul2();
        A_codomain = C_codomain - A_codomain;
        A_codomain.set_mul2();

        let t0 = kernel.X + kernel.Z;
        let t1 = kernel.X - kernel.Z;
        for P in img_points.iter_mut() {
            let mut t2 = P.X + P.Z;
            let mut t3 = P.Z - P.X;
            t3 *= t0;
            t2 *= t1;
            P.X *= t3 - t2;
            P.Z *= t3 + t2;
        }

        let mut C24_cod = C_codomain.mul2();
        let A24_cod = A_codomain + C24_cod;
        C24_cod.set_mul2();

        *A24 = A24_cod;
        *C24 = C24_cod;
    }

    pub fn velu_odd_isogeny_proj(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        degree: usize,
        img_points: &mut [PointX<Fq>],
    ) {
        // Convert from Montgomery to projective twisted Edwards (A_ed : D_ed)
        let mut A_ed = *A24; // A_ed = (A + 2*C)
        let mut D_ed = *A24 - *C24; // D_ed = (A - 2*C)

        // We precompute (X - Z) and (X + Z) for (X : Z) = [i]P for i in 0..((ell - 1)/2)
        let d = (degree - 1) >> 1;
        let mut constants = vec![(Fq::ZERO, Fq::ZERO); d];
        Self::edwards_multiples(A24, C24, &mut constants, kernel);

        // Compute the product of the edward multiples
        let mut prod_Y = Fq::ONE;
        let mut prod_Z = Fq::ONE;
        for (Y_ed, Z_ed) in constants.iter() {
            prod_Y *= *Y_ed;
            prod_Z *= *Z_ed;
        }

        // Compute prod_Y^8 and prod_Z^8
        for _ in 0..3 {
            prod_Y.set_square();
            prod_Z.set_square();
        }

        // Compute the new codomain in projective twisted Edwards
        // A_new = A_old^ell * prod_Z^8
        // D_new = D_old^ell * prod_Y^8
        A_ed.set_pow_u64_vartime(degree as u64);
        D_ed.set_pow_u64_vartime(degree as u64);
        A_ed *= prod_Z;
        D_ed *= prod_Y;

        // Evaluate each point through the isogeny
        for P in img_points.iter_mut() {
            let P_sum = P.X + P.Z;
            let P_diff = P.X - P.Z;

            let mut EY_sum;
            let mut EZ_diff;
            let mut X_new = Fq::ONE;
            let mut Z_new = Fq::ONE;
            for (Y_ed, Z_ed) in constants.iter() {
                EZ_diff = *Z_ed * P_diff;
                EY_sum = *Y_ed * P_sum;
                X_new *= EZ_diff + EY_sum;
                Z_new *= EZ_diff - EY_sum;
            }

            P.X *= X_new.square();
            P.Z *= Z_new.square();
        }

        // Convert back to Montgomery (A24 : C24)
        *A24 = A_ed;
        *C24 = A_ed - D_ed;
    }

    // ============================================================
    // Sqrt Velu functions for large ell-isogenies
    // WARNING: this is underperforming currently due to missing optimisations in some
    // specialised polynomial arithmetic to allow for fast scaled remainder trees and
    // other smaller things.

    /// Precompute the points in the three partitions I, J and K using x-only arithmetic.
    // TODO: study Algorithm 1 of ia.cr/2024/584 and see if this still will lead to
    // performance gains for our implementation.
    fn precompute_partitions(
        xI: &mut [PointX<Fq>],
        xJ: &mut [PointX<Fq>],
        xK: &mut [PointX<Fq>],
        A24: &Fq,
        C24: &Fq,
        P: &PointX<Fq>,
    ) {
        let size_I = xI.len();
        let size_J = xJ.len();
        let size_K = xK.len();

        debug_assert!(size_I >= size_J);
        debug_assert!(size_J > 1);

        let mut P2 = *P;
        Self::xdbl_proj(A24, C24, &mut P2.X, &mut P2.Z);

        // First we compute [j]P for j in {1, 3, ... 2*size_J - 1}
        xJ[0] = *P;
        xJ[1] = Self::xdiff_add(&xJ[0], &P2, &xJ[0]);
        for i in 2..size_J {
            xJ[i] = Self::xdiff_add(&xJ[i - 1], &P2, &xJ[i - 2]);
        }

        // Next we want to compute x([i]P) for i in {2*size_J * (2i + 1)}
        let mut P4 = P2;
        Self::xdbl_proj(A24, C24, &mut P4.X, &mut P4.Z);

        // First we compute [2*size_J]P which we can do efficiently from the computation
        // of xJ above. We assume the degree of the isogeny is known (and so size_J is also
        // known, so we allow the branch (we could do two conditional swaps).
        let b_half_floor = size_J / 2;
        let b_half_ceil = size_J - b_half_floor;
        let Q = if (size_J % 2) == 1 {
            Self::xdiff_add(&xJ[b_half_ceil], &xJ[b_half_floor - 1], &P4)
        } else {
            Self::xdiff_add(&xJ[b_half_ceil], &xJ[b_half_floor - 1], &P2)
        };

        // We need [2]Q as a step-size to generate xI
        let mut Q2 = Q;
        Self::xdbl_proj(A24, C24, &mut Q2.X, &mut Q2.Z);

        xI[0] = Q;
        xI[1] = Self::xdiff_add(&xI[0], &Q2, &xI[0]);
        for i in 2..size_I {
            xI[i] = Self::xdiff_add(&xI[i - 1], &Q2, &xI[i - 2]);
        }

        // Finally we compute [k]P for k in {4*sJ*sI + 1, ...,  ell - 4, ell - 2}
        // We need to be careful when xK is small, so there's a few early returns
        // below.
        if xK.is_empty() {
            return;
        }
        xK[0] = P2;
        if size_K == 1 {
            return;
        }
        xK[1] = P4;
        for i in 2..size_K {
            xK[i] = Self::xdiff_add(&xK[i - 1], &P2, &xK[i - 2]);
        }
    }

    #[inline]
    fn precompute_eJ_values(
        eJ_coeffs: &mut [(Fq, Fq, Fq)],
        hJ_points: &[PointX<Fq>],
        A24: &Fq,
        C24: &Fq,
    ) {
        // TODO: we could work projectively here with (A : C) instead of A. This saves us one inversion
        // (about 60M for CSIDH 512 or 30M for 256-bit Fp2) but adds multiplications by C throughout codomain
        // and image computations.
        //
        // I believe this ends up being 2M * J extra for codomain computations and 2M * J extra for each
        // image. We expect velu to start being better for ell = 100, so we're paying a cost of 60M - 20M(1 + n)
        // at this degree where n is the number of images we are computing. For ell = 587, the max for
        // CSIDH 512 we're paying a cost of 60M - 48M(1 + n) and so is a saving when we do one eval.
        let A = (*A24 / *C24).mul4() - Fq::TWO;
        for (i, P) in hJ_points.iter().enumerate() {
            let (X, Z) = P.coords();
            let XZ = X * Z;

            let add_sqr = (X + Z).square();
            let XZ4neg = -XZ.mul4();
            let AXZ4neg = A * XZ4neg;

            eJ_coeffs[i] = (add_sqr, XZ4neg, AXZ4neg);
        }
    }

    /// Compute the product of an array of Fq values using a product tree.
    // TODO: this could be placed into the Fp2 library.
    fn product_tree_root_fq(v: &[Fq]) -> Fq {
        if v.is_empty() {
            return Fq::ONE;
        }
        if v.len() == 1 {
            return v[0];
        }
        let half = v.len() >> 1;
        Self::product_tree_root_fq(&v[..half]) * Self::product_tree_root_fq(&v[half..])
    }

    /// Understanding the polynomial hK = prod(x * PZ - PX) for the set in hK, then
    /// evaluate this polynomial at alpha = 1 and alpha = -1
    #[inline]
    fn hK_codomain(hK_points: &[PointX<Fq>]) -> (Fq, Fq) {
        // We have the factorisation of hK into linear pieces, so we evaluate
        // each factor to get an Fq element and then compute the product of these
        // with a product tree.
        let mut h1_linear = Vec::with_capacity(hK_points.len());
        let mut h2_linear = Vec::with_capacity(hK_points.len());

        for P in hK_points.iter() {
            h1_linear.push(P.Z - P.X);
            h2_linear.push(-(P.Z + P.X));
        }

        let h1 = Self::product_tree_root_fq(&h1_linear);
        let h2 = Self::product_tree_root_fq(&h2_linear);

        (h1, h2)
    }

    /// Understanding the polynomial hK = prod(x * PZ - PX) for the set in hK, then
    /// evaluate this polynomial at alpha and 1/alpha (projectively) where we have
    /// alpha = (X : Z) and we take as input X + Z and X - Z.
    #[inline]
    fn hK_eval(hK_points: &[PointX<Fq>], XpZ: &Fq, XmZ: &Fq) -> (Fq, Fq) {
        // We have the factorisation of hK into linear pieces, so we evaluate
        // each factor to get an Fq element and then compute the product of these
        // with a product tree.
        let mut h1_linear = Vec::with_capacity(hK_points.len());
        let mut h2_linear = Vec::with_capacity(hK_points.len());

        let mut t1;
        let mut t2;
        for P in hK_points.iter() {
            t1 = P.X + P.Z;
            t1 *= *XmZ;

            t2 = P.X - P.Z;
            t2 *= *XpZ;

            h1_linear.push(t1 - t2);
            h2_linear.push(t1 + t2);
        }

        let h1 = Self::product_tree_root_fq(&h1_linear);
        let h2 = Self::product_tree_root_fq(&h2_linear);

        (h1, h2)
    }

    /// Compute an isogeny using the sqrt-velu algorithm O(\sqrt{ell}) complexity.
    /// https://velusqrt.isogeny.org
    pub fn sqrt_velu_odd_isogeny_proj<P: Poly<Fq>>(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        degree: usize,
        img_points: &mut [PointX<Fq>],
    ) {
        // baby step, giant step values
        // TODO: better sqrt?
        // TODO: I need to use degree - 1 rather than degree + 1 as in
        // the paper to avoid K < 0... Try and get to the bottom of this.
        let size_J = (((degree - 1) as f64).sqrt() as usize) / 2;
        let size_I = (degree - 1) / (4 * size_J);
        let size_K = (degree - 4 * size_J * size_I - 1) / 2;

        // Compute the points in the I, J and K partitions.
        let mut hI_points = vec![PointX::INFINITY; size_I];
        let mut hJ_points = vec![PointX::INFINITY; size_J];
        let mut hK_points = vec![PointX::INFINITY; size_K];
        Self::precompute_partitions(
            &mut hI_points,
            &mut hJ_points,
            &mut hK_points,
            A24,
            C24,
            kernel,
        );

        // Precompute the roots of the polynomial Prod(x - [i]P.x()) for i in the set I
        // In another implementation, we might instead compute the polynomial, which we
        // would use as input into a remainder tree instead...
        PointX::batch_normalise(&mut hI_points);
        let hI_roots: Vec<Fq> = hI_points.iter().map(|P| P.X).collect();

        // Precompute (X + Z)^2, (X - Z)^2, -4XZ and -4XZ*A
        let mut eJ_precomp = vec![(Fq::ZERO, Fq::ZERO, Fq::ZERO); size_J];
        Self::precompute_eJ_values(&mut eJ_precomp, &hJ_points, A24, C24); // Cost: size_J * (1S + 2M)

        let mut E0J_leaves: Vec<P> = Vec::with_capacity(size_J);
        let mut E1J_leaves: Vec<P> = Vec::with_capacity(size_J);
        for (sum_sqr, XZ4neg, AXZ4neg) in eJ_precomp.iter() {
            // (X - Z)^2 = (X + Z)^2 - 4 * X * Z
            let c0_0 = *sum_sqr + *XZ4neg;
            let c0_1 = *AXZ4neg - sum_sqr.mul2();

            let c1_0 = *sum_sqr;
            let c1_1 = c0_0.mul2() - *AXZ4neg;

            // Each quadratic factor here is a palindrome.
            // TODO: we could write specialised polynomial arithmetic which
            // abuses this for faaster multiplication.
            E0J_leaves.push(P::new_from_slice(&[c0_0, c0_1, c0_0]));
            E1J_leaves.push(P::new_from_slice(&[c1_0, c1_1, c1_0]));
        }
        let E0J = P::product_tree_root(&E0J_leaves);
        let E1J = P::product_tree_root(&E1J_leaves);
        debug_assert!(E0J.degree().unwrap() == 2 * size_J);
        debug_assert!(E1J.degree().unwrap() == 2 * size_J);

        // Compute the codomain.
        let r0 = E0J.resultant_from_roots(&hI_roots);
        let r1 = E1J.resultant_from_roots(&hI_roots);
        let (m0, m1) = Self::hK_codomain(&hK_points);

        // Compute (ri * mi)^8 * (A ∓ 2)^degree
        let mut num = r0 * m0;
        let mut den = r1 * m1;
        for _ in 0..3 {
            num.set_square();
            den.set_square();
        }

        let mut A_ed = *A24; // A_ed = (A + 2*C)
        let mut D_ed = *A24 - *C24; // D_ed = (A - 2*C)
        A_ed.set_pow_u64_vartime(degree as u64);
        D_ed.set_pow_u64_vartime(degree as u64);

        A_ed *= den;
        D_ed *= num;

        // Evaluate each point through the isogeny.
        for P in img_points.iter_mut() {
            if P.is_zero() == u32::MAX {
                continue;
            }

            // We use the sum and difference for the EJ leaves as well as when
            // evaluating hK at alpha = (X / Z) and 1/alpha.
            let XpZ = P.X + P.Z;
            let XmZ = P.X - P.Z;
            let XZ2 = (P.X * P.Z).mul2();
            let X2Z2 = XpZ.square() - XZ2; // X^2 + Z^2

            let mut E0J_leaves: Vec<P> = Vec::with_capacity(size_J);
            for (i, (sum_sqr, XZ4neg, AXZ4neg)) in eJ_precomp.iter().enumerate() {
                let Pj = hJ_points[i];

                // Precompute some multiplications for later.
                let add = Pj.X + Pj.Z;
                let sub = Pj.X - Pj.Z;

                // Constant coefficient: c0 = [2 * (X * Xj - Z * Zj)]^2
                // Quadratic coefficient: c1 = [2 * (X * Zj - Z * Xj)]^2
                let t1 = XmZ * add;
                let t2 = XpZ * sub;
                let c0 = (t1 + t2).square();
                let c2 = (t1 - t2).square();

                // Linear coefficent is given by three terms
                // - [2 * (Xj^2 + Zj^2)] * 2 X Z +
                //   (X^2 + Z^2) * (-4 * Xj * Zj) +
                //   (2 A X Z) * (-4 Xj * Zj).

                // Note that we use: 2 * (Xj + Zj)^2 - 4 Xj Zj = 2 * (Xj^2 + Zj^2)
                let mut c1 = -sum_sqr.mul2();
                c1 -= *XZ4neg;
                c1 += *AXZ4neg;
                c1 *= XZ2;
                c1 += X2Z2 * *XZ4neg;
                c1.set_mul2();

                E0J_leaves.push(P::new_from_slice(&[c0, c1, c2]));
            }

            let E0J = P::product_tree_root(&E0J_leaves);
            let E1J = E0J.reverse();

            let r0 = E0J.resultant_from_roots(&hI_roots);
            let r1 = E1J.resultant_from_roots(&hI_roots);
            let (m0, m1) = Self::hK_eval(&hK_points, &XpZ, &XmZ);

            let r0m0 = r0 * m0;
            let r1m1 = r1 * m1;

            P.X *= r1m1.square();
            P.Z *= r0m0.square();
        }

        // Convert back to Montgomery (A24 : C24)
        *A24 = A_ed;
        *C24 = A_ed - D_ed;
    }

    // ============================================================
    // Internal methods for computing isogeny chains of degree ell^e
    // and generic composite orders of degree \prod ell_i^ei

    fn velu_prime_power_isogeny_proj<P: Poly<Fq>>(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        degree: usize,
        len: usize,
        img_points: &mut [PointX<Fq>],
    ) {
        // We push the image points through at the same time as the strategy
        // points, so we need to know how many images we're computing to keep
        // track of them for the end.
        let n = img_points.len();

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - len.leading_zeros() + 1) as usize;

        // These are a set of points of order ell^i, orders keeps track of the order i
        // we also prepend the points we want to evaluate onto this vector.
        let mut stategy_points: Vec<PointX<Fq>> = vec![PointX::INFINITY; space + n];
        let mut orders: Vec<usize> = vec![0; space];

        // Set the first elements of the vector to the points we want to push
        // through the isogeny.
        stategy_points[..n].copy_from_slice(img_points);

        // Then set the next element to be the kernel input
        stategy_points[n] = *kernel;
        orders[0] = len;

        // Compute the isogeny chain with a naive balanced strategy.
        let mut k = 0;
        for _ in 0..len {
            // Get the next point of order 2
            while orders[k] != 1 {
                k += 1;
                let m = orders[k - 1] / 2;
                stategy_points[n + k] = stategy_points[n + k - 1];

                // when ell = 2 we can do repeated doubling
                if degree == 2 {
                    Self::xdbl_proj_iter(A24, C24, &mut stategy_points[n + k], m);
                } else {
                    // Otherwise we have this janky repeated multiplication.
                    stategy_points[n + k] = Self::xmul_proj_u64_iter_vartime(
                        A24,
                        C24,
                        &stategy_points[n + k],
                        degree,
                        m,
                    )
                }
                orders[k] = orders[k - 1] - m;
            }

            // Point of order ell to compute isogeny
            let ker_step = stategy_points[n + k];

            // Compute the ell-isogeny
            if degree == 2 {
                Self::velu_two_isogeny_proj(A24, C24, &ker_step, &mut stategy_points[..(k + n)]);
            } else if degree < VELU_SQRT_THRESHOLD {
                Self::velu_odd_isogeny_proj(
                    A24,
                    C24,
                    &ker_step,
                    degree,
                    &mut stategy_points[..(k + n)],
                );
            } else {
                Self::sqrt_velu_odd_isogeny_proj::<P>(
                    A24,
                    C24,
                    &ker_step,
                    degree,
                    &mut stategy_points[..(k + n)],
                );
            }

            // Reduce the order of the points we evaluated and decrease k
            for ord in orders.iter_mut().take(k) {
                *ord -= 1;
            }
            k = k.saturating_sub(1);
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&stategy_points[..n]);
    }

    fn velu_composite_isogeny_proj<P: Poly<Fq>>(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        degrees: &[(usize, usize)],
        img_points: &mut [PointX<Fq>],
    ) {
        let mut eval_points = img_points.to_vec();
        eval_points.push(*kernel);

        for i in 0..degrees.len() {
            // At each step we will compute a degree^len isogeny
            let (ell, n) = degrees[i];

            // Clear the cofactor
            // TODO: this will be slow, must be a better way to batch things into the u64
            // ker_step will be a point with order ell^n
            let mut ker_step = *eval_points.last().unwrap();

            // Collect the factorisation of the cofactor by skipping the ell already handled.
            // We don't have to do anything on the last factor of the isogeny.
            let factorisation: Vec<(usize, usize)> = degrees.iter().skip(i + 1).cloned().collect();
            if !factorisation.is_empty() {
                let cofactor: Vec<u64> = factorisation_to_bn_vartime(&factorisation);
                ker_step = Self::xmul_proj_bn_vartime(A24, C24, &ker_step, &cofactor);
            }

            Self::velu_prime_power_isogeny_proj::<P>(A24, C24, &ker_step, ell, n, &mut eval_points)
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);
    }

    // ============================================================
    // Public functions which compute isogenies given user-friendly
    // inputs and types etc.

    pub fn velu_prime_isogeny<P: Poly<Fq>>(
        self,
        kernel: &PointX<Fq>,
        degree: usize,
        img_points: &mut [PointX<Fq>],
    ) -> Self {
        let mut A24 = self.A + Fq::TWO;
        let mut C24 = Fq::FOUR;

        // 2-isogenies are handled with a special function
        if degree == 2 {
            Self::velu_two_isogeny_proj(&mut A24, &mut C24, kernel, img_points);
        } else if degree < VELU_SQRT_THRESHOLD {
            Self::velu_odd_isogeny_proj(&mut A24, &mut C24, kernel, degree, img_points);
        } else {
            Self::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, kernel, degree, img_points);
        }
        Self::curve_from_A24_proj(&A24, &C24)
    }

    pub fn velu_prime_power_isogeny<P: Poly<Fq>>(
        self,
        kernel: &PointX<Fq>,
        degree: usize,
        len: usize,
        img_points: &mut [PointX<Fq>],
    ) -> Self {
        let mut A24 = self.A + Fq::TWO;
        let mut C24 = Fq::FOUR;
        Self::velu_prime_power_isogeny_proj::<P>(
            &mut A24, &mut C24, kernel, degree, len, img_points,
        );

        Self::curve_from_A24_proj(&A24, &C24)
    }

    pub fn velu_composite_isogeny<P: Poly<Fq>>(
        self,
        kernel: &PointX<Fq>,
        degrees: &[(usize, usize)],
        img_points: &mut [PointX<Fq>],
    ) -> Self {
        let mut A24 = self.A + Fq::TWO;
        let mut C24 = Fq::FOUR;

        Self::velu_composite_isogeny_proj::<P>(&mut A24, &mut C24, kernel, degrees, img_points);

        Self::curve_from_A24_proj(&A24, &C24)
    }
}
