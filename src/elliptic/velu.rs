// TODO:
//
// Cofactor clearing. At the moment clearing the cofactor means multiplying by each ell_i e_i times,
// which seems silly. I think we should probably have some function which converts the prime factorisation into &[u64; ...]
// and then we use these limbs to do var time point multiplication.

use fp2::traits::Fp as FqTrait;

use crate::polynomial_ring::poly::Poly;

use super::{curve::Curve, point::PointX};

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
    fn set_xmul_proj_u64_vartime(A24: &Fq, C24: &Fq, P3: &mut PointX<Fq>, P: &PointX<Fq>, n: u64) {
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
                let nbitlen = 63 - n.leading_ones();

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

    /// Return n*P as a new point (x-only variant) using (A24 : C24).
    /// Integer n is encoded as a u64 which is assumed to be a public value.
    fn xmul_proj_u64_vartime(A24: &Fq, C24: &Fq, P: &PointX<Fq>, n: u64) -> PointX<Fq> {
        let mut P3 = PointX::INFINITY;
        Self::set_xmul_proj_u64_vartime(A24, C24, &mut P3, P, n);
        P3
    }

    /// Return [n^e]*P as a new point (x-only variant) using (A24 : C24).
    /// Integer n is encoded as a u64 which is assumed to be a public value.
    fn xmul_proj_u64_iter_vartime(
        A24: &Fq,
        C24: &Fq,
        P: &PointX<Fq>,
        n: u64,
        e: usize,
    ) -> PointX<Fq> {
        // If n^e fits inside a u64 then we simply send this
        let (x, overflow) = n.overflowing_pow(e as u32);
        if !overflow {
            return Self::xmul_proj_u64_vartime(A24, C24, P, x);
        }

        // Compute the largest power y such that n^y fits inside a u64
        // and represent n^e = n^(k * y + r)
        let y = (64 / (64 - n.leading_zeros())) as usize;
        let k = e / y;
        let r = e % y;
        debug_assert!(y * k + r == e);

        // Compute n^y and n^r for the multiplications below.
        // TODO: this is still wasteful, we should represent
        // the scalar as a [u64] and pass this to a single
        // function!
        let n_y = n.wrapping_pow(y as u32);
        let n_r = n.wrapping_pow(r as u32);

        let mut P3 = *P;
        for _ in 0..k {
            P3 = Self::xmul_proj_u64_vartime(A24, C24, &P3, n_y);
        }
        P3 = Self::xmul_proj_u64_vartime(A24, C24, &P3, n_r);

        P3
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
    // WARNING: this is underperforming because of inefficienies absolutely everywhere!

    /// Compute roots of the polynomial \Prod (Z - x(Q)) for Q in the set
    /// I = {2b(2i + 1) | 0 <= i < c}
    fn precompute_hi_roots(hI_roots: &mut [Fq], A24: &Fq, C24: &Fq, kernel: &PointX<Fq>, b: usize) {
        // Set Q = [2b]K
        let mut Q = Self::xmul_proj_u64_vartime(A24, C24, kernel, (b + b) as u64);

        // Initalise values for the addition chain, we repeatedly add [2]Q for each step
        let mut step = Q;
        Self::xdbl_proj(A24, C24, &mut step.X, &mut step.Z);
        let mut diff = Q;

        // TODO: work projectively here to avoid the inversion for each point!
        for r in hI_roots.iter_mut() {
            *r = Q.x();
            // TODO: we can skip the last addition
            let R = Self::xdiff_add(&Q, &step, &diff);
            diff = Q;
            Q = R;
        }
    }

    /// Given a point (X : Z) compute the three elliptic resultants on the Montgomery
    /// curve, given by (with `a = X / Z` for now...).
    /// -
    /// - `f1 = (x * QZ - QX)^2
    /// - `f2 = -2((xa + 1)(x + a) + 2Axa)`
    /// - `f3 = -2* (2*A*QX*QZ*x + (QX*x + QZ)*(QZ*x + QX))
    fn _old_elliptic_resultants<P: Poly<Fq>>(A: &Fq, Q: &PointX<Fq>) -> (P, P, P) {
        // TODO: is it better to invert A / C to be used here, or work with (A : C)
        // projectively? The cost balance is 1 inversion (~30M) or the cost of
        // ~sqrt(degree) * n multiplications. Where n is the number of multiplications
        // by C in this function.

        // TODO: do operation counting to see if this is the fastest way
        // to compute these three polynomials...
        let QX_sqr = Q.X.square();
        let QZ_sqr = Q.Z.square();
        let QXQZ = Q.X * Q.Z;
        let QXQZ_m2 = -QXQZ.mul2();

        // f1 = (x * QZ - QX)^2
        let f1 = P::new_from_slice(&[QX_sqr, QXQZ_m2, QZ_sqr]);

        // f2 = -2*A*QX*QZ*x - 2*(QX*x + QZ)*(QZ*x + QX)
        //    = -2 * (QX*QZ*x^2 + (A*QX*QZ + QX^2 + QZ^2)*x + QX*QZ)
        let f2_linear = -(*A * QXQZ + QX_sqr + QZ_sqr).mul2();
        let f2 = P::new_from_slice(&[QXQZ_m2, f2_linear, QXQZ_m2]);

        // f3 = (x * QX - QZ)^2
        let f3 = P::new_from_slice(&[QZ_sqr, QXQZ_m2, QX_sqr]);

        (f1, f2, f3)
    }

    /// Given a point (X : Z) we want to compute three elliptic resultants on the Montgomery
    /// curve, but of the three quadratic polynomials, there are only actually four unique
    /// coefficients. QX^2, QZ^2, -2*QX*QZ and -2 * (A*QX*QZ + QX^2 + QZ^2). This is what
    /// we precompute and store for later computations.
    fn elliptic_resultants(A: &Fq, Q: &PointX<Fq>) -> (Fq, Fq, Fq, Fq) {
        // TODO: is it better to invert A / C to be used here, or work with (A : C)
        // projectively? The cost balance is 1 inversion (~30M) or the cost of
        // ~sqrt(degree) * n multiplications. Where n is the number of multiplications
        // by C in this function.

        // TODO: do operation counting to see if this is the fastest way
        // to compute these three polynomials...
        let QX_sqr = Q.X.square();
        let QZ_sqr = Q.Z.square();
        let QXQZ = Q.X * Q.Z;
        let QXQZ_m2 = -QXQZ.mul2();

        let f2_linear = -(*A * QXQZ + QX_sqr + QZ_sqr).mul2();

        (QX_sqr, QZ_sqr, QXQZ_m2, f2_linear)
    }

    fn precompute_eJ_coeffs(
        eJ_coeffs: &mut [(Fq, Fq, Fq, Fq)],
        A24: &Fq,
        C24: &Fq,
        ker: &PointX<Fq>,
    ) {
        // TODO: this is dumb, but we need A rather than (A24 : C24) for the
        // current formula...
        let A = (*A24 / *C24).mul4() - Fq::TWO;

        // Initalise values for the addition chain, we repeatedly add [2]Q for each step
        // to compute [j]Q for j in {1, 3, 5, 7, ...}
        let mut Q = *ker;
        let mut step = Q;
        Self::xdbl_proj(A24, C24, &mut step.X, &mut step.Z);
        let mut diff = Q;

        for eJ_coeffs in eJ_coeffs.iter_mut() {
            *eJ_coeffs = Self::elliptic_resultants(&A, &Q);
            // TODO: we can skip the last addition
            let R = Self::xdiff_add(&Q, &step, &diff);
            diff = Q;
            Q = R;
        }
    }

    fn precompute_hK_points(hK_points: &mut [PointX<Fq>], A24: &Fq, C24: &Fq, kernel: &PointX<Fq>) {
        let mut Q = *kernel;
        Self::xdbl_proj(A24, C24, &mut Q.X, &mut Q.Z);
        let step = Q;
        let mut R = Q;
        Self::xdbl_proj(A24, C24, &mut R.X, &mut R.Z);

        for P in hK_points.iter_mut() {
            *P = Q;

            // TODO: we can skip this last add in the loop
            let S = Self::xdiff_add(&R, &step, &Q);
            (Q, R) = (R, S)
        }
    }

    fn hK_codomain(hK_points: &[PointX<Fq>]) -> (Fq, Fq) {
        let mut h1 = Fq::ONE;
        let mut h2 = Fq::ONE;
        for P in hK_points.iter() {
            h1 *= P.Z - P.X;
            h2 *= -(P.Z + P.X);
        }
        (h1, h2)
    }

    fn hK_eval(hK_points: &[PointX<Fq>], alpha: &Fq) -> (Fq, Fq) {
        let mut h1 = Fq::ONE;
        let mut h2 = Fq::ONE;
        for P in hK_points.iter() {
            h1 *= P.Z - P.X * *alpha;
            h2 *= *alpha * P.Z - P.X;
        }
        (h1, h2)
    }

    pub fn sqrt_velu_odd_isogeny_proj<P: Poly<Fq>>(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        degree: usize,
        img_points: &mut [PointX<Fq>],
    ) {
        // baby step, giant step values
        // TODO: better sqrt?
        let b = (((degree + 1) as f64).sqrt() as usize) / 2;
        let c = (degree + 1) / (4 * b);
        let stop = degree - 4 * b * c;

        // Precompute the roots of the polynomial Prod(x - [i]P.x()) for i in the set I
        let mut hI_roots = vec![Fq::ZERO; c];
        Self::precompute_hi_roots(&mut hI_roots, A24, C24, kernel, b);

        // Precompute three polynomials for each point [j]P for j in J
        let mut eJ_coeffs = vec![(Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ZERO); b];
        Self::precompute_eJ_coeffs(&mut eJ_coeffs, A24, C24, kernel);

        // Precompute the polynomial (x - [k]P.x()) for k in the set K
        let mut hK_points = vec![PointX::INFINITY; (stop - 1) / 2];
        Self::precompute_hK_points(&mut hK_points, A24, C24, kernel);

        // Compute prod(f0 + f1 + f2) and prod(f0 - f1 + f2) which we will compute
        // resultants of with hI_roots.
        let mut E0J = P::new_from_ele(&Fq::ONE);
        let mut E1J = P::new_from_ele(&Fq::ONE);
        for (QX_sqr, QZ_sqr, f1_1, f2_1) in eJ_coeffs.iter() {
            let c0_0 = *QX_sqr + *f1_1 + *QZ_sqr;
            let c0_1 = f1_1.mul2() + *f2_1;
            let c1_0 = c0_0 - f1_1.mul2();
            let c1_1 = c0_1 - f2_1.mul2();

            // t0 = f0 + f1 + f2 using that f0 = f2.reverse()
            // t1 = f0 - f1 + f2 using that f0 = f2.reverse()
            let t0 = P::new_from_slice(&[c0_0, c0_1, c0_0]);
            let t1 = P::new_from_slice(&[c1_0, c1_1, c1_0]);

            E0J *= t0;
            E1J *= t1;
        }

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
        //
        // Here we have for each point P (X : Z) the option to compute
        // the fi elliptic polynomials of degree 2 with scaling by (X/Z) twice per step
        // at a cost of b * 6M or by computing PX^2, PZ^2 and PX * PZ and then scaling
        // three times at a cost of b * (9M  + 1M + 2S).
        //
        // Additionally in hK_eval if we use alpha we need `stop` additional
        // multiplications.
        //
        // The question is for what values of b and stop should we compute alpha = X/Z
        // to save multiplications for the cost of about 30M for the inversion for alpha
        //
        // Roughly it seems that if b = 5 it seems to be worth inverting and seeing as we
        // expect to use this for degree ~ 100 I think the inversion is always good?
        for P in img_points.iter_mut() {
            if P.is_zero() == u32::MAX {
                continue;
            }

            let alpha = P.x();
            let mut E1J = P::new_from_ele(&Fq::ONE);
            for (QX_sqr, QZ_sqr, f1_1, f2_1) in eJ_coeffs.iter() {
                // TODO: can I do better here than making the polynomials and then
                // scaling? Probably...
                let mut f0 = P::new_from_slice(&[*QX_sqr, *f1_1, *QZ_sqr]);
                let f1 = P::new_from_slice(&[*f1_1, *f2_1, *f1_1]);
                let f2 = P::new_from_slice(&[*QZ_sqr, *f1_1, *QX_sqr]);

                // Compute alpha^2 * f0 + alpha f1 + f2
                f0 *= alpha;
                f0 += f1.clone();
                f0 *= alpha;
                f0 += f2.clone();

                E1J *= f0;
            }
            let E0J = E1J.reverse();

            let r0 = E0J.resultant_from_roots(&hI_roots);
            let r1 = E1J.resultant_from_roots(&hI_roots);
            let (m0, m1) = Self::hK_eval(&hK_points, &alpha);

            let r0m0 = r0 * m0;
            let r1m1 = r1 * m1;

            P.X = r0m0.square() * alpha;
            P.Z = r1m1.square();
        }

        // Convert back to Montgomery (A24 : C24)
        *A24 = A_ed;
        *C24 = A_ed - D_ed;
    }

    // ============================================================
    // Internal methods for computing isogeny chains of degree ell^e
    // and generic composite orders of degree \prod ell_i^ei

    fn velu_prime_power_isogeny_proj(
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
                        degree as u64,
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
            } else {
                Self::velu_odd_isogeny_proj(
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

    fn velu_composite_isogeny_proj(
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
            for (p, e) in degrees.iter().skip(i + 1) {
                ker_step = Self::xmul_proj_u64_iter_vartime(A24, C24, &ker_step, *p as u64, *e);
            }
            Self::velu_prime_power_isogeny_proj(A24, C24, &ker_step, ell, n, &mut eval_points)
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);
    }

    // ============================================================
    // Public functions which compute isogenies given user-friendly
    // inputs and types etc.

    pub fn velu_prime_isogeny(
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
        } else {
            Self::velu_odd_isogeny_proj(&mut A24, &mut C24, kernel, degree, img_points);
        }
        Self::curve_from_A24_proj(&A24, &C24)
    }

    pub fn velu_prime_power_isogeny(
        self,
        kernel: &PointX<Fq>,
        degree: usize,
        len: usize,
        img_points: &mut [PointX<Fq>],
    ) -> Self {
        let mut A24 = self.A + Fq::TWO;
        let mut C24 = Fq::FOUR;
        Self::velu_prime_power_isogeny_proj(&mut A24, &mut C24, kernel, degree, len, img_points);

        Self::curve_from_A24_proj(&A24, &C24)
    }

    pub fn velu_composite_isogeny(
        self,
        kernel: &PointX<Fq>,
        degrees: &[(usize, usize)],
        img_points: &mut [PointX<Fq>],
    ) -> Self {
        let mut A24 = self.A + Fq::TWO;
        let mut C24 = Fq::FOUR;

        Self::velu_composite_isogeny_proj(&mut A24, &mut C24, kernel, degrees, img_points);

        Self::curve_from_A24_proj(&A24, &C24)
    }
}
