// TODO:
//
// Cofactor clearing. At the moment clearing the cofactor means multiplying by each ell_i e_i times,
// which seems silly. I think we should probably have some function which converts the prime factorisation into &[u64; ...]
// and then we use these limbs to do var time point multiplication.

use fp2::traits::Fp as FqTrait;

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
            4 => {
                *P3 = *P;
                Self::xdbl_proj(A24, C24, &mut P3.X, &mut P3.Z);
                Self::xdbl_proj(A24, C24, &mut P3.X, &mut P3.Z);
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

    /// Compute (X + Z) and (X - Z) for [i]P for 0 < i <= (ell - 1) / 2
    fn edwards_multiples(A24: &Fq, C24: &Fq, constants: &mut [(Fq, Fq)], P: &PointX<Fq>) {
        let mut iP = PointXMultiples::new(A24, C24, P);
        for c in constants.iter_mut() {
            let (X, Z) = iP.next().unwrap().coords();
            *c = (X - Z, X + Z)
        }
    }

    fn velu_two_isogeny_proj(kernel: &PointX<Fq>, img_points: &mut [PointX<Fq>]) -> (Fq, Fq) {
        assert!(kernel.X.is_zero() != u32::MAX); // TODO: Handle this edge case

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

        (A24_cod, C24_cod)
    }

    fn velu_odd_isogeny_proj(
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

        // Convert back to Montgomery (A24 : C24)
        // let A_codomain = (A_ed + D_ed).mul2();
        // let C_codomain = A_ed - D_ed;
        *A24 = A_ed.mul4();
        *C24 = (A_ed - D_ed).mul4();
    }

    fn velu_two_power_isogeny_proj(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        len: usize,
        img_points: &mut [PointX<Fq>],
    ) {
        let mut eval_points = img_points.to_vec();
        eval_points.push(*kernel);

        // TODO: use a balanced strategy?
        for i in 0..len {
            // Compute [2^(len - i - 1)]
            let mut ker_step = *eval_points.last().unwrap();
            Self::xdbl_proj_iter(A24, C24, &mut ker_step, len - i - 1);

            // Compute the 2-isogeny
            (*A24, *C24) = Self::velu_two_isogeny_proj(&ker_step, &mut eval_points);
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);
    }

    fn velu_odd_power_isogeny_proj(
        A24: &mut Fq,
        C24: &mut Fq,
        kernel: &PointX<Fq>,
        degree: usize,
        len: usize,
        img_points: &mut [PointX<Fq>],
    ) {
        let mut eval_points = img_points.to_vec();
        eval_points.push(*kernel);

        // TODO: use a balanced strategy?
        for i in 0..len {
            // Compute [ell^(len - i - 1)] K which is a point of order ell
            let mut ker_step = *eval_points.last().unwrap();
            for _ in 0..(len - i - 1) {
                ker_step = Self::xmul_proj_u64_vartime(A24, C24, &ker_step, degree as u64);
            }
            Self::velu_odd_isogeny_proj(A24, C24, &ker_step, degree, &mut eval_points);
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);
    }

    pub fn velu_composite_isogeny_proj(
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
                for _ in 0..*e {
                    ker_step = Self::xmul_proj_u64_vartime(A24, C24, &ker_step, *p as u64);
                }
            }

            // Handle 2^n isogenies separately.
            if ell == 2 {
                Self::velu_two_power_isogeny_proj(A24, C24, &ker_step, n, &mut eval_points)
            } else {
                Self::velu_odd_power_isogeny_proj(A24, C24, &ker_step, ell, n, &mut eval_points)
            };
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);
    }

    // ============================================================

    pub fn velu_prime_isogeny(
        self,
        kernel: &PointX<Fq>,
        degree: usize,
        img_points: &mut [PointX<Fq>],
    ) -> Self {
        // 2-isogenies are handled with a special function
        if degree == 2 {
            let (A, C) = Self::velu_two_isogeny_proj(kernel, img_points);
            Self::new(&(A / C))
        } else {
            let mut A24 = self.A + Fq::TWO;
            let mut C24 = Fq::FOUR;

            Self::velu_odd_isogeny_proj(&mut A24, &mut C24, kernel, degree, img_points);
            Self::curve_from_A24_proj(&A24, &C24)
        }
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

        if degree == 2 {
            Self::velu_two_power_isogeny_proj(&mut A24, &mut C24, kernel, len, img_points)
        } else {
            Self::velu_odd_power_isogeny_proj(&mut A24, &mut C24, kernel, degree, len, img_points)
        };

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
