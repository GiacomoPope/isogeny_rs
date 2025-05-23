// TODO:
// Currently after every isogeny the curve is returned, which means computing A / C, we need a new mulx which we can
// use with (A, C) as input instead of A24. This will require a bit of refactoring too for the API etc, but this is
// "roughly" working now.
//
// Another TODO is the cofactor clearing. At the moment clearing the cofactor means multiplying by each ell_i e_i times,
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
    pub fn new(E: &Curve<Fq>, P: &PointX<Fq>) -> Self {
        Self {
            P: *P,
            Q: *P,
            R: E.xdouble(P),
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
    fn edwards_multiples(&self, constants: &mut [(Fq, Fq)], P: &PointX<Fq>) {
        let mut iP = PointXMultiples::new(self, P);
        for c in constants.iter_mut() {
            let (X, Z) = iP.next().unwrap().coords();
            *c = (X - Z, X + Z)
        }
    }

    pub fn velu_two_isogeny(
        &self,
        kernel: &PointX<Fq>,
        img_points: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
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

        Curve::new(&(A_codomain / C_codomain))
    }

    pub fn velu_odd_prime_isogeny(
        &self,
        kernel: &PointX<Fq>,
        degree: usize,
        img_points: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
        // Convert from Montgomery to projective twisted Edwards (A_ed : D_ed)
        let mut A_ed = self.A + Fq::TWO;
        let mut D_ed = self.A - Fq::TWO;

        // We precompute (X - Z) and (X + Z) for (X : Z) = [i]P for i in 0..((ell - 1)/2)
        let d = (degree - 1) >> 1;
        let mut constants = vec![(Fq::ZERO, Fq::ZERO); d];
        self.edwards_multiples(&mut constants, kernel);

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

        // Convert back to Montgomery and return the domain curve.
        let A_codomain = (A_ed + D_ed).mul2();
        let C_codomain = A_ed - D_ed;

        // TODO: if I wrote a new mulx and doublex which uses (A : C)
        // instead of A or A24 then we could keep the curve Montgomery
        // coefficient projective throughout the chain and return
        // (A : C) instead.
        Curve::new(&(A_codomain / C_codomain))
    }

    pub fn velu_prime_power_isogeny(
        self,
        kernel: &PointX<Fq>,
        degree: usize,
        len: usize,
        img_points: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
        let mut codomain = self;
        let mut eval_points = img_points.to_vec();
        eval_points.push(*kernel);

        // TODO: use a balanced strategy?
        for i in 0..len {
            // Compute [ell^(len - i - 1)] K which is a point of order ell
            let mut ker_step = *eval_points.last().unwrap();
            for _ in 0..(len - i - 1) {
                ker_step = codomain.xmul_u64_vartime(&ker_step, degree as u64);
            }
            // 2-isogenies are handled with a special function
            if degree == 2 {
                codomain = codomain.velu_two_isogeny(&ker_step, &mut eval_points);
            } else {
                codomain = codomain.velu_odd_prime_isogeny(&ker_step, degree, &mut eval_points);
            }
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);

        codomain
    }

    pub fn velu_composite_isogeny(
        self,
        kernel: &PointX<Fq>,
        degrees: &[(usize, usize)],
        img_points: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
        let mut codomain = self;
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
                    ker_step = codomain.xmul_u64_vartime(&ker_step, *p as u64);
                }
            }

            // Now ker_step should be a point of order degree^len, we compute a prime-power isogeny
            codomain = codomain.velu_prime_power_isogeny(&ker_step, ell, n, &mut eval_points);
        }

        // TODO: I don't like this copy...
        img_points.copy_from_slice(&eval_points[..eval_points.len() - 1]);

        codomain
    }
}
