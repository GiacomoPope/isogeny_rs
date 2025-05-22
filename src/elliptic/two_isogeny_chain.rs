use fp2::traits::Fp2 as FqTrait;

use super::{curve::Curve, point::PointX};

impl<Fq: FqTrait> Curve<Fq> {
    /// Compute a curve from the projective coordinates of (A + 2) / 4 = (A24 : C24)
    #[inline]
    fn curve_from_A24_proj(A24: &Fq, C24: &Fq) -> Curve<Fq> {
        // Compute A from (A24 : C24)
        let mut A = (*A24) + (*A24);
        A -= *C24;
        A += A;
        A /= *C24;

        Curve::new(&A)
    }

    /// Compute [2]P in place using projective (A + 2) / 4 = (A24 : C24)
    /// Cost: 2S + 4M
    #[inline(always)]
    fn xdbl_proj(A24: &Fq, C24: &Fq, P: &mut PointX<Fq>) {
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
    fn xdbl_proj_iter(A24: &Fq, C24: &Fq, P: &mut PointX<Fq>, n: usize) {
        for _ in 0..n {
            Self::xdbl_proj(A24, C24, P);
        }
    }

    /// Compute the codomain of the 2-isogeny E -> E/<ker> for ker != (0 : 1)
    fn two_isogeny_codomain(ker: &PointX<Fq>) -> (Fq, Fq) {
        let mut A24 = ker.X.square();
        let C24 = ker.Z.square();
        A24 = C24 - A24;
        (A24, C24)
    }

    /// Evaluate a point Q in place under the action of the 2-isogeny E -> E/<ker>
    /// for ker != (0 : 1)
    fn two_isogeny_eval(ker: &PointX<Fq>, Q: &mut PointX<Fq>) {
        let mut t0 = ker.X + ker.Z;
        let mut t1 = ker.X - ker.Z;
        let mut t2 = Q.X + Q.Z;
        let mut t3 = Q.X - Q.Z;

        t0 *= t3;
        t1 *= t2;
        t2 = t0 + t1;
        t3 = t0 - t1;

        Q.X *= t2;
        Q.Z *= t3;
    }

    /// Compute the codomain of the 2-isogeny E -> E/<ker> for ker == (0 : 1)
    fn two_isogeny_codomain_singular(A24: &mut Fq, C24: &mut Fq) -> (Fq, Fq) {
        let mut t0 = A24.mul2();
        t0 -= *C24;
        t0.set_mul2();
        t0 /= *C24;
        let c0 = t0;
        *A24 = t0.mul2();
        t0.set_square();
        t0 -= Fq::FOUR;
        let r = t0.set_sqrt();
        assert!(r == u32::MAX); // TODO
        let c1 = -t0;
        *C24 = t0.mul2();
        *A24 += *C24;
        C24.set_mul2();

        (c0, c1)
    }

    /// Evaluate a point Q in place under the action of the 2-isogeny E -> E/<ker>
    /// for ker = (0 : 1)
    fn two_isogeny_eval_singular(c0: &Fq, c1: &Fq, Q: &mut PointX<Fq>) {
        let t0 = Q.X * Q.Z;
        let mut t1 = (*c0) * Q.Z;
        t1 += Q.X;
        t1 *= Q.X;

        Q.X = Q.Z.square();
        Q.X += t1;
        Q.Z = t0 * (*c1);
    }

    /// Compute the codomain of the 4-isogeny E -> E/<ker> for [2]ker != (0 : 1)
    /// Returns the codomain (A24 : C24) together with three constants (c0, c1, c1)
    /// used for computing images.
    fn four_isogeny_codomain(ker: &PointX<Fq>) -> (Fq, Fq, Fq, Fq, Fq) {
        let mut c0 = ker.Z.square();
        let c1 = ker.X - ker.Z;
        let c2 = ker.X + ker.Z;
        let t0 = ker.X.square();
        let t1 = c0 + t0;
        let t2 = c0 - t0;
        let A24 = t1 * t2;
        let C24 = c0.square();
        c0.set_mul4();

        (A24, C24, c0, c1, c2)
    }

    /// Evaluate a point Q in place under the action of the 4-isogeny E -> E/<ker>
    /// for [2]ker != (0 : 1)
    fn four_isogeny_eval(c0: &Fq, c1: &Fq, c2: &Fq, Q: &mut PointX<Fq>) {
        let mut t0 = Q.X + Q.Z;
        let mut t1 = Q.X - Q.Z;
        Q.X = t0 * (*c1);
        Q.Z = t1 * (*c2);
        t0 *= t1;
        t0 *= *c0;
        t1 = Q.X + Q.Z;
        Q.Z = Q.X - Q.Z;
        t1.set_square();
        Q.Z.set_square();
        Q.X = t0 + t1;
        t0 -= Q.Z;
        Q.X *= t1;
        Q.Z *= t0;
    }

    /// Compute a 2^n isogeny using the naive approach. WARNING: branches on whether
    /// the kernel is of the form (0 : 1) or not, and so is not constant time.
    /// Returns the codomain as well as a `u32` which is equal to `0xFF..FF` on success
    /// or `0x00..00` on failure, which happens when either the kernel is of the wrong
    /// order, or the kernel contained the point (0 : 1) and `allow_singular` is `false`.
    pub fn two_isogeny_chain_small_vartime(
        &self,
        kernel: &PointX<Fq>,
        n: usize,
        images: &mut [PointX<Fq>],
        allow_singular: bool,
    ) -> (Curve<Fq>, u32) {
        let mut A24 = self.A24;
        let mut C24 = Fq::ONE;

        let mut ker: PointX<Fq> = *kernel;
        let mut ker_step: PointX<Fq>;

        for i in 0..n {
            // Double the kernel to get a point of order 2
            ker_step = ker;
            Self::xdbl_proj_iter(&A24, &C24, &mut ker_step, n - i - 1);

            // For the first step, we check the kernel has exactly order 2^n
            // we also allow the first step to be the singular isogeny with
            // kernel (0 : 1) only if `allow_singular` is `true`.
            if i == 0 {
                // First check if the kernel has the correct order.
                let mut inf = ker_step;
                Self::xdbl_proj(&A24, &C24, &mut inf);
                if (!ker_step.Z.is_zero() & inf.Z.is_zero()) != u32::MAX {
                    return (*self, 0);
                }

                // When the kernel is of the form (0 : 1) we need specialised isogenies.
                // These should never occur in some cases and we identify this for failure
                // cases.
                if ker_step.X.is_zero() == u32::MAX {
                    if !allow_singular {
                        return (*self, 0);
                    }
                    // Compute the codomain from ker_step for kernel (0 : 1)
                    let (c0, c1) = Self::two_isogeny_codomain_singular(&mut A24, &mut C24);

                    // Push through the kernel
                    Self::two_isogeny_eval_singular(&c0, &c1, &mut ker);

                    // Push through the points to evaluate
                    for P in images.iter_mut() {
                        Self::two_isogeny_eval_singular(&c0, &c1, P);
                    }
                    continue;
                }
            }

            // Compute the codomain from ker_step
            (A24, C24) = Self::two_isogeny_codomain(&ker_step);

            // Push through the kernel
            Self::two_isogeny_eval(&ker_step, &mut ker);

            // Push through the points to evaluate
            for P in images.iter_mut() {
                Self::two_isogeny_eval(&ker_step, P);
            }
        }
        (Self::curve_from_A24_proj(&A24, &C24), u32::MAX)
    }

    /// Compute a 2^n isogeny using a balanced strategy with 4-isogenies for
    /// every step except the last in the case when n is odd.
    /// Returns the codomain as well as a `u32` which is equal to `0xFF..FF` on success
    /// or `0x00..00` on failure, which happens when the kernel found to have the wrong
    /// order or the kernel is above the "singular" point (0 : 1).
    pub fn two_isogeny_chain(
        &self,
        kernel: &PointX<Fq>,
        n: usize,
        images: &mut [PointX<Fq>],
    ) -> (Curve<Fq>, u32) {
        // For 2-isogenies we represent (A + 2) / 4 projectively as (A24 : C24)
        let mut A24 = self.A24;
        let mut C24 = Fq::ONE;

        // Precompute constants from the codomain at each step for computing images.
        let mut c0;
        let mut c1;
        let mut c2;

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - n.leading_zeros()) as usize;

        // These are a set of points of order 2^i
        let mut stategy_points: Vec<PointX<Fq>> = vec![PointX::INFINITY; space];

        // The values i such that each point in stategy_points has order 2^i
        let mut orders: Vec<usize> = vec![0; space];

        // Initalise the first values for the strategy
        stategy_points[0] = *kernel;
        orders[0] = n;

        // Value to determine success / failure of isogeny chain
        let mut ok = u32::MAX;

        let mut k = 0;
        for i in 0..(n >> 1) {
            // Get the next point of order 4
            while orders[k] != 2 {
                k += 1;
                let m = 2 * (orders[k - 1] / 4) + (orders[k - 1] & 1);
                stategy_points[k] = stategy_points[k - 1];
                Self::xdbl_proj_iter(&A24, &C24, &mut stategy_points[k], m);
                orders[k] = orders[k - 1] - m;
            }
            // Point of order four to compute isogeny with
            let ker_step = stategy_points[k];

            // For the first step we perform a check that the kernel has the
            // exact order and that the kernel is not above the point (0 : 1)
            if i == 0 {
                let mut tmp = ker_step;

                // Ensure that the [2]ker is not (0 : 1)
                Self::xdbl_proj(&A24, &C24, &mut tmp);
                ok &= !tmp.X.is_zero();

                // Ensure that the kernel has exact order
                // [2]ker != 0 and [4]ker = 0
                ok &= !tmp.Z.is_zero();
                Self::xdbl_proj(&A24, &C24, &mut tmp);
                ok &= tmp.Z.is_zero();
            }

            // Compute the codomain from the current step
            (A24, C24, c0, c1, c2) = Self::four_isogeny_codomain(&ker_step);

            // Push through the kernel points and reduce the stored order
            for i in 0..space {
                Self::four_isogeny_eval(&c0, &c1, &c2, &mut stategy_points[i]);
                orders[i] = orders[i].saturating_sub(2);
            }

            // Push through the points to evaluate
            for P in images.iter_mut() {
                Self::four_isogeny_eval(&c0, &c1, &c2, P);
            }
            k = k.saturating_sub(1);
        }

        // The chain had odd length, so we need to do one final 2-isogeny to finish
        // the chain.
        if n & 1 == 1 {
            // Point of order two to finish isogeny with
            let ker_step = stategy_points[0];

            // Ensure that the [2]ker is not (0 : 1)
            ok &= !ker_step.X.is_zero();

            // Ensure the point has order exactly 2
            let mut tmp = ker_step;
            Self::xdbl_proj(&A24, &C24, &mut tmp);
            ok &= tmp.Z.is_zero();

            // Compute the codomain from ker_step
            (A24, C24) = Self::two_isogeny_codomain(&ker_step);

            // Push through the points to evaluate
            for P in images.iter_mut() {
                Self::two_isogeny_eval(&ker_step, P);
            }
        }

        (Self::curve_from_A24_proj(&A24, &C24), ok)
    }
}
