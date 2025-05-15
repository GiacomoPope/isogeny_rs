use fp2::fq::Fq as FqTrait;

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

    /// Compute a 2^n isogeny using the naive approach
    pub fn two_isogeny_chain_naive(
        self,
        kernel: &PointX<Fq>,
        n: usize,
        images: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
        let mut A24 = self.A24;
        let mut C24 = Fq::ONE;

        let mut ker: PointX<Fq> = *kernel;
        let mut ker_step: PointX<Fq>;

        for i in 0..n {
            // Double the kernel to get a point of order 2
            ker_step = ker;
            Self::xdbl_proj_iter(&A24, &C24, &mut ker_step, n - i - 1);

            if ker_step.X.is_zero() == u32::MAX {
                // Compute the codomain from ker_step for kernel (0 : 1)
                let (c0, c1) = Self::two_isogeny_codomain_singular(&mut A24, &mut C24);

                // Push through the kernel
                Self::two_isogeny_eval_singular(&c0, &c1, &mut ker);

                // Push through the points to evaluate
                for P in images.iter_mut() {
                    Self::two_isogeny_eval_singular(&c0, &c1, P);
                }
            } else {
                // Compute the codomain from ker_step
                (A24, C24) = Self::two_isogeny_codomain(&ker_step);

                // Push through the kernel
                Self::two_isogeny_eval(&ker_step, &mut ker);

                // Push through the points to evaluate
                for P in images.iter_mut() {
                    Self::two_isogeny_eval(&ker_step, P);
                }
            }
        }
        Self::curve_from_A24_proj(&A24, &C24)
    }

    /// Compute a 2^n isogeny using a balanced strategy
    pub fn two_isogeny_chain(
        self,
        kernel: &PointX<Fq>,
        n: usize,
        images: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
        // For 2-isogenies we represent (A + 2) / 4 projectively as (A24 : C24)
        let mut A24 = self.A24;
        let mut C24 = Fq::ONE;

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - n.leading_zeros() + 1) as usize;

        // These are a set of points of order 2^i
        let mut stategy_points: Vec<PointX<Fq>> = vec![PointX::INFINITY; space];

        // The values i such that each point in stategy_points has order 2^i
        let mut orders: Vec<usize> = vec![0; space];

        // Initalise the first values for the strategy
        stategy_points[0] = *kernel;
        orders[0] = n;

        let mut k = 0;
        for _ in 0..n {
            // Get the next point of order 2
            while orders[k] != 1 {
                k += 1;
                let m = orders[k - 1] / 2;
                stategy_points[k] = stategy_points[k - 1];
                Self::xdbl_proj_iter(&A24, &C24, &mut stategy_points[k], m);
                orders[k] = orders[k - 1] - m;
            }
            // Point of order two to compute isogeny with
            let ker_step = stategy_points[k];

            // Compute the codomain from the current step
            (A24, C24) = Self::two_isogeny_codomain(&ker_step);

            // Push through the kernel points and reduce the stored order
            for i in 0..space {
                Self::two_isogeny_eval(&ker_step, &mut stategy_points[i]);
                orders[i] = orders[i].saturating_sub(1);
            }

            // Push through the points to evaluate
            for P in images.iter_mut() {
                Self::two_isogeny_eval(&ker_step, P);
            }
            k = k.saturating_sub(1);
        }

        Self::curve_from_A24_proj(&A24, &C24)
    }
}
