use fp2::fq::Fq as FqTrait;

use super::{curve::Curve, point::PointX};

impl<Fq: FqTrait> Curve<Fq> {
    /// Compute a curve from the projective coordinates of A^±_{24} = (A + 2C : A - 2C)
    #[inline]
    fn curve_from_A_plus_minus(A24_plus: &Fq, A24_minus: &Fq) -> Curve<Fq> {
        // Compute A from (A + 2C : A - 2C)
        let num = (*A24_plus + *A24_minus).mul2();
        let den = *A24_plus - *A24_minus;
        let A = num / den;

        Curve::new(&A)
    }

    /// Compute \[3\]P given the constants A^±_{24} = (A + 2C : A - 2C) naturally computed during
    /// the 3-isogeny codomain computation.
    /// Cost: 5S + 7M
    #[inline(always)]
    fn xtpl_proj(A24_plus: &Fq, A24_minus: &Fq, P: &mut PointX<Fq>) {
        let mut t0 = P.X - P.Z;
        let mut t2 = t0.square();
        let mut t1 = P.X + P.Z;
        let mut t3 = t1.square();
        let t4 = t1 + t0;
        t0 = t1 - t0;
        t1 = t4.square();
        t1 -= t3;
        t1 -= t2;
        let t5 = t3 * (*A24_plus);
        t3 *= t5;
        let t6 = t2 * (*A24_minus);
        t2 *= t6;
        t3 = t2 - t3;
        t2 = t5 - t6;
        t1 *= t2;
        t2 = t3 + t1;
        t2.set_square();
        P.X = t2 * t4;
        t1 = t3 - t1;
        t1.set_square();
        P.Z = t1 * t0;
    }

    /// Compute \[3^n\]P in place using the constants A^±_{24} = (A + 2C : A - 2C)
    /// Cost: n * (5S + 7M)
    fn xtpl_proj_iter(A24_plus: &Fq, A24_minus: &Fq, P: &mut PointX<Fq>, n: usize) {
        for _ in 0..n {
            Self::xtpl_proj(A24_plus, A24_minus, P);
        }
    }

    /// Given a point P = (XP : ZP) of order 3, computes the
    /// 3-isogeny codomain with coefficient A represented as
    /// A^±_{24} = (A + 2C : A - 2C) along with constants c0, c1
    /// used for computing images.
    /// Cost: 3S + 2M
    #[inline]
    fn three_isogeny_codomain(P: &PointX<Fq>) -> (Fq, Fq, Fq, Fq) {
        let c0 = P.X - P.Z;
        let t0 = c0.square();
        let c1 = P.X + P.Z;
        let t1 = c1.square();
        let mut t3 = P.X + P.X;
        t3.set_square();
        let t2 = t3 - t0;
        t3 -= t1;
        let mut t4 = t0 + t3;
        t4.set_mul2();
        t4 += t1;
        let A24_minus = t2 * t4;
        t4 = t1 + t2;
        t4.set_mul2();
        t4 += t0;
        let A24_plus = t3 * t4;

        (A24_plus, A24_minus, c0, c1)
    }

    /// Given constants (c0, c1) along with the point Q = (XQ : ZQ)
    /// compute the image of this point in place
    /// Cost: 2S + 4M
    #[inline(always)]
    fn three_isogeny_eval(c0: &Fq, c1: &Fq, Q: &mut PointX<Fq>) {
        let mut t0 = Q.X + Q.Z;
        let mut t1 = Q.X - Q.Z;
        t0 *= *c0;
        t1 *= *c1;
        let mut t2 = t0 + t1;
        t0 = t1 - t0;
        t2.set_square();
        t0.set_square();
        Q.X *= t2;
        Q.Z *= t0;
    }

    /// Compute a 3^n isogeny using a balanced strategy
    pub fn three_isogeny_chain(
        self,
        kernel: &PointX<Fq>,
        n: usize,
        images: &mut [PointX<Fq>],
    ) -> Curve<Fq> {
        // For codomain computation we track the constants (A + 2C : A - 2C)
        let mut A24_plus = self.A + Fq::TWO;
        let mut A24_minus = self.A - Fq::TWO;

        // To compute images we need two constants (c0, c1)
        let mut c0;
        let mut c1;

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - n.leading_zeros() + 1) as usize;

        // These are a set of points of order 3^i
        let mut stategy_points: Vec<PointX<Fq>> = vec![PointX::INFINITY; space];

        // The values i such that each point in stategy_points has order 3^i
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
                Self::xtpl_proj_iter(&A24_plus, &A24_minus, &mut stategy_points[k], m);
                orders[k] = orders[k - 1] - m;
            }
            // Point of order two to compute isogeny with
            let ker_step = stategy_points[k];

            // Compute the codomain from ker_step
            (A24_plus, A24_minus, c0, c1) = Self::three_isogeny_codomain(&ker_step);

            // Push through the kernel points and reduce the stored order
            for i in 0..space {
                Self::three_isogeny_eval(&c0, &c1, &mut stategy_points[i]);
                orders[i] = orders[i].saturating_sub(1);
            }

            // Push through the points to evaluate
            for P in images.iter_mut() {
                Self::three_isogeny_eval(&c0, &c1, P);
            }
            k = k.saturating_sub(1);
        }
        Self::curve_from_A_plus_minus(&A24_plus, &A24_minus)
    }
}
