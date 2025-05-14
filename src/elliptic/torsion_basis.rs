use super::basis::BasisX;
use super::curve::Curve;
use super::point::PointX;

use fp2::fq::Fq as FqTrait;

impl<Fq: FqTrait> Curve<Fq> {
    /// TODO:
    fn projective_difference(self, P: &PointX<Fq>, Q: &PointX<Fq>) -> PointX<Fq> {
        let mut t0 = P.X * Q.X;
        let mut t1 = P.Z * Q.Z;
        let mut bxx = t0 - t1;
        bxx.set_square(); // bxx = (PX * QX - PZ * QZ)^2
        let mut bxz = t0 + t1;
        t0 = P.X * Q.Z;
        t1 = P.Z * Q.X;
        let mut bzz = t0 + t1;
        bxz *= bzz; // (PX * QX + PZ * QZ) * (PX * QZ + PZ * QX)
        bzz = t0 - t1;
        bzz.set_square(); // bzz = (PX * QZ + PZ * QX)^2
        t0 *= t1;
        t0 *= self.A;
        t0.set_mul2();
        bxz += t0; // bzz = (PX * QX + PZ * QZ) * (PX * QZ + PZ * QX) + 2 * A * PX * QZ * PZ * QX

        // We now normalise the result by \bar{(PZ * QZ)}^2
        // TODO: I can remove a square here, right?
        t0 = P.Z.conjugate();
        t0.set_square();
        t1 = Q.Z.conjugate();
        t1.set_square();
        t0 *= t1;
        bxx *= t0;
        bxz *= t0;
        bzz *= t0;

        // Solve the quadratic equation
        t0 = bxz.square();
        t1 = bxx * bzz;
        t0 -= t1;
        let r = t0.set_sqrt();
        assert!(r == u32::MAX);

        // Set the point (bxz + sqrt(bxz^2 - bxx*bzz) : bzz)
        PointX::new(&(bxz + t0), &bzz)
    }

    /// Given the Montgomery coefficient A which is not a square, we can
    /// compute xP as ([n]*A : 1) for some [n]
    fn full_even_torsion_point_from_A(self) -> (PointX<Fq>, u8) {
        let mut x = self.A;
        let mut h: usize = 1;
        while self.is_on_curve(&x) == 0 {
            x += self.A;
            h += 1;
        }
        let xP = PointX::new(&x, &Fq::ONE);

        // If h fits in 7 bits, we use it as a hint, otherwise we fallback
        // to zero to signify failure to fit the hint within the allocated
        // size.
        let hint = if h < 128 { h as u8 } else { 0 };

        (xP, hint)
    }

    /// When the Montgomery coefficient A is a square, we find an element
    /// of Fp2 of the form 1 + i*b for which -A / (1 + i*b) is a valid
    /// x-coordinate on the curve and (1 + b^2) is a non-quadratic residue
    /// in the base field GF(p).
    fn full_even_torsion_point_from_nqr(self) -> (PointX<Fq>, u8) {
        // We need to find a b such that (1 + b^2) is a not a square and
        // x = - A / (1 + i* b) is a valid x-coordinate on the curve.
        let mut is_nqr: u32 = 0;
        let mut on_curve: u32 = 0;

        let mut b = Fq::ZERO;
        let mut z = Fq::ONE;

        let mut h: u16 = 0;
        while on_curve == 0 {
            // Find 1 + h^2 which is not a square in GF(p)
            while is_nqr == 0 {
                b.set_x0_small((h * h + 1) as i32);
                is_nqr = !b.is_square_base_field();
                h += 1;
            }

            // We now need to determine whether -A / (1 + i*h)
            // is a point on the curve. This is the same as checking
            // whether A^2 * (z - 1) - z^2 is a non-square in GF(p^2)
            z = Fq::ONE;
            let mut t0 = Fq::ZERO;
            t0.set_x1_small((h - 1) as i32);
            z.set_x1_small((h - 1) as i32);

            // z is a valid coordinate providing that (A^2 * (z - 1) - z^2) is a NQR
            let t1 = self.A.square() * t0 - z.square();
            on_curve = !t1.is_square();

            // Reset is_nqr to find a new value on failure of the above.
            is_nqr = 0;
        }

        // Create a point
        let x = -self.A / z;
        let xP = PointX::new(&x, &Fq::ONE);

        // If h fits in 7 bits, we use it as a hint, otherwise we fallback
        // to zero to signify failure to fit the hint within the allocated
        // size.
        let hint = if h < 128 { (h - 1) as u8 } else { 0 };

        (xP, hint)
    }

    /// Compute the torsion basis E\[2^e\] together with a hint for use with
    /// the method `torsion_basis_2e_from_hint()` following the method in the
    /// SQIsign spec. For a curve with order c * 2^f, `e_diff` should be `f - e`
    /// and `cofactor` and `cofactor_bitsize` should encode `c` as an array of
    /// bytes.
    ///
    /// Warning: requires that self.A != 0
    pub fn torsion_basis_2e_with_hint(
        self,
        e_diff: usize,
        cofactor: &[u8],
        cofactor_bitsize: usize,
    ) -> (BasisX<Fq>, u8) {
        // Whether or not A is square determines how we sample P.
        // We assume that the curve itself is public, and so we can
        // branch on this outcome.
        let A_is_square: u8 = (self.A.is_square() & 1) as u8;

        // Compute the point P which has full even order 2^f.
        let (mut xP, hint_P) = if A_is_square == 0 {
            self.full_even_torsion_point_from_A()
        } else {
            self.full_even_torsion_point_from_nqr()
        };

        println!("xP: {}", xP);

        // We can compute another linearly independent point from xP
        // which also has full even torsion.
        let mut xPQ = PointX::new(&-(xP.X + self.A), &Fq::ONE);

        // First we clear the odd cofactor using multiplication
        xP = self.xmul(&xP, cofactor, cofactor_bitsize);
        xPQ = self.xmul(&xPQ, cofactor, cofactor_bitsize);

        // We clear the 2^(f - e) order with repeated doubling.
        xP = self.xmul_2e(&xP, e_diff);
        xPQ = self.xmul_2e(&xPQ, e_diff);

        // Compute the difference point to get the basis x(P), x(Q) and x(P-Q)
        let xQ = self.projective_difference(&xP, &xPQ);

        // Set the basis and the output hint of the form hint_P | hint_A
        let basis = BasisX::from_array([xP, xQ, xPQ]);
        let hint = (hint_P << 1) | A_is_square;

        (basis, hint)
    }

    /// Compute the torsion basis E\[2^e\] using a `hint` to speed up the
    /// computation of points `P` of full even order following the method in the
    /// SQIsign spec. For a curve with order c * 2^f, `e_diff` should be `f - e`
    /// and `cofactor` and `cofactor_bitsize` should encode `c` as an array of
    /// bytes.
    ///
    /// Warning: requires that self.A != 0
    pub fn torsion_basis_2e_from_hint(
        self,
        e_diff: usize,
        cofactor: &[u8],
        cofactor_bitsize: usize,
        hint: u8,
    ) -> BasisX<Fq> {
        // We can recover whether A is a square, and a value for xP from
        // the hint
        let a_is_square = hint & 1;
        let hint_P = hint >> 1;

        // Compute (xP : 1) from the hint (fall back method for very rare failure)
        let (mut xP, _) = if hint_P == 0 {
            // If the hint is zero, we need to fall back to the usual method
            if a_is_square == 0 {
                self.full_even_torsion_point_from_A()
            } else {
                self.full_even_torsion_point_from_nqr()
            }
        } else {
            let x = if a_is_square == 0 {
                // When A is not a square, x = A*hint_P
                self.A.mul_small(hint_P as i32)
            } else {
                // Otherwise x = 1 + i * hint_P
                // let mut t = Fq::ONE;
                // t.set_x1_small(hint_P as i32);
                let t = Fq::from_i32_pair(1, hint_P as i32);
                -self.A / t
            };
            (PointX::new(&x, &Fq::ONE), 0)
        };

        println!("xP: {}", xP);

        // We can compute another linearly independent point from xP
        // which also has full even torsion.
        let mut xPQ = PointX::new(&-(xP.X + self.A), &Fq::ONE);

        // First we clear the odd cofactor using multiplication
        xP = self.xmul(&xP, cofactor, cofactor_bitsize);
        xPQ = self.xmul(&xPQ, cofactor, cofactor_bitsize);

        // We clear the 2^(f - e) order with repeated doubling.
        xP = self.xmul_2e(&xP, e_diff);
        xPQ = self.xmul_2e(&xPQ, e_diff);

        // Compute the difference point to get the basis x(P), x(Q) and x(P-Q)
        let xQ = self.projective_difference(&xP, &xPQ);
        let basis = BasisX::from_array([xP, xQ, xPQ]);

        basis
    }
}
