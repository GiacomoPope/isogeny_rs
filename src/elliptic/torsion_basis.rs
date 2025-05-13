use super::basis::BasisX;
use super::curve::Curve;
use super::point::PointX;

use fp2::fq::Fq as FqTrait;

impl<Fq: FqTrait> Curve<Fq> {
    fn projective_difference(self, xP: &PointX<Fq>, xQ: &PointX<Fq>) -> PointX<Fq> {
        todo!()
    }

    /// Given the Montgomery coefficient A which is not a square, we can
    /// compute xP as ([n]*A : 1) for some [n]
    fn full_even_torsion_point_from_A(self) -> (PointX<Fq>, u8) {
        let mut x = self.A;
        let mut h: usize = 0;
        while self.is_on_curve(&x) == 0 {
            h += 1;
            x += self.A;
        }
        let xP = PointX::new(&x, &Fq::ONE);

        // If h fits in 7 bits, we use it as a hint, otherwise we fallback
        // to zero to signify failure to fit the hint within the allocated
        // size.
        let hint = if h < 128 { h as u8 } else { 0 };

        (xP, hint)
    }

    fn full_even_torsion_point_from_nqr(self) -> (PointX<Fq>, u8) {
        todo!()
    }

    /// Compute the torsion basis E\[2^e\] together with a hint for use with
    /// the method `torsion_basis_2e_from_hint()` following the method in the
    /// SQIsign spec.
    ///
    /// Warning: requires that self.A != 0
    pub fn torsion_basis_2e_with_hint(
        self,
        e: usize,
        cofactor: &[u8],
        cofactor_bitsize: usize,
    ) -> (BasisX<Fq>, u8) {
        // Whether or not A is square determines how we sample P.
        // We assume that the curve itself is public, and so we can
        // branch on this outcome.
        let A_is_square = !((self.A.legendre() as u8) >> 1) & 1;

        let (mut xP, hint_P) = if A_is_square == 0 {
            self.full_even_torsion_point_from_A()
        } else {
            self.full_even_torsion_point_from_nqr()
        };

        // We can compute another linearly independent point from xP
        // which also has full even torsion.
        let mut xPQ = PointX::new(&-(xP.X + self.A), &Fq::ONE);

        // Now we must clear the cofactor to get points of order exactly 2^e
        xP = self.xmul(&xP, cofactor, cofactor_bitsize);
        xPQ = self.xmul(&xPQ, cofactor, cofactor_bitsize);

        // Compute the difference point to get the basis x(P), x(Q) and x(P-Q)
        let xQ = self.projective_difference(&xP, &xPQ);

        // Set the basis and the output hint of the form hint_P | hint_A
        let basis = BasisX::from_array([xP, xQ, xPQ]);
        let hint = (hint_P << 1) | A_is_square;

        (basis, hint)
    }

    /// Compute the torsion basis E\[2^e\] using a `hint` to speed up the
    /// computation of points `P` of full even order following the method in the
    /// SQIsign spec.
    ///
    /// Warning: requires that self.A != 0
    pub fn torsion_basis_2e_from_hint(
        self,
        e: usize,
        cofactor: &[u8],
        cofactor_bitsize: usize,
        hint: u8,
    ) -> BasisX<Fq> {
        // We can recover whether A is a square, and a value for xP from
        // the hint
        let A_is_square = hint & 1;
        let hint_P = hint >> 1;

        // Compute (xP : 1) from the hint (fall back method for very rare failure)
        let (mut xP, _) = if hint_P == 0 {
            // If the hint is zero, we need to fall back to the usual method
            if A_is_square == 0 {
                self.full_even_torsion_point_from_A()
            } else {
                self.full_even_torsion_point_from_nqr()
            }
        } else {
            let x = if A_is_square == 0 {
                // When A is not a square, x = A*hint_P
                self.A.mul_small(hint_P as i32)
            } else {
                // Otherwise x = 1 + i * hint_P
                let mut t = Fq::ONE;
                t.set_x1_small(hint_P as i32);
                -self.A / t
            };
            (PointX::new(&x, &Fq::ONE), 0)
        };

        // We can compute another linearly independent point from xP
        // which also has full even torsion.
        let mut xPQ = PointX::new(&-(xP.X + self.A), &Fq::ONE);

        // Now we must clear the cofactor to get points of order exactly 2^e
        xP = self.xmul(&xP, cofactor, cofactor_bitsize);
        xPQ = self.xmul(&xPQ, cofactor, cofactor_bitsize);
        let xQ = self.projective_difference(&xP, &xPQ);
        let basis = BasisX::from_array([xP, xQ, xPQ]);

        basis
    }
}
