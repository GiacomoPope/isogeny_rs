// ========================================================
// Compting the gluing (2,2)-isogeny from a product of
// elliptic curves to a level 2 theta structure
//
// A lot of the code below is to compute a 4x4 matrix,
// represented as an array [Fq; 16] to compute a symplectic
// basis transformation to for the points into a form
// compatible with the isogeny formula
// ========================================================

use super::elliptic_product::{EllipticProduct, ProductPoint};
use super::theta_point::ThetaPoint;
use super::theta_structure::ThetaStructure;
use super::theta_util::{apply_base_change, to_hadamard};
use crate::elliptic::curve::Curve;
use crate::elliptic::projective_point::Point;

use fp2::fq::Fq as FqTrait;

impl<Fq: FqTrait> EllipticProduct<Fq> {
    /// Given a point in the four torsion, compute the 2x2 matrix needed
    /// for the basis change
    /// M = (1/lambda) * [[a, b], [c, d]] represented as an array [a, b, c, d] and
    /// a scalar lambda.
    /// Cost: 9M + 3S
    fn action_by_translation(E: &Curve<Fq>, T: &Point<Fq>) -> ((Fq, Fq, Fq, Fq), Fq) {
        let (x, z) = T.to_xz();
        let (u, w) = E.xdbl_coords(&x, &z); // Cost 3M 2S

        // Precompute some pieces
        let wx = w * x;
        let uz = u * z;
        let det = wx - uz;

        let lambda = z * det;

        // Compute the matrix coefficients
        // 1/(D * Z) [[-UZ^2, WZ^2, UXZ - XD, UZ^2]]
        let d = uz * z;
        let a = -d;
        let b = -w * z.square();
        let c = x * (uz - det);

        ((a, b, c, d), lambda)
    }

    /// Given the four torsion below the isogeny kernel, compute the
    /// compatible symplectic transform to allow the theta points to have
    /// a good representation for the gluing isogeny
    ///
    /// Input is expected to be K1 = (P1, P2), K2 = (Q1, Q2) in E1 x E2
    /// Inside (E1 x E2)[4].
    /// Cost 104M + 8S
    fn get_base_matrix(&self, P1P2: &ProductPoint<Fq>, Q1Q2: &ProductPoint<Fq>) -> [Fq; 16] {
        // First compute the submatrices from each point
        let (E1, E2) = self.curves();
        let (P1, P2) = P1P2.points();
        let (Q1, Q2) = Q1Q2.points();

        // NOTES: action_by_translation computes a matrix G together with a scaling factor
        // such that the coefficients we want are given by G / d.
        //
        // The idea is we want to avoid division, so we compute the matrix M with a projective
        // factor (which is fine, as we act my matrix multiplication on ThetaPoints which have
        // a projective representation).
        // Naturally, we find that the denominator of each coefficent is:
        //   - a0, b0, c0, d0 is dg_1 * dg_2 * dh_1 * dh_2
        //   - a1, b1, c1, dg is dg_1 * dg_2 * dh_1 * dh_2^2
        //   - a2, b2, c2, dh is dg_1^2 * dg_2 * dh_1 * dh_2
        //   - a3, b3, c3, d3 is dg_1^2 * dg_2 * dh_1 * dh_2^2
        // We require a common denominator for all components so we need to further
        // scale the 12 coefficients by
        //   - `a0, b0, c0, d0` by dg_1 * dh_2
        //   - `a1, b1, c1, dg` by dg_1
        //   - `a2, b2, c2, dh` by dh_2
        //
        // We can do a single inversion during action by translation for a cost of
        // 21M + 4*(2S + 11M) + 1I to compute the coefficients, then the coefficients below
        // cost 44M to compute. If we model 1I to be log(p) multiplications then this
        // is ~100M for level 1 bringing the total to about 209 M.
        //
        // If we don't do any inversions then four action_by_translation costs 4*(9M + 3S)
        // but then the computation of the coefficients grows. For example, the coefficient
        // a0 is given by:
        //    a0 = 1 + g00_1 * g00_2 + h00_1 * h00_2 + gh00_1 * gh00_2
        // with the denominator explicit this is:
        //    a0 = 1 + g00_1/dg_1 * g00_2/dg_2 + h00_1/dh_1 * h00_2/dh_2 + gh00_1 * gh00_2 / (dg_1 * dg_2 * dh_1 * dh_2)
        // and clearing the denominator we then have to compute the numerator as
        //    a0 = (dg_1 * dg_2 * dh_1 * dh_2
        //      += g00_1 * g00_2 * dh_1 * dh_2
        //      += h00_1 * h00_2 * dg_1 * dg_2
        //      += gh00_1 * gh00_2) * dg_1 * dh_2
        // Most of the work has to happen in a0, b0, c0 and d0 and the best I have so far
        // is a cost of 24M to compute all the coefficients projectively, bringing the total
        // cost to 44 + 24 = 68. Adding it all together the new cost is 104M + 12S, which
        // should be about half what it used to be for level 1 and even better for higher
        // levels where inversions cost more.

        // First we compute action by translation for each of the four points of order
        // four. This is done projectively with the denominators in dg_i and dh_i.
        // Cost: 4 x (9M + 3S) = 36M + 12S
        let ((g00_1, g01_1, g10_1, g11_1), dg_1) = Self::action_by_translation(&E1, &P1);
        let ((g00_2, g01_2, g10_2, g11_2), dg_2) = Self::action_by_translation(&E2, &P2);
        let ((h00_1, _, h10_1, _), dh_1) = Self::action_by_translation(&E1, &Q1);
        let ((h00_2, h01_2, h10_2, h11_2), dh_2) = Self::action_by_translation(&E2, &Q2);

        // Compute the product of g1 * h1 and g2 * h2 as 2x2 matricies
        // and extract out the first column

        // first col of g1 * h1 = [[gh00_1, *], [gh10_1, *]]
        let gh00_1 = g00_1 * h00_1 + g01_1 * h10_1;
        let gh10_1 = g10_1 * h00_1 + g11_1 * h10_1;

        // first col of g2 * h2 = [[gh00_2, *], [gh10_2, *]]
        let gh00_2 = g00_2 * h00_2 + g01_2 * h10_2;
        let gh10_2 = g10_2 * h00_2 + g11_2 * h10_2;

        // start the trace with the (projective) identity
        let dg_12 = dg_1 * dg_2;
        let dh_12 = dh_1 * dh_2;

        let mut a = dg_12 * dh_12;
        let mut b = Fq::ZERO;
        let mut c = Fq::ZERO;
        let mut d = Fq::ZERO;

        // T1
        a += g00_1 * g00_2 * dh_12;
        b += g00_1 * g10_2 * dh_12;
        c += g10_1 * g00_2 * dh_12;
        d += g10_1 * g10_2 * dh_12;

        // T2
        a += h00_1 * h00_2 * dg_12;
        b += h00_1 * h10_2 * dg_12;
        c += h10_1 * h00_2 * dg_12;
        d += h10_1 * h10_2 * dg_12;

        // T1+T2
        a += gh00_1 * gh00_2;
        b += gh00_1 * gh10_2;
        c += gh10_1 * gh00_2;
        d += gh10_1 * gh10_2;

        // Now we act by (0, Q2)
        let mut a1 = h00_2 * a + h01_2 * b;
        let mut b1 = h10_2 * a + h11_2 * b;
        let mut c1 = h00_2 * c + h01_2 * d;
        let mut d1 = h10_2 * c + h11_2 * d;

        // Now we act by (P1, 0)
        let mut a2 = g00_1 * a + g01_1 * c;
        let mut b2 = g00_1 * b + g01_1 * d;
        let mut c2 = g10_1 * a + g11_1 * c;
        let mut d2 = g10_1 * b + g11_1 * d;

        // Now we act by (P1, Q2)
        let a3 = g00_1 * a1 + g01_1 * c1;
        let b3 = g00_1 * b1 + g01_1 * d1;
        let c3 = g10_1 * a1 + g11_1 * c1;
        let d3 = g10_1 * b1 + g11_1 * d1;

        // Ensure that the denominator for all coefficients is dg_1^2 * dg_2 * dh_1 * dh_2^2
        let tmp = dg_1 * dh_2;
        a *= tmp;
        b *= tmp;
        c *= tmp;
        d *= tmp;

        a1 *= dg_1;
        b1 *= dg_1;
        c1 *= dg_1;
        d1 *= dg_1;

        a2 *= dh_2;
        b2 *= dh_2;
        c2 *= dh_2;
        d2 *= dh_2;

        // 68M for coefficient computation
        [a, b, c, d, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3]
    }

    /// Given a couple point as input, compute the corresponding ThetaPoint on
    /// the level two structure and then apply the basis change on this point
    /// Cost: 20M
    fn base_change_couple_point(P1P2: &ProductPoint<Fq>, M: &[Fq; 16]) -> ThetaPoint<Fq> {
        let (P1, P2) = P1P2.points();
        let (mut X1, mut Z1) = P1.to_xz();
        let (mut X2, mut Z2) = P2.to_xz();

        // If we have the point (0, 0) swap to (1, 0)
        let P1_check = X1.is_zero() & Z1.is_zero();
        X1.set_cond(&Fq::ONE, P1_check);
        Z1.set_cond(&Fq::ZERO, P1_check);

        // If we have the point (0, 0) swap to (1, 0)
        let P2_check = X2.is_zero() & Z2.is_zero();
        X2.set_cond(&Fq::ONE, P2_check);
        Z2.set_cond(&Fq::ZERO, P2_check);

        // Take all products to get level-2 theta point
        let X = X1 * X2;
        let Y = X1 * Z2;
        let Z = Z1 * X2;
        let T = Z1 * Z2;
        let mut P = ThetaPoint::from_coords(&X, &Y, &Z, &T);

        // Finally apply the base change on the point
        apply_base_change(&mut P, M);
        P
    }

    /// For a theta point which is on an elliptic product,
    /// one of the dual coordinates will be zero. We need
    /// to identify the index of this zero element for the
    /// gluing isogeny codomain and evaluation functions
    fn zero_index(dual_coords: &[Fq; 4]) -> usize {
        let mut z_idx = 0;
        for (i, el) in dual_coords.iter().enumerate() {
            // When el is zero, the result is 0xFF...FF
            // and zero otherwise, so we can use this as
            // a mask for each step.
            let el_is_zero = el.is_zero();
            z_idx |= i as u32 & el_is_zero;
        }
        z_idx as usize
    }

    /// Given the 8-torsion above the kernel of order 2, computes the
    /// codomain ThetaStructure (2,2)-isogenous from a product of elliptic
    /// curves
    ///
    /// NOTE: this function is a little fussy as we need to avoid the
    /// zero-dual coordinate. There's a chance refactoring this could make
    /// it appear more friendly
    ///
    /// Cost: 8S 4M
    fn gluing_codomain(
        T1_8: &ThetaPoint<Fq>,
        T2_8: &ThetaPoint<Fq>,
    ) -> (ThetaStructure<Fq>, Fq, Fq, Fq) {
        // First construct the dual coordinates of the kernel and look
        // for the element which is zero
        // For convenience we pack this as an array instead of a tuple:
        let (xA, xB, yC, yD) = T1_8.squared_theta();
        let (zA, _, zC, _) = T2_8.squared_theta();

        // One element for each array above will be zero. Due to the method we
        // use for gluing, this index is always 3, allowing a simplification of
        // the below code.
        assert!(Self::zero_index(&[xA, xB, yC, yD]) == 3);

        // Codomain coefficients are (x * z * A) \star (A : B : C : 0)
        let mut A = xA * zA;
        let mut B = xB * zA;
        let mut C = zC * xA;
        let mut D = Fq::ZERO;

        // We also want to have (A^-1 : B^-1 : C^-1 : -) projectively for inverses.
        // as (x * z * A * B * C) \star (A^-1 : B^-1 : C^-1 : -)
        let A_inv = xB * zC;
        let B_inv = C;
        let C_inv = B;

        (A, B, C, D) = to_hadamard(&A, &B, &C, &D);
        let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

        (codomain, A_inv, B_inv, C_inv)
    }

    /// Given a point and it's shifted value, compute the image of this
    /// point using the Fq elements precomputed during the codomain
    /// computation
    /// Note: The shift value is needed as one element of the dual coordinates
    /// is zero, so we take values from both and use linear algebra to recover
    /// the correct image.
    ///
    /// Cost: 8S + 9M
    fn gluing_theta_image(
        T: &ThetaPoint<Fq>,
        T_shift: &ThetaPoint<Fq>,
        A_inv: &Fq,
        B_inv: &Fq,
        C_inv: &Fq,
    ) -> ThetaPoint<Fq> {
        // Find dual coordinates of point to push through
        let (Ax, By, Cz, Dt) = T.squared_theta();

        // We are in the case where at most one of A, B, C, D is
        // zero, so we need to account for this
        // To recover values, we use the translated point to get
        let (Ay, _, Ct, Dz) = T_shift.squared_theta();

        // In the general case, we can always directly compute three elements
        // but the fourth requires scalaing by something non-zero and in
        // the original (2,2) paper this was kep general.
        //
        // However, because of the way we perform gluing, the zero index is
        // always three, and so we know where the zero is for (A : B : C : 0).
        assert!(Dt.is_zero() == u32::MAX);
        assert!(Dz.is_zero() == u32::MAX);

        // Despite this, we still need to handle the case where `t` itself
        // could be zero. We  first compute x, y, z from the precomputations
        // during codomain directly.
        let mut x = Ax * (*A_inv);
        let mut y = By * (*B_inv);
        let mut z = Cz * (*C_inv);

        // Now to compute `t` we need to compute a projective factor using that
        // Ay will not be zero for our implementation:
        assert!(Ay.is_zero() != u32::MAX);

        // Compute t with y as a projective factor
        let mut t = y * Ct * (*C_inv);

        // Scale x, y and z to ensure projective factors for t match.
        let lam = Ay * (*A_inv);
        x *= lam;
        y *= lam;
        z *= lam;

        // Perform the final Hadamard
        (x, y, z, t) = to_hadamard(&x, &y, &z, &t);

        ThetaPoint::from_coords(&x, &y, &z, &t)
    }

    /// Compute the image of a ProductPoint through a gluing isogeny given
    /// the point `P`, a point of 4-torsion `K4`, the change of basis matrix
    /// `M` and three precomputed elements of Fq from the codomain computation.
    ///
    /// Cost:
    ///   - 2 * (16M + 5S) for the ProductPoint addition
    ///   - 2 * 20M for the base change
    ///   - 9M + 8S for the gluing theta image
    ///
    /// Total: 81M + 18S per point
    fn gluing_eval(
        self,
        P: &ProductPoint<Fq>,
        K4: &ProductPoint<Fq>,
        M: &[Fq; 16],
        A_inv: &Fq,
        B_inv: &Fq,
        C_inv: &Fq,
    ) -> ThetaPoint<Fq> {
        // Need affine coordinates here to do an add, if we didn't we
        // could use faster x-only... Something to think about but no
        // obvious solution.
        let P_sum_T = self.add(P, K4);

        // After we have added the points, we can use the gluing formula
        // to recover theta points on the level 2 theta structure. First we
        // must compute the basis change as we did for the kernel:
        let T = Self::base_change_couple_point(P, M);
        let T_shift = Self::base_change_couple_point(&P_sum_T, M);

        // With a point and the shift value from the kernel, we can find
        // the image
        Self::gluing_theta_image(&T, &T_shift, A_inv, B_inv, C_inv)
    }

    /// Compute the gluing (2,2)-isogeny from a ThetaStructure computed
    /// from an elliptic product.
    /// Cost for codomain: 172*M + 40*S
    ///   - Compute the four-torsion: 4 * (6M + 6S)
    ///   - Computing the base matrix: 104M + 8S
    ///   - Computing the base change: 2 * 20M
    ///   - Computing the codomain: 8S + 4M
    ///
    /// Cost for images: 81M + 18S.
    pub fn gluing_isogeny(
        &self,
        P1P2_8: &ProductPoint<Fq>,
        Q1Q2_8: &ProductPoint<Fq>,
        eval_points: &[ProductPoint<Fq>],
        image_points: &mut [ThetaPoint<Fq>],
    ) -> ThetaStructure<Fq> {
        // First recover the four torsion below the 8 torsion
        // Cost: 4 * (6M + 6S)
        let P1P2_4 = self.double(P1P2_8);
        let Q1Q2_4 = self.double(Q1Q2_8);

        // Use the four torsion to deterministically find basis change
        // Cost: 104M + 8S
        let M = self.get_base_matrix(&P1P2_4, &Q1Q2_4);

        // Take the points P1, P2 in E1 x E2 and represent them
        // as a theta point on a level 2 structure and map them
        // through the above basis change.
        // Cost: 2 * 20M
        let T1_8 = Self::base_change_couple_point(P1P2_8, &M);
        let T2_8 = Self::base_change_couple_point(Q1Q2_8, &M);

        // Now it's time to compute the codomain and image of the isogeny
        // with kernel below T1, and T2.
        // We save the zero index, as we use it for the images, and we also
        // can precompute a few inverses to save time for evaluation.
        // Cost: 8S + 4M
        let (codomain, A_inv, B_inv, C_inv) = Self::gluing_codomain(&T1_8, &T2_8);

        // Push points through the gluing isogeny.
        // Cost: 81M + 18S per ProductPoint evaluated.
        for (i, P) in eval_points.iter().enumerate() {
            image_points[i] = self.gluing_eval(P, &P1P2_4, &M, &A_inv, &B_inv, &C_inv);
        }

        codomain
    }
}
