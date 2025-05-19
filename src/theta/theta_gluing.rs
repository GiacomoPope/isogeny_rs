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
    /// M = [[a, b], [c, d]] represented as an array [a, b, c, d]
    /// Cost: 14M + 2S + 1I
    fn get_base_submatrix(E: &Curve<Fq>, T: &Point<Fq>) -> (Fq, Fq, Fq, Fq) {
        let (x, z) = T.to_xz();
        let (u, w) = E.xdbl_coords(&x, &z); // Cost 3M 2S

        // Precompute some pieces
        let wx = w * x;
        let wz = w * z;
        let ux = u * x;
        let uz = u * z;
        let det = wx - uz;

        // Batch inversion
        let mut inverse = [det, z];
        Fq::batch_invert(&mut inverse);

        // Compute the matrix coefficients
        let d = uz * inverse[0]; // Computing d then a saves one negation
        let a = -d;
        let b = -(wz * inverse[0]);
        let c = ux * inverse[0] - x * inverse[1];

        (a, b, c, d)
    }

    /// Given the four torsion below the isogeny kernel, compute the
    /// compatible symplectic transform to allow the theta points to have
    /// a good representation for the gluing isogeny
    ///
    /// Input is expected to be K1 = (P1, P2), K2 = (Q1, Q2) in E1 x E2
    /// Inside (E1 x E2)[4].
    /// Cost 100M + 8S + 4I
    fn get_base_matrix(self, P1P2: &ProductPoint<Fq>, Q1Q2: &ProductPoint<Fq>) -> [Fq; 16] {
        // First compute the submatrices from each point
        let (E1, E2) = self.curves();
        let (P1, P2) = P1P2.points();
        let (Q1, Q2) = Q1Q2.points();

        // TODO: if these were submatrix computations were done together, we
        // could save 3 inversions... It would make the code harder to read
        // but would be an optimisation for the gluing.
        // Something to think about for when cost REALLY matters.
        // Cost: 4 x 14M + 2S + 1I = 56M + 8S + 4I
        let (g00_1, g01_1, g10_1, g11_1) = Self::get_base_submatrix(&E1, &P1);
        let (g00_2, g01_2, g10_2, g11_2) = Self::get_base_submatrix(&E2, &P2);
        let (h00_1, _, h10_1, _) = Self::get_base_submatrix(&E1, &Q1);
        let (h00_2, h01_2, h10_2, h11_2) = Self::get_base_submatrix(&E2, &Q2);

        // Compute the product of g1 * h1 and g2 * h2 as 2x2 matricies
        // and extract out the first column

        // first col of g1 * h1 = [[gh00_1, *], [gh10_1, *]]
        let gh00_1 = g00_1 * h00_1 + g01_1 * h10_1;
        let gh10_1 = g10_1 * h00_1 + g11_1 * h10_1;

        // first col of g2 * h2 = [[gh00_2, *], [gh10_2, *]]
        let gh00_2 = g00_2 * h00_2 + g01_2 * h10_2;
        let gh10_2 = g10_2 * h00_2 + g11_2 * h10_2;

        // start the trace with the identity
        let mut a = Fq::ONE;
        let mut b = Fq::ZERO;
        let mut c = Fq::ZERO;
        let mut d = Fq::ZERO;

        // T1
        a += g00_1 * g00_2;
        b += g00_1 * g10_2;
        c += g10_1 * g00_2;
        d += g10_1 * g10_2;

        // T2
        a += h00_1 * h00_2;
        b += h00_1 * h10_2;
        c += h10_1 * h00_2;
        d += h10_1 * h10_2;

        // T1+T2
        a += gh00_1 * gh00_2;
        b += gh00_1 * gh10_2;
        c += gh10_1 * gh00_2;
        d += gh10_1 * gh10_2;

        // Now we act by (0, Q2)
        let a1 = h00_2 * a + h01_2 * b;
        let b1 = h10_2 * a + h11_2 * b;
        let c1 = h00_2 * c + h01_2 * d;
        let d1 = h10_2 * c + h11_2 * d;

        // Now we act by (P1, 0)
        let a2 = g00_1 * a + g01_1 * c;
        let b2 = g00_1 * b + g01_1 * d;
        let c2 = g10_1 * a + g11_1 * c;
        let d2 = g10_1 * b + g11_1 * d;

        // Now we act by (P1, Q2)
        let a3 = g00_1 * a1 + g01_1 * c1;
        let b3 = g00_1 * b1 + g01_1 * d1;
        let c3 = g10_1 * a1 + g11_1 * c1;
        let d3 = g10_1 * b1 + g11_1 * d1;
        // 44M

        [a, b, c, d, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3]
    }

    /// Given a couple point as input, compute the corresponding ThetaPoint on
    /// the level two structure and then apply the basis change on this point
    /// Cost: 20M
    fn base_change_couple_point(P1P2: &ProductPoint<Fq>, M: [Fq; 16]) -> ThetaPoint<Fq> {
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
    /// Cost: 8S 13M 1I
    fn gluing_codomain(
        T1_8: &ThetaPoint<Fq>,
        T2_8: &ThetaPoint<Fq>,
    ) -> (ThetaStructure<Fq>, (Fq, Fq), usize) {
        // First construct the dual coordinates of the kernel and look
        // for the element which is zero
        // For convenience we pack this as an array instead of a tuple:
        let xAxByCyD: [Fq; 4] = T1_8.squared_theta().into();
        let zAtBzYtD: [Fq; 4] = T2_8.squared_theta().into();

        // One element for each array above will be zero. Identify this
        // element to get the right permutation below for filling arrays
        let z_idx = Self::zero_index(&xAxByCyD);
        assert!(z_idx == 3); // TODO: I think this is always three, so we can clean code up.

        // Compute intermediate values for codomain
        let t1 = zAtBzYtD[1 ^ z_idx];
        let t2 = xAxByCyD[2 ^ z_idx];
        let t3 = zAtBzYtD[3 ^ z_idx];
        let t4 = xAxByCyD[3 ^ z_idx];

        // Invert all four values for codomain and images
        let mut inverse = [t1, t2, t3, t4];
        Fq::batch_invert(&mut inverse);

        // Codomain coefficients
        let mut ABCD = [Fq::ZERO; 4];
        ABCD[z_idx] = Fq::ZERO;
        ABCD[1 ^ z_idx] = t1 * inverse[2];
        ABCD[2 ^ z_idx] = t2 * inverse[3];
        ABCD[3 ^ z_idx] = Fq::ONE;

        // Used for the image computation
        let a_inverse = t3 * inverse[0];
        let b_inverse = t4 * inverse[1];

        let (A, B, C, D) = to_hadamard(&ABCD[0], &ABCD[1], &ABCD[2], &ABCD[3]);
        let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

        (codomain, (a_inverse, b_inverse), z_idx)
    }

    /// Given a point and it's shifted value, compute the image of this
    /// point using the Fq elements precomputed during the codomain
    /// computation
    /// Note: The shift value is needed as one element of the dual coordinates
    /// is zero, so we take values from both and use linear algebra to recover
    /// the correct image.
    ///
    /// Cost: 8S + 10M + 1I
    fn gluing_image(
        T: &ThetaPoint<Fq>,
        T_shift: &ThetaPoint<Fq>,
        a_inv: &Fq,
        b_inv: &Fq,
        z_idx: usize,
    ) -> ThetaPoint<Fq> {
        // Find dual coordinates of point to push through
        let AxByCzDt: [Fq; 4] = T.squared_theta().into();

        // We are in the case where at most one of A, B, C, D is
        // zero, so we need to account for this
        // To recover values, we use the translated point to get
        let AyBxCtDz: [Fq; 4] = T_shift.squared_theta().into();

        // We can always directly compute three elements
        let y = AxByCzDt[1 ^ z_idx] * (*a_inv);
        let z = AxByCzDt[2 ^ z_idx] * (*b_inv);
        let t = AxByCzDt[3 ^ z_idx];

        // To compute the `x` value, we need to compute a scalar, lambda,
        // which we use to normalise the given `xb`. We can always do this,
        // but we are required to compute an inverse which means checking
        // whether z is zero. If z is zero, we can compute lambda by simply
        // extracting the scaled zb and dividing by b from above. However,
        // when z is zero, we instead have to use t. To ensure that this is
        // constant time, we compute both using that inverting zero just
        // gives zero and conditionally swapping lanbda with lambda_t
        let zb = AyBxCtDz[3 ^ z_idx];
        let tb = AyBxCtDz[2 ^ z_idx] * (*b_inv);

        let mut inverse = [zb, tb];
        Fq::batch_invert(&mut inverse);

        // Potentially one of these inverses are zero, but we do both
        // to avoid branching.
        let mut lam = z * inverse[0];
        let lam_t = t * inverse[1];
        lam.set_cond(&lam_t, z.is_zero());

        // Finally we recover x
        let xb = AyBxCtDz[1 ^ z_idx] * (*a_inv);
        let x = xb * lam;

        // We now have values for `x,y,z,t` but to order them we need to use
        // the xor trick as above, so we pack them into an array with the
        // right ordering and then extract them back out
        let mut xyzt = [Fq::ZERO; 4];
        xyzt[z_idx] = x;
        xyzt[1 ^ z_idx] = y;
        xyzt[2 ^ z_idx] = z;
        xyzt[3 ^ z_idx] = t;

        let (x, y, z, t) = to_hadamard(&xyzt[0], &xyzt[1], &xyzt[2], &xyzt[3]);

        ThetaPoint::from_coords(&x, &y, &z, &t)
    }

    /// Compute the gluing (2,2)-isogeny from a ThetaStructure computed
    /// from an elliptic product.
    pub fn gluing_isogeny(
        self,
        P1P2_8: &ProductPoint<Fq>,
        Q1Q2_8: &ProductPoint<Fq>,
        image_points: &[ProductPoint<Fq>],
    ) -> (ThetaStructure<Fq>, Vec<ThetaPoint<Fq>>) {
        // First recover the four torsion below the 8 torsion
        let P1P2_4 = self.double(P1P2_8);
        let Q1Q2_4 = self.double(Q1Q2_8);

        // Use the four torsion to deterministically find basis change
        let M = self.get_base_matrix(&P1P2_4, &Q1Q2_4);

        // Take the points P1, P2 in E1 x E2 and represent them
        // as a theta point on a level 2 structure and map them
        // through the above basis change
        let T1_8 = Self::base_change_couple_point(P1P2_8, M);
        let T2_8 = Self::base_change_couple_point(Q1Q2_8, M);

        // Now it's time to compute the codomain and image of the isogeny
        // with kernel below T1, and T2.
        // We save the zero index, as we use it for the images, and we also
        // can precompute a few inverses to save time for evaluation.
        let (codomain, (a_inv, b_inv), z_idx) = Self::gluing_codomain(&T1_8, &T2_8);

        // We now want to push through a set of points by evaluating each of them
        // under the action of this isogeny. As the domain is an elliptic product,
        // with elements of type CouplePoint, and the codomain is a ThetaStructure
        // with elements of type ThetaPoint, we need a new vector here and as we
        // iteratate through each CouplePoint, we can compute its image and push it
        // to the new vector.
        let mut theta_images: Vec<ThetaPoint<Fq>> = Vec::new();

        // Per image cost =
        // 2 * (16M + 5S) for the CouplePoint addition
        // 2 * 20M for the base change
        // 8S + 4M + 1I for the gluing image
        // Total:
        // 76M + 18S + 1I per point
        for P in image_points.iter() {
            // Need affine coordinates here to do an add, if we didn't we
            // could use faster x-only... Something to think about but no
            // obvious solution.
            let P_sum_T = self.add(P, &P1P2_4);

            // After we have added the points, we can use the gluing formula
            // to recover theta points on the level 2 theta structure. First we
            // must compute the basis change as we did for the kernel:
            let T = Self::base_change_couple_point(P, M);
            let T_shift = Self::base_change_couple_point(&P_sum_T, M);

            // With a point and the shift value from the kernel, we can find
            // the image
            let T_image = Self::gluing_image(&T, &T_shift, &a_inv, &b_inv, z_idx);
            theta_images.push(T_image);
        }

        (codomain, theta_images)
    }
}
