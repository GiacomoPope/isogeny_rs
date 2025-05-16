use fp2::fq::Fq as FqTrait;

use super::elliptic_product::{EllipticProduct, ProductPoint};
use super::theta_isogeny::{two_isogeny, two_isogeny_to_product};
use super::theta_point::ThetaPoint;
use super::theta_splitting::{split_to_product, splitting_isomorphism};

// ========================================================
// Compute the isogeny between elliptic products
// ========================================================

/// Compute an isogeny between elliptic products, naive method with no
/// optimised strategy. Only here for benchmarking

impl<Fq: FqTrait> EllipticProduct<Fq> {
    /// Compute an 2^n isogeny (E1 x E2 -> E3 x E4) with kernel <(P1, P2), (Q1, Q2)>
    /// Assumes points Pi, Qi have order 2^(n + 2) to allow for fast computation of
    /// the codomain for the last two steps without requiring square roots.
    /// Returns the codomain and evaluated `ProductPoint` through the isogeny together
    /// with a u32 which is `0xFF..FF` on success or `0x00..00` on failure.
    ///
    /// Naive implementation which doubles all the way to obtain the points of
    /// order 8 above the kernel. Much slower than the balanced strategy.
    pub fn elliptic_product_isogeny_no_strategy(
        self,
        P1P2: &ProductPoint<Fq>,
        Q1Q2: &ProductPoint<Fq>,
        n: usize,
        image_points: &[ProductPoint<Fq>],
    ) -> (EllipticProduct<Fq>, Vec<ProductPoint<Fq>>, u32) {
        // Values to track whether the isogeny has been computed successfully
        let mut ok = u32::MAX;
        let mut check;

        // Store the number of image points we wish to evaluate to
        // ensure we return them all from the points we push through
        let num_image_points = image_points.len();

        // We store the points we want an image for at the end of the
        // kernel points we evaluate.
        let mut kernel_couple_pts = vec![*P1P2, *Q1Q2];
        kernel_couple_pts.extend_from_slice(image_points);

        // Compute points of order 8
        let P1P2_8 = self.double_iter(&P1P2, n - 1);
        let Q1Q2_8 = self.double_iter(&Q1Q2, n - 1);

        // Compute Gluing isogeny
        let (mut domain, mut kernel_pts) =
            self.gluing_isogeny(&P1P2_8, &Q1Q2_8, &kernel_couple_pts);

        // Do all remaining steps
        let mut Tp1: ThetaPoint<Fq>;
        let mut Tp2: ThetaPoint<Fq>;
        for k in 1..n {
            // Repeatedly double to obtain points in the 8-torsion below the kernel
            Tp1 = domain.double_iter(&kernel_pts[0], n - k - 1);
            Tp2 = domain.double_iter(&kernel_pts[1], n - k - 1);

            // For the last two steps, we need to be careful because of the zero-null
            // coordinates appearing from the product structure. To avoid these, we
            // use the hadamard transform to avoid them,
            (domain, check) = if k == (n - 2) {
                two_isogeny(&Tp1, &Tp2, &mut kernel_pts, [false, false])
            } else if k == (n - 1) {
                two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
            } else {
                two_isogeny(&Tp1, &Tp2, &mut kernel_pts, [false, true])
            };
            ok &= check;
        }

        // Use a symplectic transform to first get the domain into a compatible form
        // for splitting
        let check: u32;
        (domain, check) = splitting_isomorphism(domain, &mut kernel_pts);
        ok &= check;

        // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
        // onto this product
        let img_points = &kernel_pts[kernel_pts.len() - num_image_points..];
        let (product, couple_points) = split_to_product(&domain, img_points);

        (product, couple_points, ok)
    }

    /// Compute an 2^n isogeny (E1 x E2 -> E3 x E4) with kernel <(P1, P2), (Q1, Q2)>
    /// using a balanced strategy.
    /// Assumes points Pi, Qi have order 2^(n + 2) to allow for fast computation of
    /// the codomain for the last two steps without requiring square roots.
    /// Returns the codomain and evaluated `ProductPoint` through the isogeny together
    /// with a u32 which is `0xFF..FF` on success or `0x00..00` on failure.
    ///
    /// Naive implementation which doubles all the way to obtain the points of
    /// order 8 above the kernel. Much slower than the balanced strategy.
    pub fn elliptic_product_isogeny(
        self,
        P1P2: &ProductPoint<Fq>,
        Q1Q2: &ProductPoint<Fq>,
        n: usize,
        image_points: &[ProductPoint<Fq>],
    ) -> (EllipticProduct<Fq>, Vec<ProductPoint<Fq>>, u32) {
        // Store the number of image points we wish to evaluate to
        // ensure we return them all from the points we push through
        let num_image_points = image_points.len();

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - n.leading_zeros() + 1) as usize;

        // Store points of order 2^i for the balanced strategy on E1
        let mut product_strategy_pts: Vec<ProductPoint<Fq>> =
            vec![ProductPoint::INFINITY; 2 * space];

        // The values i such that each point in stategy_points has order 2^i
        let mut orders: Vec<usize> = vec![0; space];

        // Initalise the first values for the strategy
        product_strategy_pts[0] = *P1P2;
        product_strategy_pts[1] = *Q1Q2;
        orders[0] = n;

        // Value to determine success / failure of isogeny chain
        let mut ok = u32::MAX;
        let mut check: u32;

        // Step One: Perform doubling on the ProductPoints and compute the
        // codomain from gluing. Keep intermediate doubles to push through
        // the isogeny to save on doublings later.
        let mut k = 0;
        while orders[k] != 1 {
            k += 1;
            // As gluing is expensive we offset the balanced isogeny
            // TODO: optimise choice of bound?
            let m = if orders[k - 1] >= 16 {
                orders[k - 1] >> 1
            } else {
                orders[k - 1] - 1
            };

            // Double the points m times.
            // Points are filled in pairs [2^m] P1P2 and  [2^m] Q1Q2
            product_strategy_pts[2 * k] = self.double_iter(&product_strategy_pts[2 * k - 2], m);
            product_strategy_pts[2 * k + 1] = self.double_iter(&product_strategy_pts[2 * k - 1], m);
            orders[k] = orders[k - 1].saturating_sub(m);
        }

        // We append the points we want the image for to the end of strategy points.
        // We'll extract them from the end later.
        product_strategy_pts.extend_from_slice(image_points);

        // Compute the Gluing isogeny and push through stategy_points to get a new
        // vector of ThetaPoints of the same length.
        let (mut domain, mut strategy_pts) = self.gluing_isogeny(
            &product_strategy_pts[2 * k],
            &product_strategy_pts[2 * k + 1],
            &product_strategy_pts,
        );

        // Reduce the order of the points we evaluated
        for i in 0..k {
            orders[i] -= 1;
        }

        // Step One: Perform doubling on the ProductPoints and compute the
        // codomain from gluing. Keep intermediate doubles to push through
        // the isogeny to save on doublings later.
        k -= 1;
        for i in 1..n {
            while orders[k] != 1 {
                k += 1;
                let m = orders[k - 1] >> 1;

                // Double the theta points m times.
                strategy_pts[2 * k] = domain.double_iter(&strategy_pts[2 * k - 2], m);
                strategy_pts[2 * k + 1] = domain.double_iter(&strategy_pts[2 * k - 1], m);
                orders[k] = orders[k - 1].saturating_sub(m);
            }

            // Extract out the kernel for this step.
            let T1 = strategy_pts[2 * k];
            let T2 = strategy_pts[2 * k + 1];

            // For the last two steps, we need to be careful because of the zero-null
            // coordinates appearing from the product structure. To avoid these, we
            // use the hadamard transform to avoid them,
            (domain, check) = if i == (n - 2) {
                two_isogeny(&T1, &T2, &mut strategy_pts, [false, false])
            } else if i == (n - 1) {
                two_isogeny_to_product(&T1, &T2, &mut strategy_pts)
            } else {
                two_isogeny(&T1, &T2, &mut strategy_pts, [false, true])
            };
            ok &= check;

            // Reduce the order of the points we evaluated
            for i in 0..k {
                orders[i] -= 1;
            }
            k = k.saturating_sub(1);
        }

        // Use a symplectic transform to first get the domain into a compatible form
        // for splitting
        (domain, check) = splitting_isomorphism(domain, &mut strategy_pts);
        ok &= check;

        // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
        // onto this product
        let img_points = &strategy_pts[strategy_pts.len() - num_image_points..];
        let (product, couple_points) = split_to_product(&domain, &img_points);

        (product, couple_points, ok)
    }
}
