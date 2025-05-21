use fp2::fq::Fq as FqTrait;

use crate::theta::theta_point::ThetaPoint;

use super::elliptic_product::{EllipticProduct, ProductPoint};
use super::theta_isogeny::{two_isogeny, two_isogeny_to_product};
use super::theta_splitting::{split_to_product, splitting_isomorphism};

impl<Fq: FqTrait> EllipticProduct<Fq> {
    /// Compute an 2^n isogeny (E1 x E2 -> E3 x E4) with kernel <(P1, P2), (Q1, Q2)>
    /// using a balanced strategy.
    /// Assumes points Pi, Qi have order 2^(n + 2) to allow for fast computation of
    /// the codomain for the last two steps without requiring square roots.
    /// Returns the codomain and evaluated `ProductPoint` through the isogeny together
    /// with a u32 which is `0xFF..FF` on success or `0x00..00` on failure.
    pub fn elliptic_product_isogeny(
        &self,
        P1P2: &ProductPoint<Fq>,
        Q1Q2: &ProductPoint<Fq>,
        len: usize,
        image_points: &[ProductPoint<Fq>],
    ) -> (EllipticProduct<Fq>, Vec<ProductPoint<Fq>>, u32) {
        // We push the image points through at the same time as the strategy
        // points, so we need to know how many images we're computing to keep
        // track of them.
        let n = image_points.len();

        // Compute the amount of space we need for the balanced strategy.
        let space = (usize::BITS - len.leading_zeros() + 1) as usize;

        // Store points of order 2^i for the balanced strategy. We need two
        // vectors here, as the first step computes with elements of type
        // ProductPoint, while every other step computes points of type
        // ThetaPoint.
        let mut product_pts: Vec<ProductPoint<Fq>> = vec![ProductPoint::INFINITY; 2 * space + n];
        let mut theta_pts: Vec<ThetaPoint<Fq>> = vec![ThetaPoint::ZERO; 2 * space + n];

        // The values i such that each point in stategy_points has order 2^i
        let mut orders: Vec<usize> = vec![0; space];

        // Set the first elements of the vector to the points we want to push
        // through the isogeny.
        product_pts[..n].copy_from_slice(image_points);

        // Then add the kernel points afterwards.
        product_pts[n] = *P1P2;
        product_pts[n + 1] = *Q1Q2;

        // Initalise the orders list, points in the above vectors have order
        // 2^(orders[i] + 2), as we use the 8-torsion above.
        orders[0] = len;

        // Value to determine success / failure of isogeny chain
        let mut ok = u32::MAX;
        let mut check: u32;

        // Step One: Perform doubling on the ProductPoints and compute the
        // codomain from gluing. Keep intermediate doubles to push through
        // the isogeny to save on doublings later.
        let mut k = 0;
        while orders[k] != 1 {
            k += 1;
            let m = if orders[k - 1] >= 16 {
                orders[k - 1] >> 1
            } else {
                orders[k - 1] - 1
            };

            // Double the points m times.
            // Points are filled in pairs [2^m] P1P2 and  [2^m] Q1Q2
            product_pts[n + 2 * k] = self.double_iter(&product_pts[n + 2 * k - 2], m);
            product_pts[n + 2 * k + 1] = self.double_iter(&product_pts[n + 2 * k - 1], m);
            orders[k] = orders[k - 1].saturating_sub(m);
        }

        // Compute the Gluing isogeny and push through product_strategy_pts through
        // into the vector of ThetaPoints `strategy_pts`.
        let mut domain = self.gluing_isogeny(
            &product_pts[n + 2 * k],
            &product_pts[n + 2 * k + 1],
            &product_pts[..(2 * k + n)],
            &mut theta_pts,
        );

        // Reduce the order of the points we evaluated
        for ord in orders.iter_mut().take(k) {
            *ord -= 1;
        }
        k -= 1;

        // Step One: Perform doubling on the ProductPoints and compute the
        // codomain from gluing. Keep intermediate doubles to push through
        // the isogeny to save on doublings later.
        for i in 1..len {
            while orders[k] != 1 {
                k += 1;
                let m = orders[k - 1] >> 1;

                // Double the theta points m times.
                theta_pts[2 * k + n] = domain.double_iter(&theta_pts[2 * k + n - 2], m);
                theta_pts[2 * k + n + 1] = domain.double_iter(&theta_pts[2 * k + n - 1], m);
                orders[k] = orders[k - 1].saturating_sub(m);
            }

            // Extract out the kernel for this step.
            let T1 = theta_pts[2 * k + n];
            let T2 = theta_pts[2 * k + n + 1];

            // For the last two steps, we need to be careful because of the zero-null
            // coordinates appearing from the product structure. To avoid these, we
            // use the hadamard transform to avoid them,
            (domain, check) = if i == (len - 2) {
                two_isogeny(&T1, &T2, &mut theta_pts[..(2 * k + n)], [false, false])
            } else if i == (len - 1) {
                two_isogeny_to_product(&T1, &T2, &mut theta_pts[..(2 * k + n)])
            } else {
                two_isogeny(&T1, &T2, &mut theta_pts[..(2 * k + n)], [false, true])
            };
            ok &= check;

            // Reduce the order of the points we evaluated
            for ord in orders.iter_mut().take(k) {
                *ord -= 1;
            }
            k = k.saturating_sub(1);
        }

        // Use a symplectic transform to first get the domain into a compatible form
        // for splitting
        (domain, check) = splitting_isomorphism(domain, &mut theta_pts);
        ok &= check;

        // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
        // onto this product
        let eval_points = &theta_pts[..n];
        let mut couple_points = vec![ProductPoint::INFINITY; n];
        let product = split_to_product(&domain, eval_points, &mut couple_points);

        (product, couple_points, ok)
    }
}
