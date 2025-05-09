// ========================================================
// Main Method! Compute the isogeny between elliptic
// products
// ========================================================

// A general comment about "optimal strategies" -- For the isogeny we
// have a kernel of order 2^n, and to make a step on the (2,2)-isogeny
// graph what we want to do is scale this point to get 2^(n-1) * P,
// which is a point of order two, which then allows us to compute the
// codomain and images. However, as doubling is more expensive than
// computing an image (in the general case) it's worth storing many
// values along the way while doubling and pushing them all through the
// isogeny. As these pushed through points have the same order, the
// subsequent steps will need to be doubled less (at the cost of needing
// to push through more points.)
//
// For a proper reference, see Sec 4.2 of
// https://eprint.iacr.org/2011/506.pdf
//
// Gluing: Doubling cost: 16M 16S (We have to double two elliptic curve
// points) Image cost: 76M + 18S + 1I (Very expensive!)
//
// All other steps: Doubling cost: 8S + 6M Image cost: 4S + 3M
//
// So ignoring the gluing step, we see images have 1/2 the cost
// (mathematically this is expected as our doubling formula is
// essentially just two isogeny images) and so the optimised strategy
// is computed with a weight that doubling is 2X images.
//
// For a function to see how strategies are computed, see strategy.py
// The current implementation "knows" that the gluing is more expensive
// and so has extra costs for the leftmost branch of the tree.

/// Compute an isogeny between elliptic products, naive method with no
/// optimised strategy. Only here for benchmarking
pub fn product_isogeny_no_strategy(
    E1E2: &EllipticProduct,
    P1P2: &CouplePoint,
    Q1Q2: &CouplePoint,
    image_points: &[CouplePoint],
    n: usize,
) -> (EllipticProduct, Vec<CouplePoint>) {
    // Store the number of image points we wish to evaluate to
    // ensure we return them all from the points we push through
    let num_image_points = image_points.len();

    // Convert the &[...] to a vector so we can add points to this
    // dynamically during the optimal strategy
    let mut kernel_couple_pts = image_points.to_vec();

    // Include the kernel inside the vector of points
    // to evaluate. At each step, every element of the
    // vector should be evaluated
    kernel_couple_pts.push(*P1P2);
    kernel_couple_pts.push(*Q1Q2);

    // Compute points of order 8
    let P1P2_8 = E1E2.double_iter(&P1P2, n - 1);
    let Q1Q2_8 = E1E2.double_iter(&Q1Q2, n - 1);

    // Compute Gluing isogeny
    let (mut domain, mut kernel_pts) = gluing_isogeny(&E1E2, &P1P2_8, &Q1Q2_8, &kernel_couple_pts);

    // Do all remaining steps
    let mut Tp1: ThetaPoint;
    let mut Tp2: ThetaPoint;
    for k in 1..n {
        // Repeatedly double to obtain points in the 8-torsion below the kernel
        Tp1 = domain.double_iter(&kernel_pts[num_image_points], n - k - 1);
        Tp2 = domain.double_iter(&kernel_pts[num_image_points + 1], n - k - 1);

        // For the last two steps, we need to be careful because of the zero-null
        // coordinates appearing from the product structure. To avoid these, we
        // use the hadamard transform to avoid them,
        if k == (n - 2) {
            domain = two_isogeny(&Tp1, &Tp2, &mut kernel_pts, [false, false])
        } else if k == (n - 1) {
            domain = two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
        } else {
            domain = two_isogeny(&Tp1, &Tp2, &mut kernel_pts, [false, true])
        }
    }

    // Use a symplectic transform to first get the domain into a compatible form
    // for splitting
    domain = splitting_isomorphism(domain, &mut kernel_pts);

    // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
    // onto this product
    let (product, couple_points) = split_to_product(&domain, &kernel_pts, num_image_points);

    (product, couple_points)
}

/// Compute an isogeny between elliptic products, use an optimised
/// strategy for all steps assuming doubling is always more expensive
/// that images, which is not true for gluing.
pub fn product_isogeny(
    E1E2: &EllipticProduct,
    P1P2: &CouplePoint,
    Q1Q2: &CouplePoint,
    image_points: &[CouplePoint],
    n: usize,
    strategy: &[usize],
) -> (EllipticProduct, Vec<CouplePoint>) {
    // Store the number of image points we wish to evaluate to
    // ensure we return them all from the points we push through
    let num_image_points = image_points.len();

    // Convert the &[...] to a vector so we can add points to this
    // dynamically during the optimal strategy
    let mut kernel_couple_pts = image_points.to_vec();

    // Include the kernel inside the vector of points
    // to evaluate. At each step, every element of the
    // vector should be evaluated
    kernel_couple_pts.push(*P1P2);
    kernel_couple_pts.push(*Q1Q2);

    // Bookkeeping for optimised strategy
    let mut strat_idx = 0;
    let mut level: Vec<usize> = vec![0];
    let mut prev: usize = level.iter().sum();

    // =======================================================
    // Gluing Step
    // TODO:
    // Because of type differences there's annoying code reuse
    // for the optimal strategy here and again for every step
    // in the chain thereafter. Which is bothersome. Maybe there
    // is a better way to write this...
    // =======================================================
    let mut ker1 = *P1P2;
    let mut ker2 = *Q1Q2;

    while prev != (n - 1) {
        // Add the next strategy to the level
        level.push(strategy[strat_idx]);

        // Double the points according to the strategy
        ker1 = E1E2.double_iter(&ker1, strategy[strat_idx]);
        ker2 = E1E2.double_iter(&ker2, strategy[strat_idx]);

        // Add these points to the image points
        kernel_couple_pts.push(ker1);
        kernel_couple_pts.push(ker2);

        // Update the strategy bookkeepping
        prev += strategy[strat_idx];
        strat_idx += 1;
    }

    // Clear out the used kernel point and update level
    kernel_couple_pts.pop();
    kernel_couple_pts.pop();
    level.pop();

    // Compute Gluing isogeny
    let (mut domain, mut kernel_pts) = gluing_isogeny(&E1E2, &ker1, &ker2, &kernel_couple_pts);

    // ======================================================
    // All other steps
    // Compute the (2^n-1, 2^n-1)-chain in the theta model
    // =======================================================

    let mut Tp1: ThetaPoint;
    let mut Tp2: ThetaPoint;
    let mut kernel_len: usize;

    // Do all remaining steps
    for k in 1..n {
        prev = level.iter().sum();
        kernel_len = kernel_pts.len();

        Tp1 = kernel_pts[kernel_len - 2];
        Tp2 = kernel_pts[kernel_len - 1];

        while prev != (n - 1 - k) {
            // Add the next strategy to the level
            level.push(strategy[strat_idx]);

            // Double the points according to the strategy
            Tp1 = domain.double_iter(&Tp1, strategy[strat_idx]);
            Tp2 = domain.double_iter(&Tp2, strategy[strat_idx]);

            // Add these points to the image points
            kernel_pts.push(Tp1);
            kernel_pts.push(Tp2);

            // Update the strategy bookkeepping
            prev += strategy[strat_idx];
            strat_idx += 1;
        }

        // Clear out the used kernel point and update level
        kernel_pts.pop();
        kernel_pts.pop();
        level.pop();

        // For the last two steps, we need to be careful because of the zero-null
        // coordinates appearing from the product structure. To avoid these, we
        // use the hadamard transform to avoid them,
        if k == (n - 2) {
            domain = two_isogeny(&Tp1, &Tp2, &mut kernel_pts, [false, false])
        } else if k == (n - 1) {
            domain = two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
        } else {
            domain = two_isogeny(&Tp1, &Tp2, &mut kernel_pts, [false, true])
        }
    }

    // Use a symplectic transform to first get the domain into a compatible form
    // for splitting
    domain = splitting_isomorphism(domain, &mut kernel_pts);

    // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
    // onto this product
    let (product, couple_points) = split_to_product(&domain, &kernel_pts, num_image_points);

    (product, couple_points)
}
