use fp2::fq::Fq as FqTrait;

use super::theta_gluing::gluing_isogeny;
use super::theta_isogeny::{two_isogeny, two_isogeny_to_product};
use super::theta_point::ThetaPoint;
use super::theta_splitting::{split_to_product, splitting_isomorphism};
use crate::elliptic::product::{CouplePoint, EllipticProduct};

// ========================================================
// Compute the isogeny between elliptic products
// ========================================================

/// Compute an isogeny between elliptic products, naive method with no
/// optimised strategy. Only here for benchmarking
pub fn product_isogeny_no_strategy<Fq: FqTrait>(
    E1E2: &EllipticProduct<Fq>,
    P1P2: &CouplePoint<Fq>,
    Q1Q2: &CouplePoint<Fq>,
    image_points: &[CouplePoint<Fq>],
    n: usize,
) -> (EllipticProduct<Fq>, Vec<CouplePoint<Fq>>) {
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
    let mut Tp1: ThetaPoint<Fq>;
    let mut Tp2: ThetaPoint<Fq>;
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
