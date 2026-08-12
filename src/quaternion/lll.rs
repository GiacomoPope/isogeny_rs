use crate::bigint::BigIntAlg;
use crate::bigint::rat::BigRat;
use crate::quaternion::algebra::{QuatConfig, RatQuat};
use crate::quaternion::ideal::QuatLeftIdeal;
use crate::quaternion::lattice::IntLattice;

// ========================================================================
// Core LLL Algorithm
// ========================================================================

fn compute_gs<T: BigIntAlg>(
    gram: &[[T; 4]; 4],
) -> ([[BigRat<T>; 4]; 4], [BigRat<T>; 4]) {
    let mut mu: [[BigRat<T>; 4]; 4] = [
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
        [BigRat::zero(), BigRat::zero(), BigRat::zero(), BigRat::zero()],
    ];
    let mut b_star: [BigRat<T>; 4] = [
        BigRat::zero(),
        BigRat::zero(),
        BigRat::zero(),
        BigRat::zero(),
    ];

    for i in 0..4 {
        mu[i][i] = BigRat::one();
        let mut bi_star_norm = BigRat::from_integer(gram[i][i].clone());
        for j in 0..i {
            let mut sum = BigRat::from_integer(gram[i][j].clone());
            for k in 0..j {
                sum = sum - mu[i][k].clone() * mu[j][k].clone() * b_star[k].clone();
            }
            mu[i][j] = sum.clone() / b_star[j].clone();
            bi_star_norm = bi_star_norm - mu[i][j].clone() * mu[i][j].clone() * b_star[j].clone();
        }
        b_star[i] = bi_star_norm;
    }
    (mu, b_star)
}

/// Executes the core LLL reduction on an IntLattice.
pub fn quat_lll_core<T: BigIntAlg, P: QuatConfig<T>>(
    gram: &mut [[T; 4]; 4],
    basis: &mut IntLattice<T, P>,
) {
    let delta = BigRat::new(T::from_i32(99), T::from_i32(100)); // DELTABAR = 0.99

    let mut k = 1;
    while k < 4 {
        let (mut mu, mut b_star) = compute_gs(gram);

        // Size reduction
        let mut reduced = false;
        for i in (0..k).rev() {
            let q = mu[k][i].round();
            if !q.is_zero() {
                // Update basis: b_k = b_k - q * b_i
                for c in 0..4 {
                    basis.generators[k].coords[c] =
                        basis.generators[k].coords[c].clone() - q.clone() * basis.generators[i].coords[c].clone();
                }

                let q_sq = q.clone() * q.clone();
                let two_q = q.clone() + q.clone();

                // Snapshot row k to avoid polluting geometric updates
                let mut old_gram_k = [T::zero(), T::zero(), T::zero(), T::zero()];
                for j in 0..4 {
                    old_gram_k[j] = gram[k][j].clone();
                }

                // Incremental Gram matrix update
                gram[k][k] = old_gram_k[k].clone() - two_q * old_gram_k[i].clone() + q_sq * gram[i][i].clone();

                for j in 0..4 {
                    if j != k {
                        gram[k][j] = old_gram_k[j].clone() - q.clone() * gram[i][j].clone();
                        gram[j][k] = gram[k][j].clone();
                    }
                }

                // Incrementally update mu for subsequent size reductions in the same loop
                for j in 0..i {
                    mu[k][j] = mu[k][j].clone() - BigRat::from_integer(q.clone()) * mu[i][j].clone();
                }

                reduced = true;
            }
        }

        if reduced {
            // Recompute GSO to guarantee exact mu and b_star values for the Lovasz condition
            let gs = compute_gs(gram);
            mu = gs.0;
            b_star = gs.1;
        }

        // Lovasz condition check
        let mu_k_km1 = &mu[k][k - 1];
        let lovasz_lhs = &b_star[k];

        let mu_sq = mu_k_km1.clone() * mu_k_km1.clone();
        let factor = delta.clone() - mu_sq;
        let lovasz_rhs = factor * b_star[k - 1].clone();

        if lovasz_lhs < &lovasz_rhs {
            basis.generators.swap(k, k - 1);

            for i in 0..4 {
                let temp = gram[k][i].clone();
                gram[k][i] = gram[k - 1][i].clone();
                gram[k - 1][i] = temp;
            }
            for i in 0..4 {
                let temp = gram[i][k].clone();
                gram[i][k] = gram[i][k - 1].clone();
                gram[i][k - 1] = temp;
            }

            k = std::cmp::max(1, k - 1);
        } else {
            k += 1;
        }
    }
}

// ========================================================================
// Qlapoti Support Applications
// ========================================================================

/// Computes the Gram matrix for a left ideal class.
pub fn quat_lideal_class_gram<T: BigIntAlg, P: QuatConfig<T>>(
    lideal: &QuatLeftIdeal<T, P>,
) -> [[T; 4]; 4] {
    let mut gram = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];
    let p = P::p();

    let mut divisor = lideal.lattice.denom.clone() * lideal.lattice.denom.clone();
    divisor = divisor * lideal.norm.clone();

    for i in 0..4 {
        for j in 0..=i {
            let gi = &lideal.lattice.basis.generators[i].coords;
            let gj = &lideal.lattice.basis.generators[j].coords;

            let mut sum = gi[0].clone() * gj[0].clone();
            sum = sum + gi[1].clone() * gj[1].clone();
            sum = sum + p.clone() * gi[2].clone() * gj[2].clone();
            sum = sum + p.clone() * gi[3].clone() * gj[3].clone();

        let trace = sum.clone() + sum;
            gram[i][j] = trace.clone() / divisor.clone();
            gram[j][i] = gram[i][j].clone();
        }
    }
    gram
}

/// L2-reduces the basis of the left ideal using the exact Trace form.
pub fn quat_lideal_reduce_basis<T: BigIntAlg, P: QuatConfig<T>>(
    reduced: &mut IntLattice<T, P>,
    gram: &mut [[T; 4]; 4],
    lideal: &QuatLeftIdeal<T, P>,
) {
    let mut exact_gram = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];
    let p = P::p();

    // Construct exact integer Gram matrix to prevent LLL truncation
    for i in 0..4 {
        for j in 0..=i {
            let gi = &lideal.lattice.basis.generators[i].coords;
            let gj = &lideal.lattice.basis.generators[j].coords;

            let mut sum = gi[0].clone() * gj[0].clone();
            sum = sum + gi[1].clone() * gj[1].clone();
            sum = sum + p.clone() * gi[2].clone() * gj[2].clone();
            sum = sum + p.clone() * gi[3].clone() * gj[3].clone();

            let trace = sum.clone() + sum;
            exact_gram[i][j] = trace.clone();
            exact_gram[j][i] = trace;
        }
    }

    *reduced = lideal.lattice.basis.clone();
    quat_lll_core(&mut exact_gram, reduced);

    // Recompute exact gram of the reduced basis
    for i in 0..4 {
        for j in 0..=i {
            let gi = &reduced.generators[i].coords;
            let gj = &reduced.generators[j].coords;

            let mut sum = gi[0].clone() * gj[0].clone();
            sum = sum + gi[1].clone() * gj[1].clone();
            sum = sum + p.clone() * gi[2].clone() * gj[2].clone();
            sum = sum + p.clone() * gi[3].clone() * gj[3].clone();

            let trace = sum.clone() + sum;
            exact_gram[i][j] = trace.clone();
            exact_gram[j][i] = trace;
        }
    }

    let norm = lideal.norm.clone();
    let two = T::from_i32(2);

    // return the normalized Class Gram Matrix
    for i in 0..4 {
        for j in 0..=i {
            gram[i][j] = exact_gram[i][j].clone() / norm.clone();
        }
    }

    for i in 0..4 {
        gram[i][i] = gram[i][i].clone() / two.clone();
        for j in (i + 1)..4 {
            gram[i][j] = T::zero();
        }
    }
}

/// Extracts the shortest equivalent ideal
pub fn quat_lideal_shortest_equivalent<T: BigIntAlg, P: QuatConfig<T>>(
    lideal: &QuatLeftIdeal<T, P>,
) -> (QuatLeftIdeal<T, P>, RatQuat<T, P>) {
    let mut red = IntLattice::zero();
    let mut gram = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];

    quat_lideal_reduce_basis(&mut red, &mut gram, lideal);

    let mut new_alpha = RatQuat::new(red.generators[0].clone(), lideal.lattice.denom.clone());
    let elem = new_alpha.clone();

    new_alpha.denom = new_alpha.denom * lideal.norm.clone();
    new_alpha.num.coords[0] = -new_alpha.num.coords[0].clone();

    let equiv = QuatLeftIdeal::mul(lideal, &new_alpha);

    (equiv, elem)
}

/// Multiplies two left ideals and returns the LLL reduced result and Gram matrix.
pub fn quat_lideal_lideal_mul_reduced<T: BigIntAlg, P: QuatConfig<T>>(
    prod: &mut QuatLeftIdeal<T, P>,
    gram: &mut [[T; 4]; 4],
    lideal1: &QuatLeftIdeal<T, P>,
    lideal2: &QuatLeftIdeal<T, P>,
) {
    *prod = QuatLeftIdeal::mul_ideal(lideal1, lideal2);

    let mut red_basis = crate::quaternion::lattice::IntLattice::zero();
    quat_lideal_reduce_basis(&mut red_basis, gram, prod);
    prod.lattice.basis = red_basis;
}
