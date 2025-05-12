use fp2::fq::Fq as FqTrait;

use super::{curve::Curve, point::PointX};

/// Compute [2]P in place using projective (A + 2) / 4 = (A24 : C24)
/// Cost: 2S + 4M
#[inline(always)]
fn xdbl_proj<Fq: FqTrait>(A24: &Fq, C24: &Fq, P: &mut PointX<Fq>) {
    let mut t0 = P.X + P.Z;
    t0.set_square();
    let mut t1 = P.X - P.Z;
    t1.set_square();
    let t2 = t0 - t1;
    t1 *= *C24;
    P.X = t0 * t1;
    t0 = t2 * (*A24);
    t0 = t0 + t1;
    P.Z = t0 * t2;
}

/// Compute [2^n]P in place using projective (A + 2) / 4 = (A24 : C24)
/// Cost: n * (2S + 4M)
fn xdbl_proj_iter<Fq: FqTrait>(A24: &Fq, C24: &Fq, P: &mut PointX<Fq>, n: usize) {
    for _ in 0..n {
        xdbl_proj(A24, C24, P);
    }
}

/// Compute the codomain of the 2-isogeny E -> E/<ker> for ker != (0 : 1)
fn two_isogeny_codomain<Fq: FqTrait>(ker: &PointX<Fq>) -> (Fq, Fq) {
    let mut A24 = ker.X.square();
    let C24 = ker.Z.square();
    A24 = C24 - A24;
    (A24, C24)
}

/// Evaluate a point Q in place under the action of the 2-isogeny E -> E/<ker>
/// for ker != (0 : 1)
fn two_isogeny_eval<Fq: FqTrait>(ker: &PointX<Fq>, Q: &mut PointX<Fq>) {
    let mut t0 = ker.X + ker.Z;
    let mut t1 = ker.X - ker.Z;
    let mut t2 = Q.X + Q.Z;
    let mut t3 = Q.X - Q.Z;

    t0 *= t3;
    t1 *= t2;
    t2 = t0 + t1;
    t3 = t0 - t1;

    Q.X *= t2;
    Q.Z *= t3;
}

/// Compute the codomain of the 2-isogeny E -> E/<ker> for ker == (0 : 1)
fn two_isogeny_codomain_singular<Fq: FqTrait>(A24: &mut Fq, C24: &mut Fq) -> (Fq, Fq) {
    let mut t0 = A24.mul2();
    t0 -= *C24;
    t0.set_mul2();
    t0 /= *C24;
    let c0 = t0;
    *A24 = t0.mul2();
    t0.set_square();
    t0 -= Fq::FOUR;
    let r = t0.set_sqrt();
    assert!(r == u32::MAX); // TODO
    let c1 = -t0;
    *C24 = t0.mul2();
    *A24 += *C24;
    C24.set_mul2();

    (c0, c1)
}

/// Evaluate a point Q in place under the action of the 2-isogeny E -> E/<ker>
/// for ker = (0 : 1)
fn two_isogeny_eval_singular<Fq: FqTrait>(c0: &Fq, c1: &Fq, Q: &mut PointX<Fq>) {
    let t0 = Q.X * Q.Z;
    let mut t1 = (*c0) * Q.Z;
    t1 += Q.X;
    t1 *= Q.X;

    Q.X = Q.Z.square();
    Q.X += t1;
    Q.Z = t0 * (*c1);
}

/// Compute a 2^n isogeny using the naive approach
pub fn two_isogeny_chain_naive<Fq: FqTrait>(
    domain: &Curve<Fq>,
    kernel: &PointX<Fq>,
    n: usize,
    images: &mut [PointX<Fq>],
) -> Curve<Fq> {
    let mut A24 = domain.A24;
    let mut C24 = Fq::ONE;

    let mut ker: PointX<Fq> = *kernel;
    let mut ker_step: PointX<Fq>;

    for i in 0..n {
        // Double the kernel to get a point of order 2
        ker_step = ker;
        xdbl_proj_iter(&A24, &C24, &mut ker_step, n - i - 1);

        if ker_step.X.is_zero() == u32::MAX {
            // Compute the codomain from ker_step for kernel (0 : 1)
            let (c0, c1) = two_isogeny_codomain_singular(&mut A24, &mut C24);

            // Push through the kernel
            two_isogeny_eval_singular(&c0, &c1, &mut ker);

            // Push through the points to evaluate
            for P in images.iter_mut() {
                two_isogeny_eval_singular(&c0, &c1, P);
            }
        } else {
            // Compute the codomain from ker_step
            (A24, C24) = two_isogeny_codomain(&ker_step);

            // Push through the kernel
            two_isogeny_eval(&ker_step, &mut ker);

            // Push through the points to evaluate
            for P in images.iter_mut() {
                two_isogeny_eval(&ker_step, P);
            }
        }
    }

    // Compute A from (A24 : C24)
    let mut A = A24 + A24;
    A -= C24;
    A += A;
    A /= C24;

    Curve::new(&A)
}
