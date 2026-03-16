#![allow(non_snake_case)]

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};

const MODULUS: [u64; 4] = [
    0x65396406224d1bc7,
    0x60892bae5986e2f3,
    0x8f2e6d4fb53f8cc7,
    0x0003e368df883c71,
];
fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = MODULUS,);

use isogeny::{
    elliptic::{curve::Curve, point::PointX},
    polynomial_ring::poly::{EvalTree, Polynomial},
    utilities::test_utils::drng::DRNG,
};

type PR = Polynomial<Fp2>;

const P_X_RE_BYTES: [u8; 31] = [
    161, 11, 202, 60, 238, 183, 30, 198, 113, 154, 4, 46, 251, 33, 116, 68, 50, 187, 152, 163, 215,
    172, 206, 196, 15, 96, 128, 81, 141, 27, 2,
];
const P_X_IM_BYTES: [u8; 31] = [
    244, 36, 105, 80, 212, 91, 28, 64, 95, 113, 78, 111, 165, 12, 140, 76, 43, 91, 43, 245, 47, 7,
    52, 200, 71, 10, 182, 128, 100, 211, 1,
];
const PX_163: Fp2 = Fp2::const_decode_no_check(&P_X_RE_BYTES, &P_X_IM_BYTES);

fn partition(degree: usize) -> (usize, usize, usize) {
    let size_j = (((degree - 1) as f64).sqrt() as usize) / 2;
    let size_i = (degree - 1) / (4 * size_j);
    let size_k = (degree - 4 * size_j * size_i - 1) / 2;
    (size_j, size_i, size_k)
}

fn random_fps(rng: &mut DRNG, n: usize) -> Vec<Fp2> {
    (0..n).map(|_| Fp2::rand(rng)).collect()
}

fn random_quadratic_leaves_palindrome(rng: &mut DRNG, n: usize) -> Vec<[Fp2; 3]> {
    (0..n)
        .map(|_| {
            let a = Fp2::rand(rng);
            let b = Fp2::rand(rng);
            [a, b, a] // palindrome, as in the codomain path
        })
        .collect()
}

fn random_quadratic_leaves_general(rng: &mut DRNG, n: usize) -> Vec<[Fp2; 3]> {
    (0..n)
        .map(|_| [Fp2::rand(rng), Fp2::rand(rng), Fp2::rand(rng)])
        .collect()
}

// ===============================================================================
// Stage 1: product tree from quadratic leaves
//
// Called twice for the codomain (E0J / E1J, palindrome leaves) and once per
// image point (general leaves).  sJ leaves -> degree-2·sJ polynomial.
// ===============================================================================
fn bench_product_tree(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage1_product_tree");
    let mut rng = DRNG::from_seed(b"product_tree_bench");

    for &ell in &[103usize, 163, 211] {
        let (size_j, size_i, size_k) = partition(ell);
        eprintln!(
            "ell={ell}: sJ={size_j}  sI={size_i}  sK={size_k}  \
             check=4·sI·sJ+2·sK+1={}",
            4 * size_i * size_j + 2 * size_k + 1
        );

        // Codomain path: palindrome leaves [a, b, a]
        let pal = random_quadratic_leaves_palindrome(&mut rng, size_j);
        group.bench_with_input(
            BenchmarkId::new("palindrome_codomain", ell),
            &pal,
            |b, leaves| b.iter(|| PR::root_from_quadratic_leaves(black_box(leaves))),
        );

        // Image-eval path: general leaves [c0, c1, c2]
        let g = random_quadratic_leaves_general(&mut rng, size_j);
        group.bench_with_input(
            BenchmarkId::new("general_image_eval", ell),
            &g,
            |b, leaves| b.iter(|| PR::root_from_quadratic_leaves(black_box(leaves))),
        );
    }
    group.finish();
}

// ===============================================================================
// Stage 2: EvalTree construction over the I-partition roots
//
// Paid once per isogeny call, not per image point.
// ===============================================================================
fn bench_eval_tree_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage2_eval_tree_construction");
    let mut rng = DRNG::from_seed(b"eval_tree_bench");

    for &ell in &[103usize, 163, 211] {
        let (_, size_i, _) = partition(ell);
        let roots = random_fps(&mut rng, size_i);

        group.bench_with_input(BenchmarkId::new("build", ell), &roots, |b, r| {
            b.iter(|| EvalTree::new(black_box(r)))
        });
    }
    group.finish();
}

// ===============================================================================
// Stage 3: Resultant: remainder-tree vs Horner
//
// NOTE: sqrt_velu computes TWO resultants per isogeny (E0J and E1J), plus two
// more per image point.
// ===============================================================================
fn bench_resultant(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage3_resultant");
    let mut rng = DRNG::from_seed(b"resultant_bench");

    for &ell in &[103usize, 163, 211] {
        let (size_j, size_i, _) = partition(ell);

        let ej_poly = PR::rand(&mut rng, 2 * size_j + 1);
        let roots = random_fps(&mut rng, size_i);
        let eval_tree = EvalTree::new(&roots);

        let n = usize::BITS - roots.len().leading_zeros() - 1;
        let two_n = 1 << n;
        let balanced_roots = &roots[..two_n];
        let rem_roots = &roots[two_n..];
        let balanced_tree = EvalTree::new(&balanced_roots);

        // (a) prebuilt eval tree
        group.bench_with_input(
            BenchmarkId::new("remainder_tree_prebuilt", ell),
            &(&ej_poly, &eval_tree),
            |b, (poly, tree)| {
                b.iter(|| black_box(poly).resultant_from_roots_with_tree(black_box(tree)))
            },
        );

        // (b) Horner O(sI · deg) field multiplications, no allocations
        group.bench_with_input(
            BenchmarkId::new("horner", ell),
            &(&ej_poly, &roots),
            |b, (poly, r)| b.iter(|| black_box(poly).resultant_from_roots_horner(black_box(r))),
        );

        // (c) Use a 2^n balanced tree and absorb the last roots with Horner
        group.bench_with_input(
            BenchmarkId::new("mixed strategy", ell),
            &(&ej_poly, &balanced_tree, &rem_roots),
            |b, (poly, t, r)| {
                b.iter(|| {
                    black_box(poly).resultant_mixed_strategy_with_tree(black_box(t), black_box(r))
                })
            },
        );

        // (d) rebuild eval tree on every call expected slowest
        group.bench_with_input(
            BenchmarkId::new("remainder_tree_rebuild", ell),
            &(&ej_poly, &roots),
            |b, (poly, r)| b.iter(|| black_box(poly).resultant_from_roots(black_box(r))),
        );

        // (e) rebuild eval tree on every call but use horner for remainder
        group.bench_with_input(
            BenchmarkId::new("remainder_tree_mixed_rebuild", ell),
            &(&ej_poly, &roots),
            |b, (poly, r)| b.iter(|| black_box(poly).resultant_mixed_strategy(black_box(r))),
        );
    }
    group.finish();
}

// ===============================================================================
// Stage 4 hK evaluation
//
// Simple O(sK) product loop; expected negligible cost.  Benchmarked to rule
// out any surprises (e.g. if sK ends up larger than expected at some ell).
// ===============================================================================
fn bench_hk(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage4_hK");
    let mut rng = DRNG::from_seed(b"hk_bench");

    for &ell in &[103usize, 163, 211] {
        let (_, _, size_k) = partition(ell);

        let hk_points: Vec<PointX<Fp2>> = (0..size_k)
            .map(|_| PointX::new(&Fp2::rand(&mut rng), &Fp2::rand(&mut rng)))
            .collect();
        let xpz = Fp2::rand(&mut rng);
        let xmz = Fp2::rand(&mut rng);

        group.bench_with_input(
            BenchmarkId::new("hK_codomain", ell),
            &hk_points,
            |b, pts| b.iter(|| Curve::<Fp2>::hK_codomain(black_box(pts))),
        );

        group.bench_with_input(
            BenchmarkId::new("hK_eval", ell),
            &(&hk_points, xpz, xmz),
            |b, (pts, p, m)| {
                b.iter(|| Curve::<Fp2>::hK_eval(black_box(pts), black_box(p), black_box(m)))
            },
        );
    }
    group.finish();
}

// ===============================================================================
// Stage 5 full end-to-end isogeny (real kernel point, ell = 163)
//
// Four sub-benchmarks:
//   velu_codomain_only      traditional O(ell) velu, 0 images
//   sqrt_velu_codomain_only sqrt velu, 0 images
//   velu_one_image          traditional velu, 1 image point
//   sqrt_velu_one_image     sqrt velu, 1 image point
//
// Comparing (codomain_only) vs (one_image) gives the marginal cost of one
// image evaluation.  Comparing velu vs sqrt_velu at each row shows where the
// crossover actually sits for your implementation.
//
// The dummy image point is the kernel itself; it maps to infinity, which is a
// valid projective output we are measuring throughput, not correctness.
// ===============================================================================
fn bench_full_isogeny(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage5_full_isogeny");

    let ell = 163usize;
    let ker = PointX::new(&PX_163, &Fp2::ONE);
    let dummy_img = PointX::new(&PX_163, &Fp2::ONE);

    // velu, 0 images
    group.bench_function(BenchmarkId::new("velu_codomain_only", ell), |b| {
        b.iter(|| {
            let mut A24 = Fp2::TWO;
            let mut C24 = Fp2::FOUR;
            Curve::velu_odd_isogeny_proj(
                &mut A24,
                &mut C24,
                black_box(&ker),
                black_box(ell),
                &mut [],
            );
            (A24, C24)
        })
    });

    // sqrt velu, 0 images
    group.bench_function(BenchmarkId::new("sqrt_velu_codomain_only", ell), |b| {
        b.iter(|| {
            let mut A24 = Fp2::TWO;
            let mut C24 = Fp2::FOUR;
            Curve::sqrt_velu_odd_isogeny_proj::<PR>(
                &mut A24,
                &mut C24,
                black_box(&ker),
                black_box(ell),
                &mut [],
            );
            (A24, C24)
        })
    });

    // velu, 1 image
    group.bench_function(BenchmarkId::new("velu_one_image", ell), |b| {
        b.iter(|| {
            let mut A24 = Fp2::TWO;
            let mut C24 = Fp2::FOUR;
            let mut imgs = [dummy_img];
            Curve::velu_odd_isogeny_proj(
                &mut A24,
                &mut C24,
                black_box(&ker),
                black_box(ell),
                &mut imgs,
            );
            (A24, C24, imgs)
        })
    });

    // sqrt velu, 1 image
    group.bench_function(BenchmarkId::new("sqrt_velu_one_image", ell), |b| {
        b.iter(|| {
            let mut A24 = Fp2::TWO;
            let mut C24 = Fp2::FOUR;
            let mut imgs = [dummy_img];
            Curve::sqrt_velu_odd_isogeny_proj::<PR>(
                &mut A24,
                &mut C24,
                black_box(&ker),
                black_box(ell),
                &mut imgs,
            );
            (A24, C24, imgs)
        })
    });

    group.finish();
}

// ===============================================================================
// Stage 6 precompute_partitions
//
// Building the I/J/K point sets via the Montgomery ladder.  Should be a small
// fraction of the total but worth isolating from the polynomial work.
// ===============================================================================
fn bench_precompute_partitions(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage6_precompute_partitions");

    let ell = 163usize;
    let (size_j, size_i, size_k) = partition(ell);
    let ker = PointX::new(&PX_163, &Fp2::ONE);
    let A24 = Fp2::TWO;
    let C24 = Fp2::FOUR;

    group.bench_function(BenchmarkId::new("partition_points", ell), |b| {
        b.iter(|| {
            let mut hi = vec![PointX::INFINITY; size_i];
            let mut hj = vec![PointX::INFINITY; size_j];
            let mut hk = vec![PointX::INFINITY; size_k];
            Curve::<Fp2>::precompute_partitions(
                &mut hi,
                &mut hj,
                &mut hk,
                black_box(&A24),
                black_box(&C24),
                black_box(&ker),
            );
            (hi, hj, hk)
        })
    });

    group.finish();
}

// ===============================================================================
// Stage 7 precompute_eJ_values
//
// The per-isogeny setup that converts J points into the (sum_sqr, XZ4neg,
// AXZ4neg) triples reused across every image evaluation.  Should be cheap
// but is called once per isogeny, not once per image.
// ===============================================================================
fn bench_eJ_precompute(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage7_eJ_precompute");
    let mut rng = DRNG::from_seed(b"ej_precompute_bench");

    for &ell in &[103usize, 163, 211] {
        let (size_j, _, _) = partition(ell);

        let hj_points: Vec<PointX<Fp2>> = (0..size_j)
            .map(|_| PointX::new(&Fp2::rand(&mut rng), &Fp2::rand(&mut rng)))
            .collect();
        let A24 = Fp2::rand(&mut rng);
        let C24 = Fp2::rand(&mut rng);

        group.bench_with_input(
            BenchmarkId::new("precompute_eJ_values", ell),
            &(&hj_points, A24, C24),
            |b, (pts, a, cc)| {
                b.iter(|| {
                    let mut eJ = vec![(Fp2::ZERO, Fp2::ZERO, Fp2::ZERO); size_j];
                    Curve::<Fp2>::precompute_eJ_values(
                        &mut eJ,
                        black_box(pts),
                        black_box(a),
                        black_box(cc),
                    );
                    eJ
                })
            },
        );
    }
    group.finish();
}

// ===============================================================================

criterion_group!(
    benches,
    // bench_product_tree,
    // bench_eval_tree_construction,
    bench_resultant,
    // bench_hk,
    // bench_full_isogeny,
    // bench_precompute_partitions,
    // bench_eJ_precompute,
);
criterion_main!(benches);
