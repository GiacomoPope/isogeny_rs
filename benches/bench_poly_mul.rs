use criterion::{Criterion, black_box, criterion_group, criterion_main};

use isogeny::fields::csidh::CsidhFp as Fp;
use isogeny::{
    polynomial_ring::poly::{EvalTree, Polynomial},
    utilities::test_utils::drng::DRNG,
};
type PR = Polynomial<Fp>;

// const F_DEGREES: &[usize] = &[23];
// const G_DEGREES: &[usize] = &[12];

// fn bench_basic_mul(c: &mut Criterion) {
//     let mut group = c.benchmark_group("basic_mul");
//     let mut rng = DRNG::from_seed("bench_poly_mul".as_bytes());

//     for &f_deg in F_DEGREES {
//         for &g_deg in G_DEGREES {
//             let f = PR::rand(&mut rng, f_deg);
//             let g = PR::rand(&mut rng, g_deg);

//             group.bench_with_input(
//                 BenchmarkId::new("f_deg_g_deg", format!("{f_deg}x{g_deg}")),
//                 &(f_deg, g_deg),
//                 |b, _| b.iter(|| black_box(f.basic_mul(black_box(&g)))),
//             );
//         }
//     }

//     group.finish();
// }

// fn bench_karatsuba_mul(c: &mut Criterion) {
//     let mut group = c.benchmark_group("karatsuba_mul");
//     let mut rng = DRNG::from_seed("bench_poly_mul".as_bytes());

//     for &f_deg in F_DEGREES {
//         for &g_deg in G_DEGREES {
//             let f = PR::rand(&mut rng, f_deg);
//             let g = PR::rand(&mut rng, g_deg);

//             group.bench_with_input(
//                 BenchmarkId::new("f_deg_g_deg", format!("{f_deg}x{g_deg}")),
//                 &(f_deg, g_deg),
//                 |b, _| b.iter(|| black_box(f.karatsuba_mul(black_box(&g)))),
//             );
//         }
//     }

//     group.finish();
// }

// fn bench_mul_low(c: &mut Criterion) {
//     let mut group = c.benchmark_group("poly_mul_low");
//     let mut rng = DRNG::from_seed("bench_poly_mul_low".as_bytes());

//     for &f_deg in F_DEGREES {
//         for &g_deg in G_DEGREES {
//             let f = PR::rand(&mut rng, f_deg);
//             let g = PR::rand(&mut rng, g_deg);
//             let n = (f_deg + g_deg) / 2;

//             group.bench_with_input(
//                 BenchmarkId::new("f_deg_g_deg", format!("{f_deg}x{g_deg}")),
//                 &(f_deg, g_deg),
//                 |b, _| b.iter(|| black_box(f.mul_low_test(black_box(&g), n))),
//             );
//         }
//     }

//     group.finish();
// }

// fn bench_mul_middle(c: &mut Criterion) {
//     let mut group = c.benchmark_group("poly_mul_middle");
//     let mut rng = DRNG::from_seed("bench_poly_mul_middle".as_bytes());

//     for &f_deg in F_DEGREES {
//         for &g_deg in G_DEGREES {
//             assert!(2 * g_deg - 1 == f_deg);
//             let f = PR::rand(&mut rng, f_deg);
//             let g = PR::rand(&mut rng, g_deg);

//             group.bench_with_input(
//                 BenchmarkId::new("f_deg_g_deg", format!("{f_deg}x{g_deg}")),
//                 &(f_deg, g_deg),
//                 |b, _| b.iter(|| black_box(f.mul_middle_test(black_box(&g)))),
//             );
//         }
//     }

//     group.finish();
// }

// // inv_mod_xn for n_pad sizes you actually hit
// fn bench_inv_mod_xn(c: &mut Criterion) {
//     let mut rng = DRNG::from_seed("bench_inv".as_bytes());
//     let mut group = c.benchmark_group("inv_mod_xn");
//     for &n in &[16usize, 32, 64] {
//         let f = PR::rand(&mut rng, n + 3);
//         group.bench_with_input(BenchmarkId::new("n", n), &n, |b, &n| {
//             let mut h = vec![Fp::ZERO; n];
//             b.iter(|| PR::inv_mod_xn(&mut h, &f.coeffs, n));
//         });
//     }
//     group.finish();
// }

fn bench_sqrt_velu_comparison(c: &mut Criterion) {
    let mut rng = DRNG::from_seed("sqrt_velu_comparison".as_bytes());
    let mut group = c.benchmark_group("sqrt_velu_comparison");

    let deg = 10;
    let n_roots = 10;

    let f = PR::rand(&mut rng, deg + 1);
    let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
    let eval_tree = EvalTree::new(&roots);

    group.bench_function("product_tree", |b: &mut criterion::Bencher<'_>| {
        b.iter(|| PR::product_tree_from_roots_simple(&roots))
    });

    // Naive O(n^2): evaluate polynomial at each root using Horner, done 4 times
    // (codomain r0, r1 + one image point r0, r1)
    group.bench_function("naive_evaluation", |b| {
        b.iter(|| {
            let mut res: Fp = Fp::ONE;
            for r in roots.iter() {
                res *= black_box(&f).evaluate(black_box(r));
            }

            res
        })
    });

    // New method part 1: build EvalTree once
    group.bench_function("eval_tree_padded", |b| b.iter(|| EvalTree::new(&roots)));

    // New method part 2: single remainder tree eval
    group.bench_function("remainder_tree_eval_x1", |b| {
        b.iter(|| black_box(&f).resultant_from_roots_with_tree(black_box(&eval_tree)))
    });

    group.finish();
}

criterion_group!(
    benches,
    // bench_basic_mul,
    // bench_karatsuba_mul,
    // bench_mul_low,
    // bench_mul_middle,
    // bench_inv_mod_xn
    bench_sqrt_velu_comparison
);
criterion_main!(benches);
