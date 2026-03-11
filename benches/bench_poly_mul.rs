use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};

use isogeny::fields::csidh::CsidhFp as Fp;
use isogeny::{polynomial_ring::poly::Polynomial, utilities::test_utils::drng::DRNG};
type PR = Polynomial<Fp>;

const F_DEGREES: &[usize] = &[4, 8, 12, 16, 20, 40, 80];
const G_DEGREES: &[usize] = &[1, 2, 4, 8, 16, 32];

fn bench_basic_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("basic_mul");
    let mut rng = DRNG::from_seed("bench_poly_mul".as_bytes());

    for &f_deg in F_DEGREES {
        for &g_deg in G_DEGREES {
            let f = PR::rand(&mut rng, f_deg);
            let g = PR::rand(&mut rng, g_deg);

            group.bench_with_input(
                BenchmarkId::new("f_deg_g_deg", format!("{f_deg}x{g_deg}")),
                &(f_deg, g_deg),
                |b, _| b.iter(|| black_box(f.basic_mul(black_box(&g)))),
            );
        }
    }

    group.finish();
}

fn bench_karatsuba_mul(c: &mut Criterion) {
    let mut group = c.benchmark_group("karatsuba_mul");
    let mut rng = DRNG::from_seed("bench_poly_mul".as_bytes());

    for &f_deg in F_DEGREES {
        for &g_deg in G_DEGREES {
            let f = PR::rand(&mut rng, f_deg);
            let g = PR::rand(&mut rng, g_deg);

            group.bench_with_input(
                BenchmarkId::new("f_deg_g_deg", format!("{f_deg}x{g_deg}")),
                &(f_deg, g_deg),
                |b, _| b.iter(|| black_box(f.karatsuba_mul(black_box(&g)))),
            );
        }
    }

    group.finish();
}

criterion_group!(benches, bench_basic_mul, bench_karatsuba_mul);
criterion_main!(benches);
