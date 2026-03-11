#![allow(non_snake_case)]

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use isogeny::{
    elliptic::{curve::Curve, point::PointX},
    fields::csidh::CsidhFp as Fp,
    polynomial_ring::poly::Polynomial,
};
use std::time::Duration;

type P = Polynomial<Fp>;

// Kernel x-coordinates

const P103_X: Fp = Fp::const_decode_no_check(&[
    10, 101, 12, 18, 20, 58, 190, 183, 31, 106, 251, 46, 117, 208, 205, 0, 107, 11, 165, 131, 209,
    236, 94, 28, 131, 191, 90, 27, 169, 183, 177, 232, 184, 253, 81, 112, 173, 171, 3, 40, 82, 4,
    35, 58, 6, 72, 122, 207, 250, 143, 143, 158, 115, 63, 137, 255, 167, 228, 53, 34, 239, 74, 57,
    81,
]);
const P211_X: Fp = Fp::const_decode_no_check(&[
    238, 189, 22, 253, 65, 106, 138, 139, 190, 131, 156, 59, 59, 190, 21, 233, 174, 176, 192, 195,
    153, 179, 182, 105, 153, 207, 46, 81, 87, 20, 18, 219, 59, 234, 94, 76, 42, 161, 219, 64, 125,
    18, 249, 176, 192, 106, 185, 52, 186, 120, 3, 43, 102, 34, 142, 158, 116, 209, 128, 206, 38,
    220, 39, 95,
]);
const P311_X: Fp = Fp::const_decode_no_check(&[
    91, 231, 70, 194, 23, 155, 11, 243, 133, 79, 126, 163, 30, 195, 233, 66, 90, 73, 130, 99, 201,
    169, 156, 140, 207, 94, 104, 147, 102, 42, 34, 22, 17, 17, 125, 54, 214, 38, 252, 116, 52, 29,
    111, 43, 189, 39, 0, 97, 220, 143, 44, 183, 30, 50, 224, 4, 91, 217, 50, 188, 28, 248, 49, 96,
]);
const P587_X: Fp = Fp::const_decode_no_check(&[
    208, 223, 49, 61, 194, 244, 239, 9, 36, 20, 150, 189, 4, 98, 68, 49, 2, 82, 184, 250, 135, 61,
    188, 99, 169, 246, 32, 23, 64, 24, 96, 82, 187, 189, 77, 176, 221, 63, 218, 65, 29, 44, 31,
    159, 35, 25, 201, 79, 58, 228, 140, 42, 34, 42, 194, 48, 164, 47, 125, 226, 159, 124, 215, 19,
]);

// Expected Montgomery coefficients

const CODOMAIN_103: Fp = Fp::const_decode_no_check(&[
    176, 146, 168, 142, 34, 65, 67, 205, 144, 174, 85, 98, 64, 194, 29, 197, 131, 160, 208, 132,
    137, 89, 172, 179, 151, 224, 125, 200, 254, 84, 35, 248, 104, 58, 86, 252, 40, 50, 122, 63,
    188, 178, 136, 47, 150, 0, 168, 217, 2, 40, 73, 127, 138, 116, 51, 170, 208, 74, 245, 115, 18,
    190, 69, 62,
]);
const CODOMAIN_211: Fp = Fp::const_decode_no_check(&[
    40, 57, 173, 8, 5, 128, 178, 175, 193, 182, 157, 167, 18, 153, 250, 112, 96, 241, 172, 104,
    243, 64, 191, 104, 102, 71, 35, 164, 60, 231, 125, 213, 223, 178, 56, 133, 62, 134, 110, 39, 3,
    42, 145, 108, 148, 32, 164, 66, 79, 238, 83, 95, 74, 107, 129, 11, 190, 74, 178, 175, 60, 176,
    24, 64,
]);
const CODOMAIN_311: Fp = Fp::const_decode_no_check(&[
    89, 43, 146, 254, 39, 28, 73, 108, 102, 15, 144, 157, 34, 159, 168, 205, 141, 33, 0, 161, 47,
    38, 172, 204, 249, 232, 3, 79, 20, 13, 223, 160, 6, 158, 83, 181, 132, 161, 250, 131, 171, 44,
    4, 59, 150, 88, 131, 142, 156, 114, 77, 53, 244, 28, 62, 23, 59, 55, 56, 23, 63, 103, 151, 8,
]);
const CODOMAIN_587: Fp = Fp::const_decode_no_check(&[
    99, 164, 168, 164, 123, 19, 25, 132, 44, 91, 235, 107, 139, 228, 68, 154, 5, 32, 226, 199, 207,
    162, 164, 67, 6, 236, 167, 158, 121, 221, 59, 182, 25, 113, 68, 137, 43, 193, 177, 154, 93,
    238, 25, 71, 120, 131, 205, 202, 105, 110, 85, 248, 120, 170, 49, 163, 112, 192, 163, 235, 212,
    111, 68, 35,
]);

// Test data for benchmarks

struct IsogenyCase {
    ell: usize,
    ker_x: Fp,
    expected_codomain: Fp,
}

const CASES: &[IsogenyCase] = &[
    IsogenyCase {
        ell: 103,
        ker_x: P103_X,
        expected_codomain: CODOMAIN_103,
    },
    IsogenyCase {
        ell: 211,
        ker_x: P211_X,
        expected_codomain: CODOMAIN_211,
    },
    IsogenyCase {
        ell: 311,
        ker_x: P311_X,
        expected_codomain: CODOMAIN_311,
    },
    IsogenyCase {
        ell: 587,
        ker_x: P587_X,
        expected_codomain: CODOMAIN_587,
    },
];

/// Verify correctness of sqrt-vélu for a given case and return (ker, A24, C24).
fn setup_and_verify(case: &IsogenyCase) -> (PointX<Fp>, Fp, Fp) {
    let ker = PointX::new(&case.ker_x, &Fp::ONE);
    let mut images = [ker];
    let mut A24 = Fp::TWO;
    let mut C24 = Fp::FOUR;

    Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, case.ell, &mut images);

    let codomain = Curve::curve_from_A24_proj(&A24, &C24);
    let expected = Curve::new(&case.expected_codomain);
    assert!(codomain.j_invariant().equals(&expected.j_invariant()) == u32::MAX);
    assert!(images[0].Z.is_zero() == u32::MAX);

    (ker, A24, C24)
}

fn bench_velu(c: &mut Criterion) {
    let mut group = c.benchmark_group("velu");

    for case in CASES {
        let (ker, A24, C24) = setup_and_verify(case);

        group.bench_function(format!("ell={}", case.ell), |b| {
            b.iter(|| {
                Curve::velu_odd_isogeny_proj(
                    &mut black_box(A24),
                    &mut black_box(C24),
                    &black_box(ker),
                    black_box(case.ell),
                    &mut black_box([]),
                )
            })
        });
    }

    group.finish();
}

fn bench_sqrt_velu(c: &mut Criterion) {
    let mut group = c.benchmark_group("sqrt_velu");

    for case in CASES {
        let (ker, A24, C24) = setup_and_verify(case);

        group.bench_function(format!("ell={}", case.ell), |b| {
            b.iter(|| {
                Curve::sqrt_velu_odd_isogeny_proj::<P>(
                    &mut black_box(A24),
                    &mut black_box(C24),
                    &black_box(ker),
                    black_box(case.ell),
                    &mut black_box([]),
                )
            })
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(10));
    targets = bench_velu, bench_sqrt_velu
}
criterion_main!(benches);
