#![allow(non_snake_case)]

mod benchmark_velu {
    use criterion::{Criterion, black_box, criterion_group};
    use isogeny::{
        elliptic::{curve::Curve, point::PointX},
        polynomial_ring::poly::Polynomial,
    };
    use std::time::Duration;

    // Toy prime, with smooth p + 1
    // p + 1 = 2^3 * 3^2 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97 * 163^16
    const MODULUS: [u64; 4] = [
        0x65396406224d1bc7,
        0x60892bae5986e2f3,
        0x8f2e6d4fb53f8cc7,
        0x0003e368df883c71,
    ];
    fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = MODULUS,);

    type P = Polynomial<Fp2>;

    // Kernel of order 163
    const P_X_RE_BYTES: [u8; 31] = [
        161, 11, 202, 60, 238, 183, 30, 198, 113, 154, 4, 46, 251, 33, 116, 68, 50, 187, 152, 163,
        215, 172, 206, 196, 15, 96, 128, 81, 141, 27, 2,
    ];
    const P_X_IM_BYTES: [u8; 31] = [
        244, 36, 105, 80, 212, 91, 28, 64, 95, 113, 78, 111, 165, 12, 140, 76, 43, 91, 43, 245, 47,
        7, 52, 200, 71, 10, 182, 128, 100, 211, 1,
    ];
    const AP_RE_BYTES: [u8; 31] = [
        81, 48, 37, 31, 25, 41, 102, 139, 211, 247, 169, 227, 160, 80, 111, 30, 149, 218, 159, 209,
        161, 29, 168, 254, 63, 197, 31, 62, 144, 178, 0,
    ];
    const AP_IM_BYTES: [u8; 31] = [
        124, 122, 247, 73, 156, 36, 71, 39, 84, 47, 89, 103, 72, 193, 196, 240, 28, 139, 178, 131,
        241, 80, 164, 147, 32, 114, 247, 55, 117, 48, 3,
    ];

    const PX: Fp2 = Fp2::const_decode_no_check(&P_X_RE_BYTES, &P_X_IM_BYTES);
    const CODOMAIN_AP: Fp2 = Fp2::const_decode_no_check(&AP_RE_BYTES, &AP_IM_BYTES);

    fn bench_velu(c: &mut Criterion) {
        let ker = PointX::new(&PX, &Fp2::ONE);
        let codomain_test = Curve::new(&CODOMAIN_AP);
        let mut images = [ker];

        let mut A24 = Fp2::TWO;
        let mut C24 = Fp2::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 163, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);

        let bench_id = format!("Benchmarking 163-isogeny using velu",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| {
                Curve::velu_odd_isogeny_proj(
                    &mut black_box(A24),
                    &mut black_box(C24),
                    &black_box(ker),
                    black_box(163),
                    &mut black_box([]),
                );
            })
        });
    }

    fn bench_sqrt_velu(c: &mut Criterion) {
        let ker = PointX::new(&PX, &Fp2::ONE);
        let codomain_test = Curve::new(&CODOMAIN_AP);
        let mut images = [ker];

        let mut A24 = Fp2::TWO;
        let mut C24 = Fp2::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 163, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);

        let bench_id = format!("Benchmarking 163-isogeny using sqrt velu",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| {
                Curve::sqrt_velu_odd_isogeny_proj::<P>(
                    &mut black_box(A24),
                    &mut black_box(C24),
                    &black_box(ker),
                    black_box(163),
                    &mut black_box([]),
                );
            })
        });
    }

    fn bench_inversion(c: &mut Criterion) {
        let bench_id = format!("Benchmarking inversion...",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| {
                black_box(PX).invert();
            })
        });
    }

    fn bench_varmul(c: &mut Criterion) {
        let ker = PointX::new(&PX, &Fp2::ONE);
        let A24 = Fp2::TWO;
        let C24 = Fp2::FOUR;
        let v = 6;

        let bench_id = format!("Benchmarking multiplication...",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| {
                Curve::xmul_proj_u64_vartime(
                    black_box(&A24),
                    black_box(&C24),
                    black_box(&ker),
                    black_box((v + v) as u64),
                )
            })
        });
    }

    criterion_group! {
        name = benchmark_velu;
        config = Criterion::default().measurement_time(Duration::from_secs(10));
        targets = bench_velu, bench_sqrt_velu, bench_inversion, bench_varmul
    }
}

fn main() {
    benchmark_velu::benchmark_velu();
}
