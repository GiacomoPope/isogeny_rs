use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use isogeny::quaternion::algebra::{IntQuat, QuatConfig};
use num_bigint::{BigInt, Sign};
use rand::{rng, RngCore};
use std::sync::LazyLock;

// ========================================================================
// SQISign Moduli Configurations
// ========================================================================

fn u64_slice_to_bigint(limbs: &[u64]) -> BigInt {
    let mut bytes = Vec::with_capacity(limbs.len() * 8);
    for &limb in limbs {
        bytes.extend_from_slice(&limb.to_le_bytes());
    }
    BigInt::from_bytes_le(Sign::Plus, &bytes)
}

// Level I
pub struct PSqisign1;
static P1_VAL: LazyLock<BigInt> = LazyLock::new(|| {
    u64_slice_to_bigint(&[
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x04FFFFFFFFFFFFFF,
    ])
});
impl QuatConfig<BigInt> for PSqisign1 {
    fn p() -> &'static BigInt {
        &P1_VAL
    }
}

// Level III
pub struct PSqisign3;
static P3_VAL: LazyLock<BigInt> = LazyLock::new(|| {
    u64_slice_to_bigint(&[
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x40FFFFFFFFFFFFFF,
    ])
});
impl QuatConfig<BigInt> for PSqisign3 {
    fn p() -> &'static BigInt {
        &P3_VAL
    }
}

// Level V
pub struct PSqisign5;
static P5_VAL: LazyLock<BigInt> = LazyLock::new(|| {
    u64_slice_to_bigint(&[
        0xFFFFFFFFFFFFFFFFu64,
        0xFFFFFFFFFFFFFFFFu64,
        0xFFFFFFFFFFFFFFFFu64,
        0xFFFFFFFFFFFFFFFFu64,
        0xFFFFFFFFFFFFFFFFu64,
        0xFFFFFFFFFFFFFFFFu64,
        0xFFFFFFFFFFFFFFFFu64,
        0x01AFFFFFFFFFFFFFu64,
    ])
});
impl QuatConfig<BigInt> for PSqisign5 {
    fn p() -> &'static BigInt {
        &P5_VAL
    }
}

// ========================================================================
// Multiplication Implementations
// ========================================================================

fn mul_naive<P: QuatConfig<BigInt>>(
    a: &IntQuat<BigInt, P>,
    b: &IntQuat<BigInt, P>,
) -> IntQuat<BigInt, P> {
    let p = P::p();
    let ac = &a.coords;
    let bc = &b.coords;

    let mut c0 = ac[0].clone() * bc[0].clone() - ac[1].clone() * bc[1].clone();
    c0 = c0 - p.clone() * (ac[2].clone() * bc[2].clone() + ac[3].clone() * bc[3].clone());

    let mut c1 = p.clone() * (ac[2].clone() * bc[3].clone() - ac[3].clone() * bc[2].clone());
    c1 = c1 + ac[0].clone() * bc[1].clone() + ac[1].clone() * bc[0].clone();

    let c2 = ac[0].clone() * bc[2].clone() + ac[2].clone() * bc[0].clone()
        - ac[1].clone() * bc[3].clone() + ac[3].clone() * bc[1].clone();

    let c3 = ac[0].clone() * bc[3].clone() + ac[3].clone() * bc[0].clone()
        - ac[2].clone() * bc[1].clone() + ac[1].clone() * bc[2].clone();

    IntQuat::new(c0, c1, c2, c3)
}

fn mul_karatsuba<P: QuatConfig<BigInt>>(
    a: &IntQuat<BigInt, P>,
    b: &IntQuat<BigInt, P>,
) -> IntQuat<BigInt, P> {
    let p = P::p();
    let ac = &a.coords;
    let bc = &b.coords;

    let m1 = ac[0].clone() * bc[0].clone();
    let m2 = ac[1].clone() * bc[1].clone();
    let m3 = (ac[0].clone() + ac[1].clone()) * (bc[0].clone() + bc[1].clone());

    let m4 = ac[2].clone() * bc[2].clone();
    let m5 = ac[3].clone() * bc[3].clone();
    let m6 = (ac[2].clone() + ac[3].clone()) * (bc[2].clone() - bc[3].clone());

    let m7 = ac[0].clone() * bc[2].clone();
    let m8 = ac[1].clone() * bc[3].clone();
    let m9 = (ac[0].clone() + ac[1].clone()) * (bc[2].clone() + bc[3].clone());

    let m10 = ac[2].clone() * bc[0].clone();
    let m11 = ac[3].clone() * bc[1].clone();
    let m12 = (ac[2].clone() + ac[3].clone()) * (bc[0].clone() - bc[1].clone());

    let c0 = (m1.clone() - m2.clone()) - p.clone() * (m4.clone() + m5.clone());
    let c1 = (m3 - m1.clone() - m2.clone()) - p.clone() * (m6 - m4 + m5);
    let c2 = (m7.clone() - m8.clone()) + (m10.clone() + m11.clone());
    let c3 = (m9 - m7 - m8) + (m12 - m10 + m11);

    IntQuat::new(c0, c1, c2, c3)
}

// ========================================================================
// Benchmarking Logic
// ========================================================================

fn generate_random_quat<P: QuatConfig<BigInt>>(bits: u64) -> IntQuat<BigInt, P> {
    let bytes_len = (bits / 8) as usize;
    let mut c0 = vec![0u8; bytes_len];
    let mut c1 = vec![0u8; bytes_len];
    let mut c2 = vec![0u8; bytes_len];
    let mut c3 = vec![0u8; bytes_len];

    let mut r = rng();
    r.fill_bytes(&mut c0);
    r.fill_bytes(&mut c1);
    r.fill_bytes(&mut c2);
    r.fill_bytes(&mut c3);

    IntQuat::new(
        BigInt::from_bytes_le(Sign::Plus, &c0),
        BigInt::from_bytes_le(Sign::Plus, &c1),
        BigInt::from_bytes_le(Sign::Plus, &c2),
        BigInt::from_bytes_le(Sign::Plus, &c3),
    )
}

fn bench_for_level<P: QuatConfig<BigInt>>(c: &mut Criterion, level_name: &str) {
    let mut group = c.benchmark_group(format!("Quaternion Multiplication ({})", level_name));

    let bit_sizes = [256, 1024, 2048, 4096, 8192];

    for &size in bit_sizes.iter() {
        let q1: IntQuat<BigInt, P> = generate_random_quat(size);
        let q2: IntQuat<BigInt, P> = generate_random_quat(size);

        let res_naive = mul_naive(&q1, &q2);
        let res_kara = mul_karatsuba(&q1, &q2);
        assert!(
            res_naive.coords == res_kara.coords,
            "[DEBUG] Karatsuba multiplication diverged from Naive for bit size {} in {}!",
            size, level_name
        );

        group.bench_with_input(BenchmarkId::new("Naive (16 Mul)", size), &size, |b, _| {
            b.iter(|| mul_naive(black_box(&q1), black_box(&q2)))
        });

        group.bench_with_input(BenchmarkId::new("Karatsuba (12 Mul)", size), &size, |b, _| {
            b.iter(|| mul_karatsuba(black_box(&q1), black_box(&q2)))
        });
    }

    group.finish();
}

fn bench_quat_mul(c: &mut Criterion) {
    bench_for_level::<PSqisign1>(c, "Level I");
    bench_for_level::<PSqisign3>(c, "Level III");
    bench_for_level::<PSqisign5>(c, "Level V");
}

criterion_group!(benches, bench_quat_mul);
criterion_main!(benches);
