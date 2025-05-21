#![allow(non_snake_case)]

mod benchmark_biscalar {
    use std::time::Duration;

    use criterion::{Criterion, black_box, criterion_group};
    use fp2::fq::Fq as FqTrait;
    use rand_core::RngCore;

    use isogeny::elliptic::basis::BasisX;
    use isogeny::elliptic::curve::Curve;
    use isogeny::fields::sqisign::SqiField248 as Fp2;
    use isogeny::utilities::test_utils::drng::DRNG;

    fn benchmark_ladder_biscalar(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("test_biscalar_ladder".as_bytes());

        let A = Fp2::from_i32(6);
        let E = Curve::new(&A);

        let mut a: [u8; 32] = [0; 32];
        let mut b: [u8; 32] = [0; 32];
        rng.fill_bytes(&mut a);
        rng.fill_bytes(&mut b);

        // Compute [a]P + [b]Q with projective points
        let P = E.rand_point(&mut rng);
        let Q = E.rand_point(&mut rng);
        let PQ = E.sub(&P, &Q);

        let aP = E.mul(&P, &a, 32 << 3);
        let bQ = E.mul(&Q, &b, 32 << 3);
        let aPbQ = E.add(&aP, &bQ);

        // Compute  [a]P + [b]Q with x-only points
        let xP = P.to_pointx();
        let xQ = Q.to_pointx();
        let xPQ = PQ.to_pointx();
        let basis = BasisX::from_points(&xP, &xQ, &xPQ);
        let xaPbQ = E.ladder_biscalar(&basis, &a, &b, 32 << 3, 32 << 3);

        // Ensure they're the same.
        assert!(xaPbQ.equals(&aPbQ.to_pointx()) == u32::MAX);

        let bench_id = format!("Benchmarking biscalar ladder with 128 bit scalar",);
        c.bench_function(&bench_id, |bb| {
            bb.iter(|| {
                black_box(E).ladder_biscalar(
                    &black_box(basis),
                    &black_box(a),
                    &black_box(b),
                    black_box(32 << 3),
                    black_box(32 << 3),
                )
            })
        });
    }

    criterion_group! {
        name = benchmark_biscalar;
        config = Criterion::default().measurement_time(Duration::from_secs(10));
        targets = benchmark_ladder_biscalar,
    }
}

fn main() {
    benchmark_biscalar::benchmark_biscalar();
}
