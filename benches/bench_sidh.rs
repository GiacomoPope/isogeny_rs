mod benchmark_sidh {
    use isogeny::protocols::sidh_parameters::SIDH_434;
    use isogeny::utilities::drng::DRNG;

    use criterion::{Criterion, black_box, criterion_group};
    use std::time::Duration;

    fn benchmark_keygen_alice(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("keygen_alice".as_bytes());
        let bench_id = format!("Benchmarking Alice Keygen for SIKE434 Parameters",);
        c.bench_function(&bench_id, |b| b.iter(|| SIDH_434.keygen_alice(&mut rng)));
    }

    fn benchmark_keygen_bob(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("keygen_bob".as_bytes());
        let bench_id = format!("Benchmarking Bob Keygen for SIKE434 Parameters",);
        c.bench_function(&bench_id, |b| b.iter(|| SIDH_434.keygen_bob(&mut rng)));
    }

    fn benchmark_secret_alice(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("alice shared_secret".as_bytes());
        let (_, alice_priv) = SIDH_434.keygen_alice(&mut rng);
        let (bob_pub, _) = SIDH_434.keygen_bob(&mut rng);

        let bench_id = format!("Benchmarking Alice Secret Generation for SIKE434 Parameters",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| black_box(alice_priv).shared_secret(black_box(&bob_pub)))
        });
    }

    fn benchmark_secret_bob(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("bob shared_secret".as_bytes());
        let (alice_pub, _) = SIDH_434.keygen_alice(&mut rng);
        let (_, bob_priv) = SIDH_434.keygen_bob(&mut rng);

        let bench_id = format!("Benchmarking Bob Secret Generation for SIKE434 Parameters",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| black_box(bob_priv).shared_secret(black_box(&alice_pub)))
        });
    }

    criterion_group! {
        name = sidh_benchmarks;
        config = Criterion::default().measurement_time(Duration::from_secs(3));
        targets = benchmark_keygen_alice, benchmark_keygen_bob, benchmark_secret_alice, benchmark_secret_bob
    }
}

fn main() {
    benchmark_sidh::sidh_benchmarks();
}
