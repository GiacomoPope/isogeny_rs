mod benchmark_csidh {
    use isogeny::protocols::csidh_parameters::CSIDH_512;
    use isogeny::utilities::test_utils::drng::DRNG;

    use criterion::{Criterion, criterion_group};
    use std::time::Duration;

    fn benchmark_action(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("csidh_action".as_bytes());
        let bench_id = format!("Action for CSIDH-512 Parameters",);

        let (alice_sk, alice_pk) = CSIDH_512.keygen(&mut rng);

        // for simplicity, we just use alice pk and sk
        c.bench_function(&bench_id, |b| {
            b.iter(|| CSIDH_512.action(&alice_pk, &alice_sk, &mut rng))
        });
    }

    fn benchmark_verify(c: &mut Criterion) {
        let mut rng = DRNG::from_seed("csidh_action".as_bytes());
        let bench_id = format!("Verify for CSIDH-512 Parameters",);

        let (_, alice_pk) = CSIDH_512.keygen(&mut rng);

        // for simplicity, we just use alice pk and sk
        c.bench_function(&bench_id, |b| {
            b.iter(|| CSIDH_512.verify(&alice_pk, &mut rng))
        });
    }

    criterion_group! {
        name = csidh_benchmarks;
        config = Criterion::default().measurement_time(Duration::from_secs(15));
        targets = benchmark_action, benchmark_verify
    }
}

fn main() {
    benchmark_csidh::csidh_benchmarks();
}
