use criterion::{Criterion, black_box, criterion_group, criterion_main};
use isogeny::protocols::sidh_parameters::SIDH_434;
use isogeny::utilities::test_utils::drng::DRNG;
use std::time::Duration;

fn bench_keygen(c: &mut Criterion) {
    let mut group = c.benchmark_group("sidh_keygen");

    let mut rng = DRNG::from_seed("keygen_alice".as_bytes());
    group.bench_function("alice", |b| b.iter(|| SIDH_434.keygen_alice(&mut rng)));

    let mut rng = DRNG::from_seed("keygen_bob".as_bytes());
    group.bench_function("bob", |b| b.iter(|| SIDH_434.keygen_bob(&mut rng)));

    group.finish();
}

fn bench_shared_secret(c: &mut Criterion) {
    let mut group = c.benchmark_group("sidh_shared_secret");

    let mut rng = DRNG::from_seed("alice_shared_secret".as_bytes());
    let (_, alice_priv) = SIDH_434.keygen_alice(&mut rng);
    let (bob_pub, _) = SIDH_434.keygen_bob(&mut rng);
    group.bench_function("alice", |b| {
        b.iter(|| black_box(alice_priv).shared_secret(black_box(&bob_pub)))
    });

    let mut rng = DRNG::from_seed("bob_shared_secret".as_bytes());
    let (alice_pub, _) = SIDH_434.keygen_alice(&mut rng);
    let (_, bob_priv) = SIDH_434.keygen_bob(&mut rng);
    group.bench_function("bob", |b| {
        b.iter(|| black_box(bob_priv).shared_secret(black_box(&alice_pub)))
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(3));
    targets = bench_keygen, bench_shared_secret
}
criterion_main!(benches);
