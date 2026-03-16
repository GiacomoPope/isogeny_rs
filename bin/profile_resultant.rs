use isogeny::fields::csidh::CsidhFp as Fp;
use isogeny::polynomial_ring::poly::{EvalTree, Polynomial};
use isogeny::utilities::test_utils::drng::DRNG;

fn main() {
    let mut rng = DRNG::from_seed("profile".as_bytes());

    // Use the sizes you actually benchmark
    let n_roots = 103;
    let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
    let et = EvalTree::new(&roots);

    // Run enough iterations to get a meaningful profile
    for _ in 0..10_000 {
        let f = Polynomial::rand(&mut rng, n_roots + 1);
        let _ = f.resultant_from_roots_with_tree(&et);
    }
}
