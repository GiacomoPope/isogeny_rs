use isogeny::fields::csidh::CsidhFp as Fp;
use isogeny::polynomial_ring::poly::{EvalTree, Polynomial};
use isogeny::utilities::test_utils::drng::DRNG;

fn main() {
    let mut rng = DRNG::from_seed("profile".as_bytes());

    let n_roots: usize = 10;
    let roots: Vec<Fp> = (0..n_roots).map(|_| Fp::rand(&mut rng)).collect();
    let n = usize::BITS - n_roots.leading_zeros() - 1;
    let two_n = 1 << n;
    let balanced_tree = EvalTree::new(&roots[..two_n]);

    let f = Polynomial::rand(&mut rng, n_roots + 1);
    let mut resultant = Fp::ONE;

    // Run enough iterations to get a meaningful profile
    for _ in 0..10_000 {
        for r in f.multieval_from_balanced_tree(&balanced_tree) {
            resultant *= r
        }

        for a in &roots[two_n..] {
            resultant *= f.evaluate(a)
        }
    }

    println!("{}", resultant);
}
