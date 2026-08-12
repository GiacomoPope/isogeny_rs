use isogeny::bigint::numtheory::{cornacchia, is_probable_prime};
use isogeny::bigint::BigIntAlg;
use num_bigint::BigInt;


#[inline]
fn b(val: i32) -> BigInt {
    <BigInt as BigIntAlg>::from_i32(val)
}

#[test]
fn test_cornacchia_prime() {
    // cases where a solution to x^2 + n * y^2 = p exists
    let solvable_pairs = vec![
        (b(1), b(5)),   // n=1, p=5
        (b(1), b(2)),   // n=1, p=2
        (b(1), b(41)),  // n=1, p=41
        (b(2), b(3)),   // n=2, p=3
        (b(3), b(7)),   // n=3, p=7
        (b(3), b(3)),   // n=3, p=3
    ];

    for (n, p) in solvable_pairs {
        let solution = cornacchia(&n, &p);
        assert!(
            solution.is_some(),
            "Cornacchia failed to find an existing solution for n={}, p={}",
            n, p
        );

        let (x, y) = solution.unwrap();

        // x^2 + n * y^2 == p
        let x_sq = x.clone() * x.clone();
        let y_sq = y.clone() * y.clone();
        let prod = y_sq * n.clone();
        let res = x_sq + prod;
        assert_eq!(res, p, "Cornacchia output does not satisfy the Diophantine equation");
    }

    // cases where no solution exists.
    let unsolvable_pairs = vec![
        (b(1), b(7)),   // n=1, p=7
        (b(1), b(3)),   // n=1, p=3
        (b(3), b(5)),   // n=3, p=5
    ];

    for (n, p) in unsolvable_pairs {
        let solution = cornacchia(&n, &p);
        assert!(
            solution.is_none(),
            "Cornacchia falsely found a solution for unsolvable case n={}, p={}",
            n, p
        );
    }
}

#[test]
fn test_is_probable_prime() {
    let primes = [b(2), b(3), b(5), b(7), b(11), b(103), b(1009)];
    let composites = [b(-5), b(0), b(1), b(4), b(6), b(9), b(100), b(1000)];

    for p in primes {
        assert!(is_probable_prime(&p, 20), "Failed to identify prime: {}", p);
    }

    for c in composites {
        assert!(!is_probable_prime(&c, 20), "Falsely identified composite as prime: {}", c);
    }
}

#[test]
fn test_cornacchia_robustness() {
    assert!(cornacchia(&b(1), &b(-5)).is_none(), "Cornacchia must reject negative targets");

    // TODO: maybe not
    assert!(cornacchia(&b(1), &b(15)).is_none(), "Cornacchia must reject composite targets");

    let solution = cornacchia(&b(1), &b(5));
    assert!(solution.is_some());
}
