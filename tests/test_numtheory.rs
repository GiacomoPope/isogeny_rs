use isogeny::bigint::numtheory::cornacchia;
use isogeny::bigint::BigIntAlg;
use num_bigint::BigInt;

#[inline]
fn b(val: i32) -> BigInt {
    <BigInt as BigIntAlg>::from_i32(val)
}

#[test]
fn test_cornacchia_prime() {
    // Tests for cases where a solution to x^2 + n * y^2 = p exists
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

        // Verify x^2 + n * y^2 == p
        let x_sq = x.clone() * x.clone();
        let y_sq = y.clone() * y.clone();
        let prod = y_sq * n.clone();
        let res = x_sq + prod;
        assert_eq!(res, p, "Cornacchia output does not satisfy the Diophantine equation");
    }

    // Tests for cases where no solution exists.
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
