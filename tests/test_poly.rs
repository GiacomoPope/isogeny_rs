#[cfg(test)]
mod test_polynomial_arithmetic {

    use isogeny::fields::sqisign::SqiField248Base as Fp;
    use isogeny::{polynomial_ring::poly::Polynomial, utilities::test_utils::drng::DRNG};

    type PR = Polynomial<Fp>;

    #[test]
    fn test_addition() {
        let mut rng = DRNG::from_seed("polynomial_addition".as_bytes());

        let f = PR::rand(&mut rng, 5);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 5);

        let t1 = &(&f + &g) + &h;
        let t2 = &f + &(&g + &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 += g;
        t1 += h;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = PR::new_from_ele(&Fp::ZERO);
        let t1 = &f + &zero;
        assert!(t1.equals(&f) == u32::MAX);
    }

    #[test]
    fn test_subtraction() {
        let mut rng = DRNG::from_seed("polynomial_subtraction".as_bytes());

        let f = PR::rand(&mut rng, 5);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 5);

        let t1 = &(&f - &g) - &h;
        let t2 = &f - &(&g + &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 -= g;
        t1 -= h;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = &f - &f;
        assert!(zero.is_zero() == u32::MAX);
    }
}
