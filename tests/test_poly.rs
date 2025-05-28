#[cfg(test)]
mod test_polynomial_arithmetic {

    use isogeny::fields::sqisign::SqiField248Base as Fp;
    use isogeny::{polynomial_ring::poly::Polynomial, utilities::test_utils::drng::DRNG};

    type PR = Polynomial<Fp>;

    #[test]
    fn test_addition() {
        let mut rng = DRNG::from_seed("polynomial_addition".as_bytes());

        let f = PR::rand(&mut rng, 4);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 6);

        let t1 = &(&f + &g) + &h;
        let t2 = &f + &(&g + &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 += &g;
        t1 += &h;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = PR::new_from_ele(&Fp::ZERO);
        let t1 = &f + &zero;
        assert!(t1.equals(&f) == u32::MAX);

        let f_neg = -&f;
        assert!((&f + &f_neg).is_zero() == u32::MAX);
    }

    #[test]
    fn test_subtraction() {
        let mut rng = DRNG::from_seed("polynomial_subtraction".as_bytes());

        let f = PR::rand(&mut rng, 4);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 6);

        let t1 = &(&f - &g) - &h;
        let t2 = &f - &(&g + &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 -= &g;
        t1 -= &h;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = &f - &f;
        assert!(zero.is_zero() == u32::MAX);
    }

    #[test]
    fn test_multiplication() {
        let mut rng = DRNG::from_seed("polynomial_multiplication".as_bytes());

        let f = PR::rand(&mut rng, 4);
        let g = PR::rand(&mut rng, 5);
        let h = PR::rand(&mut rng, 6);

        let t1 = &(&f * &g) * &h;
        let t2 = &f * &(&g * &h);

        assert!(t1.equals(&t2) == u32::MAX);

        let mut t1 = f.clone();
        t1 *= &g;
        t1 *= &h;
        assert!(t1.equals(&t2) == u32::MAX);

        let one = PR::new_from_ele(&Fp::ONE);
        t1 *= &one;
        assert!(t1.equals(&t2) == u32::MAX);

        let zero = PR::new_from_ele(&Fp::ZERO);
        t1 *= &zero;
        assert!(t1.is_zero() == u32::MAX);
    }

    #[test]
    fn test_squaring() {
        // TODO: do specalised squaring method?
        let mut rng = DRNG::from_seed("polynomial_squaring".as_bytes());

        let f = PR::rand(&mut rng, 5);

        let t1 = &f * &f;
        let mut t2 = f.clone();
        t2 *= &f;

        assert!(t1.equals(&t2) == u32::MAX);
    }

    #[test]
    fn test_scaling() {
        let mut rng = DRNG::from_seed("polynomial_scaling".as_bytes());

        let f = PR::rand(&mut rng, 5);

        // Compute [4]f in different ways
        let mut t1 = f.clone();
        t1 += &f;
        t1 += &f;
        t1 += &f;

        let t2 = f.scale(&Fp::FOUR);
        assert!(t1.equals(&t2) == u32::MAX);

        let t2 = f.scale_small(4);
        assert!(t1.equals(&t2) == u32::MAX);

        let g = f.scale(&Fp::ONE);
        assert!(f.equals(&g) == u32::MAX);

        let g = f.scale_small(1);
        assert!(f.equals(&g) == u32::MAX);

        let g = f.scale(&Fp::ZERO);
        assert!(g.is_zero() == u32::MAX);

        let g = f.scale_small(0);
        assert!(g.is_zero() == u32::MAX);
    }

    #[test]
    fn test_evaluation() {
        let mut rng = DRNG::from_seed("polynomial_eval".as_bytes());

        // Evaluating a constant polynomial should always give the same result
        let a = Fp::rand(&mut rng);
        let b = Fp::rand(&mut rng);
        let f = PR::rand(&mut rng, 1);
        let ea = f.evaluate(&a);
        let eb = f.evaluate(&b);
        assert!(ea.equals(&eb) == u32::MAX);

        // Random degree four polynomial
        let f = PR::rand(&mut rng, 5);

        // Evaluating at zero should give you the constant coefficient
        let e0 = f.evaluate(&Fp::ZERO);
        let c0 = f.constant_coefficient().unwrap();
        assert!(e0.equals(&c0) == u32::MAX);

        // Evaluating at one should give you the sum of the coefficients
        let e1 = f.evaluate(&Fp::ONE);
        let c1 = f[0] + f[1] + f[2] + f[3] + f[4];
        assert!(e1.equals(&c1) == u32::MAX);

        // Dumb evaluation to compare against the type method.
        let a = Fp::rand(&mut rng);
        let ea = f.evaluate(&a);
        let mut test_ea = f[0];
        test_ea += f[1] * a;
        test_ea += f[2] * (a * a);
        test_ea += f[3] * (a * a * a);
        test_ea += f[4] * (a * a * a * a);

        assert!(ea.equals(&test_ea) == u32::MAX);
    }

    #[test]
    fn test_degree() {
        let f = PR::new_from_slice(&[Fp::ZERO]);
        assert!(f.degree() == None);

        let f = PR::new_from_slice(&[Fp::ONE]);
        assert!(f.degree().unwrap() == 0);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ONE]);
        assert!(f.degree().unwrap() == 1);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ONE, Fp::ZERO]);
        assert!(f.degree().unwrap() == 4);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO]);
        assert!(f.degree() == None);

        let f = PR::new_from_slice(&[Fp::ONE, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO]);
        assert!(f.degree().unwrap() == 0);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ONE, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO]);
        assert!(f.degree().unwrap() == 1);

        let f = PR::new_from_slice(&[Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ZERO, Fp::ONE, Fp::ZERO]);
        assert!(f.degree().unwrap() == 4);
    }
}
