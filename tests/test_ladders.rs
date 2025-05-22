#![allow(non_snake_case)]

#[cfg(test)]
mod test_ladders {
    use fp2::traits::Fp as _;
    use isogeny::elliptic::basis::BasisX;
    use isogeny::elliptic::curve::Curve;
    use isogeny::fields::sike::SikeOne as Fp2;
    use isogeny::utilities::test_utils::drng::DRNG;
    use rand_core::RngCore;

    #[test]
    fn test_mul() {
        let mut rng = DRNG::from_seed("test_mul".as_bytes());

        let A = Fp2::from_i32(6);
        let E = Curve::new(&A);

        let mut scalar: [u8; 32] = [0; 32];
        rng.fill_bytes(&mut scalar);

        // Try truncating the scalar by a few bits by setting the bit length.
        for i in 0..20 {
            let P = E.rand_point(&mut rng);
            let Q = E.mul(&P, &scalar, (32 << 3) - i);

            let xP = P.to_pointx();
            let xQ = E.xmul(&xP, &scalar, (32 << 3) - i);

            assert!(xQ.equals(&Q.to_pointx()) == u32::MAX);
        }
    }

    #[test]
    fn test_3pt_ladder() {
        let mut rng = DRNG::from_seed("test_3pt_ladder".as_bytes());

        let A = Fp2::from_i32(6);
        let E = Curve::new(&A);

        let mut scalar: [u8; 32] = [0; 32];
        rng.fill_bytes(&mut scalar);

        // Try truncating the scalar by a few bits by setting the bit length.
        for i in 0..20 {
            // Compute P + [n]Q with projective points
            let P = E.rand_point(&mut rng);
            let Q = E.rand_point(&mut rng);
            let PQ = E.sub(&P, &Q);
            let nQ = E.mul(&Q, &scalar, (32 << 3) - i);
            let PnQ = E.add(&P, &nQ);

            // Compute P + [n]Q with x-only points
            let xP = P.to_pointx();
            let xQ = Q.to_pointx();
            let xPQ = PQ.to_pointx();
            let basis = BasisX::from_points(&xP, &xQ, &xPQ);
            let xPnQ = E.three_point_ladder(&basis, &scalar, (32 << 3) - i);

            // Ensure they're the same.
            assert!(xPnQ.equals(&PnQ.to_pointx()) == u32::MAX);
        }
    }

    #[test]
    fn test_ladder_biscalar() {
        let mut rng = DRNG::from_seed("test_biscalar_ladder".as_bytes());

        let A = Fp2::from_i32(6);
        let E = Curve::new(&A);

        // Try truncating the scalar by a few bits by setting the bit length.
        for i in 0..20 {
            let mut a: [u8; 32] = [0; 32];
            let mut b: [u8; 32] = [0; 32];
            rng.fill_bytes(&mut a);
            rng.fill_bytes(&mut b);

            // Compute [a]P + [b]Q with projective points
            let P = E.rand_point(&mut rng);
            let Q = E.rand_point(&mut rng);
            let PQ = E.sub(&P, &Q);

            let aP = E.mul(&P, &a, (32 << 3) - i);
            let bQ = E.mul(&Q, &b, (32 << 3) - 2 * i);
            let aPbQ = E.add(&aP, &bQ);

            // Compute  [a]P + [b]Q with x-only points
            let xP = P.to_pointx();
            let xQ = Q.to_pointx();
            let xPQ = PQ.to_pointx();
            let basis = BasisX::from_points(&xP, &xQ, &xPQ);
            let xaPbQ = E.ladder_biscalar(&basis, &a, &b, (32 << 3) - i, (32 << 3) - 2 * i);

            // Ensure they're the same.
            assert!(xaPbQ.equals(&aPbQ.to_pointx()) == u32::MAX);
        }
    }
}
