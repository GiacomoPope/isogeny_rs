#![allow(non_snake_case)]

#[cfg(test)]
mod test_product_isogeny {
    use fp2::fq::Fq;
    use isogeny::elliptic::curve::Curve;
    use isogeny::elliptic::point::PointX;
    use isogeny::elliptic::product::{CouplePoint, EllipticProduct};
    use isogeny::fields::sqisign::SqiSqignI as Fp2;
    use isogeny::theta::theta_chain::product_isogeny_no_strategy;

    static A1_STR: &str = "baa79b0dc07508bb6fea4685db4b48f237686ad1e12964c985814261bcee97015c28136967c8faa77df1d28ffe4f81c68369091bc503d27ab21c459ce88eb101";
    static A2_STR: &str = "84228651f271b0f39f2f19f2e8718f31ed3365ac9e5cb303afe663d0cfc11f0455d891b0ca6c7e653f9ba2667730bb77befe1b1a31828404284af8fd7baacc01";

    static X_P1_STR: &str = "d1bbcf7e9058a9c789ccbd1f536171efa1cc3c02e400b56c7b217ad7fbf8f700980fc332a818e8eac1e09000d474da7db0f1fd01ebb60fa64fc98da3647ff304";

    static X_P2_STR: &str = "f2c19421bd409b6a83d039a6f15026c8335be84603ec7fc8f01f65ca7ce416033606f5b2cc72c95a366e3e9fc346278ecf32aadfe898ad08f40971162e1b5701";

    static X_Q1_STR: &str = "2d1a92124027ec0e42d5e2a52ea7eff883c112d1ed20ec6bb343952de47cb602b992487573c376416df6de77de9ecdd2ff2099c81ca3cc48fc6c5146a83e4d01";

    static X_Q2_STR: &str = "b8a50598336df416c591022f842f40be51018073f82279a99ecaa95cbad91200477cf67ba9fb4d129b50022de7952a6f5e7648aca8ca68daabf46881754fd601";

    static A3_STR: &str = "c9ae315618b89d82c3ba2c2f6ba312b3fd7e67caf4ed47ff22f826f775b7af0178d3a7399d5e086b1913535fe9931d94511c79715e0f9dffbba43a01059a2d00";
    static A4_STR: &str = "2b14b4831c280c523790d20739de9502f35d4d44b4d5abad52d6d15b4bc21a0342a8c6abdf9a7f2fee3be93ca5f71f425051664645eb54af4391682d650bd100";

    #[test]
    fn test_codomain_only() {
        let (A1, check) = Fp2::decode(&hex::decode(A1_STR).unwrap());
        assert!(check == u32::MAX);

        let (A2, check) = Fp2::decode(&hex::decode(A2_STR).unwrap());
        assert!(check == u32::MAX);

        let (X_P1, check) = Fp2::decode(&hex::decode(X_P1_STR).unwrap());
        assert!(check == u32::MAX);

        let (X_P2, check) = Fp2::decode(&hex::decode(X_P2_STR).unwrap());
        assert!(check == u32::MAX);

        let (X_Q1, check) = Fp2::decode(&hex::decode(X_Q1_STR).unwrap());
        assert!(check == u32::MAX);

        let (X_Q2, check) = Fp2::decode(&hex::decode(X_Q2_STR).unwrap());
        assert!(check == u32::MAX);

        let (A3, check) = Fp2::decode(&hex::decode(A3_STR).unwrap());
        assert!(check == u32::MAX);

        let (A4, check) = Fp2::decode(&hex::decode(A4_STR).unwrap());
        assert!(check == u32::MAX);

        // Domain
        let E1 = Curve::new(&A1);
        let E2 = Curve::new(&A2);
        let E1E2 = EllipticProduct::new(&E1, &E2);

        // Kernel
        let xP1 = PointX::new(&X_P1, &Fp2::ONE);
        let xP2 = PointX::new(&X_P2, &Fp2::ONE);
        let xQ1 = PointX::new(&X_Q1, &Fp2::ONE);
        let xQ2 = PointX::new(&X_Q2, &Fp2::ONE);

        let (P1, _) = E1.complete_pointX(&xP1);
        let (P2, _) = E2.complete_pointX(&xP2);
        let (Q1, _) = E1.complete_pointX(&xQ1);
        let (Q2, _) = E2.complete_pointX(&xQ2);

        let P1P2 = CouplePoint::new(&P1, &P2);
        let Q1Q2 = CouplePoint::new(&Q1, &Q2);

        // Codomain
        let E3 = Curve::new(&A3);
        let E4 = Curve::new(&A4);

        // Compute chain
        let (E3E4, _) = product_isogeny_no_strategy(&E1E2, &P1P2, &Q1Q2, &[], 125);

        let (F3, F4) = E3E4.curves();

        println!("{}", E3.j_invariant());
        println!("{}", F3.j_invariant());
        println!("{}", E4.j_invariant());
        println!("{}", F4.j_invariant());

        assert!(1 == 2);
    }
}
