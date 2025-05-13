#![allow(non_snake_case)]

#[cfg(test)]
mod test_three_isogeny_chain {
    use fp2::fq::Fq;
    use isogeny::elliptic::curve::Curve;
    use isogeny::elliptic::isomorphism::Isomorphism;
    use isogeny::elliptic::point::PointX;
    use isogeny::elliptic::three_isogeny_chain::{three_isogeny_chain, three_isogeny_chain_naive};
    use isogeny::fields::sike::SikeOne as Fp2;

    // Kernel for phi_1 : E -> E/<P> with P not above (0 : 1)
    static XP_STR: &str = "e8878b984636db6628d4d438b1418c92be68ba5dd27136422a10594358bf4592470fb762c70398302e560324c9941451484945a30164001241408c9f6284b5c0a4995cd1ea6546e2603cfebb86192249a8813152a53df62fab3164af9a949b8188c1baaf9450fe1e34673e04c200";

    // Kernel for phi_2 : E -> E/<P> with Q above (0 : 1)
    static XQ_STR: &str = "fbe8a2d4e409b26f4462fab6e439f0cf0e3a6f24b398810a862508ba93d636930250811d1a0d3e4cf36a3726256cdb1fda1d5bf63a9e01f3071c57a99401309cdba3d62c94c66f8c02bd4e85257da778d63750ab5ff7f841462a8a7e081b7cd5ada2b986e709e1df9f9b08d41701";

    // Codomain for phi_1 : E -> E/<P>
    static A1_STR: &str = "c127942edd111648ab5abb229451eed6b24b967c2a813e0fa7c2ff4ea5bba35d8406135fe1ddd70beb07f48a4683c4bdd5bf6accc68a0140b5953c1478bf68741d1aa9fe4b4788b120dd75ff23c57c8c8d633eafb89a69360b652009bcdf39ac4ccc247cce89b850d3115a680002";

    // Codomain for phi_2 : E -> E/<Q>
    static A2_STR: &str = "9d36ba828e2e085052060a0397c31f7683d4677c3425837d7ac0a01fd91f5b994cc997f5b2363eaca8ee19dcb928160645858148c9d2019ae4fd8cdb08ff2583491d1fc39aefed7024f7390be49aa28a6b40c0842054fe1ab308eb35f8c7b7298938b639368b5b72ae24260b5e01";

    // Point R on E0 which we will evaluate
    static XR_STR: &str = "15f61615d8ef33b779cbfb12fec09a9be7f7ff7df31b7cdafc92cfa3e258ac60ba623794d659ecda7a3825e9eac8cbc8fc7b079e5beb00853a45f48ca403ca8e03d42813e1e1beb68c94df44b83bc0b270e795df42d69612914447e297993b3a60c9b084def71e5ab04cce880d01";

    // phi_1(R)
    static XR1_STR: &str = "4aa1bcc0d40b60e10e86000c1a1f33bd6c00eb16d927b60f2e97b9909085c099938881c1e37295cadd1f06f726ee384a7b6f7e9235fa001b108b7c3c0c6a142a433fc5c900e124843ea9194488096c6943eb7d20479156320ded580782437179202cf7157cb7466a4f05f01fb900";

    // phi_2(R)
    static XR2_STR: &str = "1734a4317c91be012cd4fc497b77b88adbb4b4f0121261fc381aa119c9bffc6431ce72559c9c31dc50ad8ae59c98b4986cab3b7873d90041d05b59ba616e9da9887dab7adf1b02550136a3cc54f2dca0537e163830d201a83f6633db40e8fb588741cf5eddb07af9bede0b77f501";

    #[test]
    fn test_three_isogeny_naive() {
        // Ensure all Fp2 elements from Sage decode correctly.
        let (A_test, check) = Fp2::decode(&hex::decode(A1_STR).unwrap());
        assert!(check == u32::MAX);

        let (xP, check) = Fp2::decode(&hex::decode(XP_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR, check) = Fp2::decode(&hex::decode(XR_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR_img_test, check) = Fp2::decode(&hex::decode(XR1_STR).unwrap());
        assert!(check == u32::MAX);

        // Domain, kernel and point to eval
        let six = Fp2::THREE.mul2();
        let E = Curve::new(&six);
        let ker = PointX::new(&xP, &Fp2::ONE);
        let R = PointX::new(&xR, &Fp2::ONE);

        // Values from Sage to compare against
        let E1_test = Curve::new(&A_test);
        let j1_test = E1_test.j_invariant();
        let mut R_img_test = PointX::new(&xR_img_test, &Fp2::ONE);

        // Compute chain
        let images = &mut [R];
        let E1 = three_isogeny_chain_naive(&E, &ker, 137, images);
        let R_img = images[0];

        // Assert that the codomains are isomorphic
        let j1 = E1.j_invariant();
        assert!(j1.equals(&j1_test) == u32::MAX);

        // The output from Sage might not be exactly equal, so to
        // compare points we compute an isomorphism and evaluate.
        // Compute an isomorphism between our output and Sage's
        let isomorphism = Isomorphism::new(&E1_test, &E1);
        isomorphism.isomorphism_eval(&mut R_img_test);

        // Assert the two points are equal after isomorphism.
        assert!(R_img.x().equals(&R_img_test.x()) == u32::MAX);
    }

    #[test]
    fn test_three_isogeny() {
        // Ensure all Fp2 elements from Sage decode correctly.
        let (A_test, check) = Fp2::decode(&hex::decode(A2_STR).unwrap());
        assert!(check == u32::MAX);

        let (xQ, check) = Fp2::decode(&hex::decode(XQ_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR, check) = Fp2::decode(&hex::decode(XR_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR_img_test, check) = Fp2::decode(&hex::decode(XR2_STR).unwrap());
        assert!(check == u32::MAX);

        // Domain, kernel and point to eval
        let six = Fp2::THREE.mul2();
        let E = Curve::new(&six);
        let ker = PointX::new(&xQ, &Fp2::ONE);
        let R = PointX::new(&xR, &Fp2::ONE);

        // Values from Sage to compare against
        let E2_test = Curve::new(&A_test);
        let j2_test = E2_test.j_invariant();
        let mut R_img_test = PointX::new(&xR_img_test, &Fp2::ONE);

        // Compute chain
        let images = &mut [R];
        let E2 = three_isogeny_chain(&E, &ker, 137, images);
        let R_img = images[0];

        // Assert that the codomains are isomorphic
        let j2 = E2.j_invariant();
        assert!(j2.equals(&j2_test) == u32::MAX);

        // The output from Sage might not be exactly equal, so to
        // compare points we compute an isomorphism and evaluate.
        // Compute an isomorphism between our output and Sage's
        let isomorphism = Isomorphism::new(&E2_test, &E2);
        isomorphism.isomorphism_eval(&mut R_img_test);

        // Assert the two points are equal after isomorphism.
        assert!(R_img.x().equals(&R_img_test.x()) == u32::MAX);
    }
}
