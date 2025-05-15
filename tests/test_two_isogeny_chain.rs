#![allow(non_snake_case)]

#[cfg(test)]
mod test_two_isogeny_chain {
    use fp2::fq::Fq;
    use isogeny::elliptic::curve::Curve;
    use isogeny::elliptic::isomorphism::Isomorphism;
    use isogeny::elliptic::point::PointX;
    use isogeny::fields::sqisign::SqiSignI as Fp2;

    // Domain of isogenies for testing
    static A0_STR: &str = "8be929811c642e0148f0b6fed8f00317116795d649b0656cafb510c8ab0b450283cf1721e68b40e818039ce4c5d8d7f6fd89f4a6e35acffab2568eba06827802";

    // Kernel for phi_1 : E -> E/<P> with P not above (0 : 1)
    static XP_STR: &str = "b0aa03076df1724a4b4a733594e540e9f28db6068f7a2bdffccc05ade7220f00836b8014c5214075eba997a77dff02ad36d27aebf210a4c92c0bba382918cd04";

    // Kernel for phi_2 : E -> E/<P> with Q above (0 : 1)
    static XQ_STR: &str = "f796a5fbc2673ba9535d6e9636057ebe0a1c6c5d28bf493dd2a8bfb016d56f035676ed68cdb2ede1fc3e597e94a4dc1106c87f5b507f8a6a4368e5bb35d70e02";

    // Codomain for phi_1 : E -> E/<P>
    static A1_STR: &str = "b5cb96db586b64e7dec2170e27671f4c0ca78cb09739a5b30cd4d36d19cdb9012b8085cba5d50b072fc6591433e149a12e191dc4940ff5bb97c3c6c5af877803";

    // Codomain for phi_2 : E -> E/<Q>
    static A2_STR: &str = "a034ad0ad67ed880d562df615c5efaef6b4218fd3cc3cb3f9ecdb1e351a76401e37494cb944612d38f1f438fde68d5ca6a905dd940db5dbd157616a8b5c29902";

    // Point R on E0 which we will evaluate
    static XR_STR: &str = "6b8d10ed63062eccb819cbf3eae88f36037e250221610d3bf4a117a61177e7043b3bbaaf18e3ddcd9a109fec07598f59eff47754ef5a0b69ac236345cc88b001";

    // phi_1(R)
    static XR1_STR: &str = "d8d21cbc0783c185d77f771236cbef17f53cc2a84d0d7fcca88eec7c820bb8035b9a3cba6eee6f94be7981bfaaaa0467db675b4a88f9dea1ba3ab10c85595e01";

    // phi_2(R)
    static XR2_STR: &str = "d28e10f44aaacd2a8f4da6cb12e7be05513c6d144b16b145e78549c0cb4717049013072d85765c8661d9a1d9d3f6e24a06330ef1ec1a9065bcd9c48aeac2ae04";

    #[test]
    fn test_two_isogeny_naive() {
        // Ensure all Fp2 elements from Sage decode correctly.
        let (A, check) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(check == u32::MAX);

        let (A_test, check) = Fp2::decode(&hex::decode(A1_STR).unwrap());
        assert!(check == u32::MAX);

        let (xP, check) = Fp2::decode(&hex::decode(XP_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR, check) = Fp2::decode(&hex::decode(XR_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR_img_test, check) = Fp2::decode(&hex::decode(XR1_STR).unwrap());
        assert!(check == u32::MAX);

        // Domain, kernel and point to eval
        let E = Curve::new(&A);
        let ker = PointX::new(&xP, &Fp2::ONE);
        let R = PointX::new(&xR, &Fp2::ONE);

        // Values from Sage to compare against
        let E1_test = Curve::new(&A_test);
        let j1_test = E1_test.j_invariant();
        let mut R_img_test = PointX::new(&xR_img_test, &Fp2::ONE);

        // Compute chain
        let images = &mut [R];
        let E1 = E.two_isogeny_chain_naive(&ker, 248, images);
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
    fn test_two_isogeny_singular_naive() {
        // Ensure all Fp2 elements from Sage decode correctly.
        let (A, check) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(check == u32::MAX);

        let (A_test, check) = Fp2::decode(&hex::decode(A2_STR).unwrap());
        assert!(check == u32::MAX);

        let (xQ, check) = Fp2::decode(&hex::decode(XQ_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR, check) = Fp2::decode(&hex::decode(XR_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR_img_test, check) = Fp2::decode(&hex::decode(XR2_STR).unwrap());
        assert!(check == u32::MAX);

        // Domain, kernel and point to eval
        let E = Curve::new(&A);
        let ker = PointX::new(&xQ, &Fp2::ONE);
        let R = PointX::new(&xR, &Fp2::ONE);

        // Values from Sage to compare against
        let E2_test = Curve::new(&A_test);
        let j2_test = E2_test.j_invariant();
        let mut R_img_test = PointX::new(&xR_img_test, &Fp2::ONE);

        // Compute chain
        let images = &mut [R];
        let E2 = E.two_isogeny_chain_naive(&ker, 248, images);
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

    #[test]
    fn test_two_isogeny() {
        // Ensure all Fp2 elements from Sage decode correctly.
        let (A, check) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(check == u32::MAX);

        let (A_test, check) = Fp2::decode(&hex::decode(A1_STR).unwrap());
        assert!(check == u32::MAX);

        let (xP, check) = Fp2::decode(&hex::decode(XP_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR, check) = Fp2::decode(&hex::decode(XR_STR).unwrap());
        assert!(check == u32::MAX);

        let (xR_img_test, check) = Fp2::decode(&hex::decode(XR1_STR).unwrap());
        assert!(check == u32::MAX);

        // Domain, kernel and point to eval
        let E = Curve::new(&A);
        let ker = PointX::new(&xP, &Fp2::ONE);
        let R = PointX::new(&xR, &Fp2::ONE);

        // Values from Sage to compare against
        let E1_test = Curve::new(&A_test);
        let j1_test = E1_test.j_invariant();
        let mut R_img_test = PointX::new(&xR_img_test, &Fp2::ONE);

        // Compute chain
        let images = &mut [R];
        let E1 = E.two_isogeny_chain(&ker, 248, images);
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
}
