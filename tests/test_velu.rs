#![allow(non_snake_case)]

#[cfg(test)]
mod velu_test {
    use isogeny::{
        elliptic::{curve::Curve, point::PointX},
        polynomial_ring::poly::Polynomial,
    };

    // Toy prime, with smooth p + 1
    // p + 1 = 2^3 * 3^2 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 59 * 61 * 67 * 71 * 73 * 79 * 83 * 89 * 97 * 163^16
    const MODULUS: [u64; 4] = [
        0x65396406224d1bc7,
        0x60892bae5986e2f3,
        0x8f2e6d4fb53f8cc7,
        0x0003e368df883c71,
    ];
    fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = MODULUS,);

    fp2::define_fp_tests!(Fp);
    fp2::define_fp2_tests!(Fp2, MODULUS, 10);

    // We need the poly type for calling SqrtVelu
    type P = Polynomial<Fp2>;

    // Kernel of order 163
    const P_X_RE_BYTES: [u8; 31] = [
        161, 11, 202, 60, 238, 183, 30, 198, 113, 154, 4, 46, 251, 33, 116, 68, 50, 187, 152, 163,
        215, 172, 206, 196, 15, 96, 128, 81, 141, 27, 2,
    ];
    const P_X_IM_BYTES: [u8; 31] = [
        244, 36, 105, 80, 212, 91, 28, 64, 95, 113, 78, 111, 165, 12, 140, 76, 43, 91, 43, 245, 47,
        7, 52, 200, 71, 10, 182, 128, 100, 211, 1,
    ];
    const AP_RE_BYTES: [u8; 31] = [
        81, 48, 37, 31, 25, 41, 102, 139, 211, 247, 169, 227, 160, 80, 111, 30, 149, 218, 159, 209,
        161, 29, 168, 254, 63, 197, 31, 62, 144, 178, 0,
    ];
    const AP_IM_BYTES: [u8; 31] = [
        124, 122, 247, 73, 156, 36, 71, 39, 84, 47, 89, 103, 72, 193, 196, 240, 28, 139, 178, 131,
        241, 80, 164, 147, 32, 114, 247, 55, 117, 48, 3,
    ];

    const PX: Fp2 = Fp2::const_decode_no_check(&P_X_RE_BYTES, &P_X_IM_BYTES);
    const CODOMAIN_AP: Fp2 = Fp2::const_decode_no_check(&AP_RE_BYTES, &AP_IM_BYTES);

    // Kernel of order 163^16
    const Q_X_RE_BYTES: [u8; 31] = [
        118, 241, 155, 23, 243, 208, 140, 252, 174, 149, 216, 233, 237, 74, 71, 198, 246, 209, 86,
        117, 101, 162, 70, 220, 239, 163, 232, 36, 184, 7, 2,
    ];
    const Q_X_IM_BYTES: [u8; 31] = [
        239, 37, 95, 109, 173, 69, 148, 192, 159, 248, 224, 72, 52, 57, 70, 90, 155, 175, 69, 221,
        163, 161, 224, 164, 118, 187, 148, 3, 81, 15, 0,
    ];
    const AQ_RE_BYTES: [u8; 31] = [
        251, 125, 75, 0, 66, 82, 0, 155, 19, 210, 252, 56, 51, 126, 166, 77, 100, 207, 193, 115,
        187, 176, 243, 68, 252, 223, 190, 40, 253, 120, 1,
    ];
    const AQ_IM_BYTES: [u8; 31] = [
        39, 78, 118, 239, 124, 92, 224, 130, 201, 147, 114, 202, 203, 42, 135, 47, 20, 79, 150,
        211, 204, 221, 220, 117, 136, 62, 254, 171, 107, 98, 3,
    ];

    const QX: Fp2 = Fp2::const_decode_no_check(&Q_X_RE_BYTES, &Q_X_IM_BYTES);
    const CODOMAIN_AQ: Fp2 = Fp2::const_decode_no_check(&AQ_RE_BYTES, &AQ_IM_BYTES);

    // Kernel of order 8
    const P2_X_RE_BYTES: [u8; 31] = [
        66, 10, 195, 200, 119, 121, 247, 215, 37, 201, 174, 64, 203, 82, 174, 148, 208, 66, 52, 88,
        14, 235, 78, 87, 157, 225, 72, 132, 184, 236, 0,
    ];
    const P2_X_IM_BYTES: [u8; 31] = [
        5, 196, 47, 123, 38, 72, 69, 73, 24, 117, 23, 52, 64, 212, 176, 205, 84, 65, 190, 52, 185,
        241, 164, 176, 7, 216, 51, 17, 69, 78, 1,
    ];
    const A2_RE_BYTES: [u8; 31] = [
        58, 161, 125, 74, 214, 139, 30, 71, 195, 170, 165, 52, 16, 69, 78, 244, 230, 4, 53, 208,
        29, 146, 222, 112, 5, 38, 150, 181, 152, 26, 2,
    ];
    const A2_IM_BYTES: [u8; 31] = [
        183, 224, 108, 169, 150, 90, 167, 84, 236, 71, 3, 171, 108, 54, 91, 90, 108, 250, 7, 48,
        125, 18, 89, 228, 223, 203, 17, 68, 87, 24, 0,
    ];
    const P2X: Fp2 = Fp2::const_decode_no_check(&P2_X_RE_BYTES, &P2_X_IM_BYTES);
    const CODOMAIN_A2: Fp2 = Fp2::const_decode_no_check(&A2_RE_BYTES, &A2_IM_BYTES);

    // factorisation of (p + 1)
    const DEGREE_FACTORS: [(usize, usize); 26] = [
        (2, 3),
        (3, 2),
        (5, 1),
        (7, 1),
        (11, 1),
        (13, 1),
        (17, 1),
        (19, 1),
        (23, 1),
        (29, 1),
        (31, 1),
        (37, 1),
        (41, 1),
        (43, 1),
        (47, 1),
        (53, 1),
        (59, 1),
        (61, 1),
        (67, 1),
        (71, 1),
        (73, 1),
        (79, 1),
        (83, 1),
        (89, 1),
        (97, 1),
        (163, 16),
    ];
    const K_ODD_X_RE_BYTES: [u8; 31] = [
        110, 242, 212, 112, 152, 164, 46, 167, 149, 200, 141, 65, 17, 122, 160, 245, 223, 235, 182,
        4, 122, 3, 169, 43, 226, 43, 208, 39, 105, 187, 2,
    ];
    const K_ODD_X_IM_BYTES: [u8; 31] = [
        198, 5, 186, 228, 235, 131, 86, 59, 104, 71, 222, 191, 80, 75, 254, 246, 102, 181, 120,
        101, 156, 87, 191, 182, 86, 193, 211, 92, 123, 165, 1,
    ];
    const A_ODD_RE_BYTES: [u8; 31] = [
        254, 234, 195, 162, 121, 115, 45, 207, 89, 29, 243, 8, 228, 194, 29, 191, 134, 245, 124,
        110, 56, 196, 236, 199, 189, 251, 96, 130, 103, 25, 2,
    ];
    const A_ODD_IM_BYTES: [u8; 31] = [
        218, 193, 183, 0, 14, 153, 127, 4, 105, 173, 88, 85, 208, 249, 68, 233, 25, 137, 64, 0,
        219, 130, 147, 51, 43, 89, 195, 95, 154, 196, 3,
    ];
    const K_ODD_X: Fp2 = Fp2::const_decode_no_check(&K_ODD_X_RE_BYTES, &K_ODD_X_IM_BYTES);
    const CODOMAIN_A_ODD: Fp2 = Fp2::const_decode_no_check(&A_ODD_RE_BYTES, &A_ODD_IM_BYTES);

    // Kernel of order (p + 1)
    const K_X_RE_BYTES: [u8; 31] = [
        4, 33, 54, 182, 226, 189, 216, 94, 198, 117, 116, 30, 169, 108, 115, 77, 57, 121, 7, 20,
        95, 77, 67, 42, 153, 6, 30, 107, 133, 230, 1,
    ];
    const K_X_IM_BYTES: [u8; 31] = [
        73, 62, 50, 218, 94, 160, 254, 125, 191, 162, 21, 222, 217, 124, 198, 88, 91, 75, 111, 71,
        70, 30, 208, 212, 108, 107, 230, 218, 56, 164, 1,
    ];
    const AK_RE_BYTES: [u8; 31] = [
        185, 79, 100, 141, 25, 254, 25, 71, 197, 162, 74, 250, 231, 208, 77, 153, 138, 31, 139,
        192, 160, 30, 229, 160, 33, 14, 128, 74, 72, 151, 1,
    ];
    const AK_IM_BYTES: [u8; 31] = [
        30, 47, 108, 97, 72, 75, 108, 194, 147, 111, 108, 236, 242, 151, 159, 225, 54, 114, 24, 40,
        246, 126, 168, 12, 215, 204, 193, 70, 90, 243, 2,
    ];

    const KX: Fp2 = Fp2::const_decode_no_check(&K_X_RE_BYTES, &K_X_IM_BYTES);
    const CODOMAIN_AK: Fp2 = Fp2::const_decode_no_check(&AK_RE_BYTES, &AK_IM_BYTES);

    #[test]
    fn test_odd_prime_velu_codomain() {
        let E = Curve::new(&Fp2::ZERO);
        let ker = PointX::new(&PX, &Fp2::ONE);

        let codomain_test = Curve::new(&CODOMAIN_AP);
        let mut images = [ker];
        let codomain = E.velu_prime_isogeny::<P>(&ker, 163, &mut images);

        // Ensure the codomain matches the expected result.
        assert!(codomain.A.equals(&codomain_test.A) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn test_two_prime_power_velu_codomain() {
        let E = Curve::new(&Fp2::ZERO);
        let ker = PointX::new(&P2X, &Fp2::ONE);

        let codomain_test = Curve::new(&CODOMAIN_A2);
        let mut images = [ker];
        let codomain = E.velu_prime_power_isogeny::<P>(&ker, 2, 3, &mut images);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn test_odd_prime_power_velu_codomain() {
        let E = Curve::new(&Fp2::ZERO);
        let ker = PointX::new(&QX, &Fp2::ONE);

        let codomain_test = Curve::new(&CODOMAIN_AQ);
        let mut images = [ker];
        let codomain = E.velu_prime_power_isogeny::<P>(&ker, 163, 16, &mut images);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn test_odd_composite_velu_codomain() {
        let E = Curve::new(&Fp2::ZERO);
        let ker = PointX::new(&K_ODD_X, &Fp2::ONE);

        let codomain_test = Curve::new(&CODOMAIN_A_ODD);
        let mut images = [ker];
        let codomain = E.velu_composite_isogeny::<P>(&ker, &DEGREE_FACTORS[1..], &mut images);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn test_composite_velu_codomain() {
        let E = Curve::new(&Fp2::ZERO);
        let ker = PointX::new(&KX, &Fp2::ONE);

        let codomain_test = Curve::new(&CODOMAIN_AK);
        let mut images = [ker];
        let codomain = E.velu_composite_isogeny::<P>(&ker, &DEGREE_FACTORS, &mut images);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn test_velu_projective_codomain() {
        let ker = PointX::new(&PX, &Fp2::ONE);
        let codomain_test = Curve::new(&CODOMAIN_AP);
        let mut images = [ker];

        let mut A24 = Fp2::TWO;
        let mut C24 = Fp2::FOUR;
        Curve::velu_odd_isogeny_proj(&mut A24, &mut C24, &ker, 163, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn test_sqrt_velu_projective_codomain() {
        let ker = PointX::new(&PX, &Fp2::ONE);
        let codomain_test = Curve::new(&CODOMAIN_AP);
        let mut images = [ker];

        let mut A24 = Fp2::TWO;
        let mut C24 = Fp2::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 163, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }
}
