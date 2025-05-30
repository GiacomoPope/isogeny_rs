#![allow(non_snake_case)]

#[cfg(test)]
mod test_csidh_velu {
    use isogeny::{
        elliptic::{curve::Curve, point::PointX},
        polynomial_ring::poly::Polynomial,
    };

    // CSIDH 512 prime.
    const MODULUS: [u64; 8] = [
        0x1B81B90533C6C87B,
        0xC2721BF457ACA835,
        0x516730CC1F0B4F25,
        0xA7AAC6C567F35507,
        0x5AFBFCC69322C9CD,
        0xB42D083AEDC88C42,
        0xFC8AB0D15E3E4C4A,
        0x65B48E8F740F89BF,
    ];
    fp2::define_fp_core!(typename = Fp, modulus = MODULUS,);

    type P = Polynomial<Fp>;

    const P103_X_BYTES: [u8; 64] = [
        10, 101, 12, 18, 20, 58, 190, 183, 31, 106, 251, 46, 117, 208, 205, 0, 107, 11, 165, 131,
        209, 236, 94, 28, 131, 191, 90, 27, 169, 183, 177, 232, 184, 253, 81, 112, 173, 171, 3, 40,
        82, 4, 35, 58, 6, 72, 122, 207, 250, 143, 143, 158, 115, 63, 137, 255, 167, 228, 53, 34,
        239, 74, 57, 81,
    ];
    const P211_X_BYTES: [u8; 64] = [
        238, 189, 22, 253, 65, 106, 138, 139, 190, 131, 156, 59, 59, 190, 21, 233, 174, 176, 192,
        195, 153, 179, 182, 105, 153, 207, 46, 81, 87, 20, 18, 219, 59, 234, 94, 76, 42, 161, 219,
        64, 125, 18, 249, 176, 192, 106, 185, 52, 186, 120, 3, 43, 102, 34, 142, 158, 116, 209,
        128, 206, 38, 220, 39, 95,
    ];
    const P311_X_BYTES: [u8; 64] = [
        91, 231, 70, 194, 23, 155, 11, 243, 133, 79, 126, 163, 30, 195, 233, 66, 90, 73, 130, 99,
        201, 169, 156, 140, 207, 94, 104, 147, 102, 42, 34, 22, 17, 17, 125, 54, 214, 38, 252, 116,
        52, 29, 111, 43, 189, 39, 0, 97, 220, 143, 44, 183, 30, 50, 224, 4, 91, 217, 50, 188, 28,
        248, 49, 96,
    ];
    const P587_X_BYTES: [u8; 64] = [
        208, 223, 49, 61, 194, 244, 239, 9, 36, 20, 150, 189, 4, 98, 68, 49, 2, 82, 184, 250, 135,
        61, 188, 99, 169, 246, 32, 23, 64, 24, 96, 82, 187, 189, 77, 176, 221, 63, 218, 65, 29, 44,
        31, 159, 35, 25, 201, 79, 58, 228, 140, 42, 34, 42, 194, 48, 164, 47, 125, 226, 159, 124,
        215, 19,
    ];

    const A103_BYTES: [u8; 64] = [
        176, 146, 168, 142, 34, 65, 67, 205, 144, 174, 85, 98, 64, 194, 29, 197, 131, 160, 208,
        132, 137, 89, 172, 179, 151, 224, 125, 200, 254, 84, 35, 248, 104, 58, 86, 252, 40, 50,
        122, 63, 188, 178, 136, 47, 150, 0, 168, 217, 2, 40, 73, 127, 138, 116, 51, 170, 208, 74,
        245, 115, 18, 190, 69, 62,
    ];
    const A211_BYTES: [u8; 64] = [
        40, 57, 173, 8, 5, 128, 178, 175, 193, 182, 157, 167, 18, 153, 250, 112, 96, 241, 172, 104,
        243, 64, 191, 104, 102, 71, 35, 164, 60, 231, 125, 213, 223, 178, 56, 133, 62, 134, 110,
        39, 3, 42, 145, 108, 148, 32, 164, 66, 79, 238, 83, 95, 74, 107, 129, 11, 190, 74, 178,
        175, 60, 176, 24, 64,
    ];
    const A311_BYTES: [u8; 64] = [
        89, 43, 146, 254, 39, 28, 73, 108, 102, 15, 144, 157, 34, 159, 168, 205, 141, 33, 0, 161,
        47, 38, 172, 204, 249, 232, 3, 79, 20, 13, 223, 160, 6, 158, 83, 181, 132, 161, 250, 131,
        171, 44, 4, 59, 150, 88, 131, 142, 156, 114, 77, 53, 244, 28, 62, 23, 59, 55, 56, 23, 63,
        103, 151, 8,
    ];
    const A587_BYTES: [u8; 64] = [
        99, 164, 168, 164, 123, 19, 25, 132, 44, 91, 235, 107, 139, 228, 68, 154, 5, 32, 226, 199,
        207, 162, 164, 67, 6, 236, 167, 158, 121, 221, 59, 182, 25, 113, 68, 137, 43, 193, 177,
        154, 93, 238, 25, 71, 120, 131, 205, 202, 105, 110, 85, 248, 120, 170, 49, 163, 112, 192,
        163, 235, 212, 111, 68, 35,
    ];

    const P103_X: Fp = Fp::const_decode_no_check(&P103_X_BYTES);
    const P211_X: Fp = Fp::const_decode_no_check(&P211_X_BYTES);
    const P311_X: Fp = Fp::const_decode_no_check(&P311_X_BYTES);
    const P587_X: Fp = Fp::const_decode_no_check(&P587_X_BYTES);

    const CODOMAIN_103: Fp = Fp::const_decode_no_check(&A103_BYTES);
    const CODOMAIN_211: Fp = Fp::const_decode_no_check(&A211_BYTES);
    const CODOMAIN_311: Fp = Fp::const_decode_no_check(&A311_BYTES);
    const CODOMAIN_587: Fp = Fp::const_decode_no_check(&A587_BYTES);

    #[test]
    fn sqrt_velu_103() {
        let ker = PointX::new(&P103_X, &Fp::ONE);
        let codomain_test = Curve::new(&CODOMAIN_103);
        let mut images = [ker];

        let mut A24 = Fp::TWO;
        let mut C24 = Fp::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 103, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn sqrt_velu_211() {
        let ker = PointX::new(&P211_X, &Fp::ONE);
        let codomain_test = Curve::new(&CODOMAIN_211);
        let mut images = [ker];

        let mut A24 = Fp::TWO;
        let mut C24 = Fp::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 211, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn sqrt_velu_311() {
        let ker = PointX::new(&P311_X, &Fp::ONE);
        let codomain_test = Curve::new(&CODOMAIN_311);
        let mut images = [ker];

        let mut A24 = Fp::TWO;
        let mut C24 = Fp::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 311, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }

    #[test]
    fn sqrt_velu_587() {
        let ker = PointX::new(&P587_X, &Fp::ONE);
        let codomain_test = Curve::new(&CODOMAIN_587);
        let mut images = [ker];

        let mut A24 = Fp::TWO;
        let mut C24 = Fp::FOUR;

        Curve::sqrt_velu_odd_isogeny_proj::<P>(&mut A24, &mut C24, &ker, 587, &mut images);
        let codomain = Curve::curve_from_A24_proj(&A24, &C24);

        // Ensure the codomain matches the expected result.
        assert!(codomain.j_invariant().equals(&codomain_test.j_invariant()) == u32::MAX);

        // Pushing the kernel through should give the point at infinity.
        assert!(images[0].Z.is_zero() == u32::MAX);
    }
}
