#![allow(non_snake_case)]

#[cfg(test)]
mod test_product_isogeny {
    use isogeny::{
        elliptic::{curve::Curve, projective_point::Point},
        theta::elliptic_product::{EllipticProduct, ProductPoint},
    };

    // Modulus used in test
    // https://github.com/ThetaIsogenies/two-isogenies/
    static MODULUS: [u64; 4] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xB52F88A2BBB638F2,
        0x300C882522D1C193,
    ];
    fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = MODULUS,);

    // Curve coefficients
    static A2_STR: &str = "2687f041b47d1ea9c00ae938b2a761aac6bad9ea80dd9dbe1c24d9ef697d7d0475898661998dd3a7b186e2558d1cf0dd771fb49483988c2ff578547815e8f00e";

    // Kernel Data
    static P1_X_STR: &str = "63385515142ea2f7a0c3747c22668aef99b31f43987354116da1915c24ff2f0a279051962feace976834986fc955b11bb8e1ffea47d8ce994ad22ee86c7f7a00";
    static P1_Y_STR: &str = "38243631804d307b72ecd037da591d1f06ca606bf1bb71d77ce10467b00f7b082472068bfea9eb9d80b68e04bb194a23e5214ba41625915d8e590024e5dcf611";
    static P2_X_STR: &str = "3e218d8b18cf09ce29cefaa35467225134910411fe33625136e50f7b9c59e51c517629d8786a98603cc06470a10dea83f7eca03b9b378297c21755bf0aee1324";
    static P2_Y_STR: &str = "06d56830cea82c91cf4059078566145d4b90d992177916ffe1380060057d75277610f27fc5558e5d028699493a300d84521f0c077c6e52c6adc1820c8f53f11a";
    static Q1_X_STR: &str = "3b569c320301eb5aa1aa078b7399d31e2e0e6c70e91223d1d3346be7145c100b16e4c042048452157133174122c04b1c9a17f38c28e959828933a95eebf6a305";
    static Q1_Y_STR: &str = "04a2a994555e67a3ddbb5f87dbef7903a9fc9724cb36e51924d28522222c3d2bd07a4d44b06b176b63d78733a1a4839606276b8c523bd6dc8f23de01e2817e17";
    static Q2_X_STR: &str = "4076f17bf841c5d6acc194e75a4cd020fe3ea03b0914cf3f3db9cc882f8b4724a18ce4bb13c2b3c46abd8e6cdb502dd7f48e58ee49d6c1d632532f6ed995e12b";
    static Q2_Y_STR: &str = "5b88c588b1f27496800acbef34817c5fa5cbcd728de00e31a46fc7aa6ef9af0e173d1ba96e1b2c2ebc5bc3dd3f980344b508b9df1863fb624855dc1a8cc17b12";

    // Points to push through isogeny
    static PA_X_STR: &str = "95a2c9abfff13d34db8c9d79ecdeb8cf8cbdb2de119f5fc3c51e69dcffab602e79fc55b299bdc279f8c3c20eac06322b43cb71718e367e1f795b4c8fc59fbb0b";
    static PA_Y_STR: &str = "b5b5fdfca0cc051994cbeb4338a5e629714df6cc7496cdd98900fdf281ff342990f21fa62afcf762a1c6f3e65635ab87f5c90269722434f2479482b2ec089a05";

    // Images to check against (from Sage script)
    static PA_X_IM_STR: &str = "0b2588ba81d363ee3d7c11e52e8e9db94986820a433eefb9c294de6393a1fa09eaa2baf6945de410af1ab4c189a1d5b885257ad7d2133400b03a24e4aef1f603";
    static PA_Y_IM_STR: &str = "bc5f46a11f3497cdbf9c05497ad8fa57104bdfb9431c3c327b702d164618d8201cf44b3eebc8a5d7c7214668523ec2e8fd4ad135faed66074e605abb098d1409";
    static PB_X_IM_STR: &str = "28f000f5b02a6d3e96bc1e53c3d3eb63c77d9ac8b01386da6a1ff31dff0df101103ce4b37cccfd5a9d6a67fb207e69d12317c5eb2a94a09ae156f343c6c5f30d";
    static PB_Y_IM_STR: &str = "12a4e94423744ae02c830a5188e345b9eac47823cfd7a70062995a0517137a028863d606013f8d4c42018cd7b540ad5b8a7135bbbfa3caf94019832346732010";

    #[test]
    fn test_elliptic_product_isogeny_with_torsion() {
        let A1 = Fp2::ZERO;
        let (A2, _) = Fp2::decode(&hex::decode(A2_STR).unwrap());

        let (P1_X, _) = Fp2::decode(&hex::decode(P1_X_STR).unwrap());
        let (P1_Y, _) = Fp2::decode(&hex::decode(P1_Y_STR).unwrap());
        let (P2_X, _) = Fp2::decode(&hex::decode(P2_X_STR).unwrap());
        let (P2_Y, _) = Fp2::decode(&hex::decode(P2_Y_STR).unwrap());
        let (Q1_X, _) = Fp2::decode(&hex::decode(Q1_X_STR).unwrap());
        let (Q1_Y, _) = Fp2::decode(&hex::decode(Q1_Y_STR).unwrap());
        let (Q2_X, _) = Fp2::decode(&hex::decode(Q2_X_STR).unwrap());
        let (Q2_Y, _) = Fp2::decode(&hex::decode(Q2_Y_STR).unwrap());

        let (PA_X, _) = Fp2::decode(&hex::decode(PA_X_STR).unwrap());
        let (PA_Y, _) = Fp2::decode(&hex::decode(PA_Y_STR).unwrap());

        let (im_PA_X, _) = Fp2::decode(&hex::decode(PA_X_IM_STR).unwrap());
        let (im_PA_Y, _) = Fp2::decode(&hex::decode(PA_Y_IM_STR).unwrap());
        let (im_PB_X, _) = Fp2::decode(&hex::decode(PB_X_IM_STR).unwrap());
        let (im_PB_Y, _) = Fp2::decode(&hex::decode(PB_Y_IM_STR).unwrap());

        // Curves which define elliptic product
        let E1 = Curve::new(&A1);
        let E2 = Curve::new(&A2);
        let E1E2 = EllipticProduct::new(&E1, &E2);

        // Points on E1 x E2
        let P1 = Point::new_xy(&P1_X, &P1_Y);
        let P2 = Point::new_xy(&P2_X, &P2_Y);
        let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
        let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
        let P1P2 = ProductPoint::new(&P1, &P2);
        let Q1Q2 = ProductPoint::new(&Q1, &Q2);

        // Point to push through isogeny
        let PA = Point::INFINITY;
        let PB = Point::new_xy(&PA_X, &PA_Y);
        let PAPB = ProductPoint::new(&PA, &PB);
        let image_points = [PAPB];

        // Points to compare against
        let im_PA = Point::new_xy(&im_PA_X, &im_PA_Y);
        let im_PB = Point::new_xy(&im_PB_X, &im_PB_Y);

        // Length of isogeny chain
        let n = 126;

        // Compute chain
        let (E3E4, images, ok) =
            E1E2.elliptic_product_isogeny_with_torsion(&P1P2, &Q1Q2, n, &image_points, true);
        assert!(ok == u32::MAX);

        let (_, E4) = E3E4.curves();
        assert!(E4.A.is_zero() == u32::MAX);

        let (P1, P2) = images[0].points();
        let P1_check: bool = P1.equals(&im_PA) | P1.equals(&(-im_PA)) == u32::MAX;
        let P2_check: bool = P2.equals(&im_PB) | P2.equals(&(-im_PB)) == u32::MAX;

        assert!(P1_check);
        assert!(P2_check);
    }

    #[test]
    fn test_elliptic_product_isogeny() {
        let A1 = Fp2::ZERO;
        let (A2, _) = Fp2::decode(&hex::decode(A2_STR).unwrap());
        let (P1_X, _) = Fp2::decode(&hex::decode(P1_X_STR).unwrap());
        let (P1_Y, _) = Fp2::decode(&hex::decode(P1_Y_STR).unwrap());
        let (P2_X, _) = Fp2::decode(&hex::decode(P2_X_STR).unwrap());
        let (P2_Y, _) = Fp2::decode(&hex::decode(P2_Y_STR).unwrap());
        let (Q1_X, _) = Fp2::decode(&hex::decode(Q1_X_STR).unwrap());
        let (Q1_Y, _) = Fp2::decode(&hex::decode(Q1_Y_STR).unwrap());
        let (Q2_X, _) = Fp2::decode(&hex::decode(Q2_X_STR).unwrap());
        let (Q2_Y, _) = Fp2::decode(&hex::decode(Q2_Y_STR).unwrap());
        let (PA_X, _) = Fp2::decode(&hex::decode(PA_X_STR).unwrap());
        let (PA_Y, _) = Fp2::decode(&hex::decode(PA_Y_STR).unwrap());

        let E1 = Curve::new(&A1);
        let E2 = Curve::new(&A2);
        let E1E2 = EllipticProduct::new(&E1, &E2);

        let P1 = Point::new_xy(&P1_X, &P1_Y);
        let P2 = Point::new_xy(&P2_X, &P2_Y);
        let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
        let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
        let P1P2 = ProductPoint::new(&P1, &P2);
        let Q1Q2 = ProductPoint::new(&Q1, &Q2);

        let PA = Point::INFINITY;
        let PB = Point::new_xy(&PA_X, &PA_Y);
        let PAPB = ProductPoint::new(&PA, &PB);
        let image_points = [PAPB];

        let n = 126;

        // Compute the original chain using elliptic_product_isogeny_with_torsion
        let (E3E4_orig, _, ok_orig) = E1E2.elliptic_product_isogeny_with_torsion(
            &P1P2,
            &Q1Q2,
            n,
            &image_points,
            false,
        );

        // Adjust the kernel by doubling twice to drop the 8-torsion padding
        let P1P2_strict = E1E2.double_iter(&P1P2, 2);
        let Q1Q2_strict = E1E2.double_iter(&Q1Q2, 2);

        // Compute the chain without 8-torsion
        let (E3E4_strict, _, ok_strict) = E1E2.elliptic_product_isogeny(
            &P1P2_strict,
            &Q1Q2_strict,
            n,
            &image_points,
            false,
        );

        assert_eq!(ok_orig, u32::MAX, "original elliptic_product_isogeny_with_torsion failed");
        assert_eq!(ok_strict, u32::MAX, "elliptic_product_isogeny failed");

        // Compare J-Invariants of codomain against elliptic_product_isogeny_with_torsion
        let (E3_orig, E4_orig) = E3E4_orig.curves();
        let (E3_strict, E4_strict) = E3E4_strict.curves();

        let j_E3_orig = &E3_orig.j_invariant();
        let j_E4_orig = &E4_orig.j_invariant();
        let j_E3_strict = &E3_strict.j_invariant();
        let j_E4_strict = &E4_strict.j_invariant();

        let strict_matches_orig = (j_E3_orig.equals(&j_E3_strict) & j_E4_orig.equals(&j_E4_strict)) |
                                  (j_E3_orig.equals(&j_E4_strict) & j_E4_orig.equals(&j_E3_strict));

        assert_eq!(
            strict_matches_orig,
            u32::MAX,
            "codomain mismatch between elliptic_product_isogeny_with_torsion and elliptic_product_isogeny"
        );
    }
}

#[cfg(test)]
mod test_product_isogeny_castryck_decru {
    use isogeny::{
        elliptic::{curve::Curve, projective_point::Point},
        theta::elliptic_product::{EllipticProduct, ProductPoint},
    };

    static MODULUS: [u64; 4] = [
        0xffffffffffffffff,
        0xdcdfffffffffffff,
        0xa0e11417aba2a421,
        0x00004937aae3413f,
    ];
    fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = MODULUS,);

    const A: u64 = 90503613537036670;
    const B: u64 = 300628336613559557;

    fn iota(p: &Point<Fp2>) -> Point<Fp2> {
        let mut res = *p;
        res.X = -res.X;
        res.Y *= Fp2::ZETA;
        res
    }

    fn psi_action(curve: &Curve<Fp2>, p: &Point<Fp2>) -> Point<Fp2> {
        let a_bytes = A.to_le_bytes();
        let b_bytes = B.to_le_bytes();
        let a_p = curve.mul(p, &a_bytes, 64);
        let b_p = curve.mul(p, &b_bytes, 64);
        let iota_b_p = iota(&b_p);
        curve.add(&a_p, &iota_b_p)
    }

    fn decode_fp2(x0_hex: &str, x1_hex: &str) -> Fp2 {
        let (x0, _) = Fp::decode(&hex::decode(x0_hex).unwrap());
        let (x1, _) = Fp::decode(&hex::decode(x1_hex).unwrap());
        Fp2::new(&x0, &x1)
    }

    fn decode_fp2_flat(h: &str) -> Fp2 {
        let bytes = hex::decode(h).unwrap();
        let (x0, _) = Fp::decode(&bytes[..30]);
        let (x1, _) = Fp::decode(&bytes[30..]);
        Fp2::new(&x0, &x1)
    }

    #[test]
    fn test_elliptic_product_isogeny_point_eval() {
        let E0 = Curve::<Fp2> {
            A: Fp2::ZERO,
            A24: Fp2::new(&Fp::new([0, 0, 0, 9223372036854775808]), &Fp::ZERO),
        };

        let EB = Curve::<Fp2> {
            A: decode_fp2(
                "5a32b8fee4cb360d1934edd47b54a665e2ff853724d2c6e0467eb01dae0f",
                "633d6ec480a858b3ffb3c09d369fa82f5b79a343bbf7ab0d8a726a841e46"
            ),
            A24: decode_fp2(
                "970cae3ff9b24d43064d3bf51e956999f87fe10d89b431b8911f6c87eb03",
                "588f1b31202ad6ecff2c70a7cd27923170d9a29100cd133c920d45617148"
            ),
        };

        let P2 = Point::<Fp2> {
            X: decode_fp2(
                "5e5c6e8e8aafe8a017472bf00902dd634a405014ea027ad57bc011182106",
                "7566814d948e0f527f8229300052560cc06d3dae298d32096c921e790806"
            ),
            Y: decode_fp2(
                "0d81d61a4a834a5c9edba4de8cf64a2975efd16db778cdffce2aaa375a27",
                "5d463d3eab5e0f81bde9aace44c5224ea5ab5da6a46284ebd2348f5e7c0f"
            ),
            Z: Fp2::ONE,
        };

        let Q2 = Point::<Fp2> {
            X: decode_fp2(
                "bf32e2f5f90e97a5f98b418836961d92ffd4223fe630af5d1230c0270a32",
                "5d2aa79a442cce1ec1fc913f00112a1395a9f4dd8afefb0decae0446f91d"
            ),
            Y: decode_fp2(
                "3eab3190e868a399ffcfc51303fdf0ef162b077f114e39df05505511c70b",
                "c9404ff59d64dfbc0ea045d05a04b9bede3abd81de29764ccad0f7434420"
            ),
            Z: Fp2::ONE,
        };

        let PB = Point::<Fp2> {
            X: decode_fp2(
                "ca88750e54a3167c2a778cb2d02dbcb9b3ce617e306d61c6da2385e50517",
                "b23bd2393cf7746e443fe8f727bed37e1fd26d51738129fa0943ad160f2d"
            ),
            Y: decode_fp2(
                "4e6204f764d176f69c5be3c073b43c0495c86d5c1dd9cb909acc671d3a23",
                "5b24e713f71ece6a3f5778cb7eb3d31fba4dbfae06d914e65339a4cb4928"
            ),
            Z: Fp2::ONE,
        };

        let QB = Point::<Fp2> {
            X: decode_fp2(
                "637264529c09212772bcc03e55b5cf512f80e00c4b24d998432db4ac4848",
                "42857cec556b0777f4a6daa2554e01be6900d5b6d7d7b95775238466263e"
            ),
            Y: decode_fp2(
                "d53b250cc535a206f737dba16a29bb1417b166e2589a8dfdd84ffb891f45",
                "d036af8478a99652a230ba4700d9cb8a8047a3b693d5fa0d18d2c50fc010"
            ),
            Z: Fp2::ONE,
        };

        let P3 = Point::<Fp2> {
            X: decode_fp2(
                "db53810117e7bfe64139ae22521d9e1b6f3b8e2f1c39f467c5d43803f207",
                "6f6fc5872440e6c89a43859944e7e050f8909ccb69c8dee98cd537611429"
            ),
            Y: decode_fp2(
                "6bc75a9ec1727c309ef14e23ca5e64296fa022f18f5d119d53025c068d2c",
                "07be41c87e42cefe2945987b55707bd0f5577b85ab021a69981c3767e031"
            ),
            Z: Fp2::ONE,
        };

        let Q3 = Point::<Fp2> {
            X: decode_fp2(
                "4faf804cebf8e45516ac63fefcf5f274dfd9544699d2c654e1e50e70f30f",
                "242c02924b6dffef81075bd5374154c3a1a63634a1e6fbbf745d44f3ab18"
            ),
            Y: decode_fp2(
                "88cad573d1bd074e1b2309757959c5812f265faa73a6e36f6704b078df23",
                "48e6d9bafa0da2fe3b726e71853655db0e2b4c588e889dc91726bd396d2c"
            ),
            Z: Fp2::ONE,
        };

        let psi_p2 = psi_action(&E0, &P2);
        let psi_q2 = psi_action(&E0, &Q2);

        let eb_e0 = EllipticProduct::new(&EB, &E0);

        let k1 = ProductPoint::new(&PB, &psi_p2);
        let k2 = ProductPoint::new(&QB, &psi_q2);

        let psi_p3 = psi_action(&E0, &P3);
        let psi_q3 = psi_action(&E0, &Q3);

        let eval_pts = [
            ProductPoint::new(&Point::INFINITY, &psi_p3),
            ProductPoint::new(&Point::INFINITY, &psi_q3),
        ];

        // Evaluate chain with sqrt strategy
        let (e_out, images, ok) = eb_e0.elliptic_product_isogeny(&k1, &k2, 117, &eval_pts, false);
        assert_eq!(ok, u32::MAX, "Isogeny chain evaluation failed");

        let (e3, e4) = e_out.curves();

        let j3 = e3.j_invariant();
        let j_e0 = E0.j_invariant();

        // Identify which curve corresponds to EB_prime in the product structure
        let (eb_prime, p_idx) = if j3.equals(&j_e0) == u32::MAX {
            (e4, 1)
        } else {
            (e3, 0)
        };

        let (p3_e3, p3_e4) = images[0].points();
        let phi_b_p3 = if p_idx == 1 { p3_e4 } else { p3_e3 };

        let (q3_e3, q3_e4) = images[1].points();
        let phi_b_q3 = if p_idx == 1 { q3_e4 } else { q3_e3 };

        let inv_z_p = phi_b_p3.Z.invert();
        let p_x = phi_b_p3.X * inv_z_p;
        let p_y = phi_b_p3.Y * inv_z_p;

        let inv_z_q = phi_b_q3.Z.invert();
        let q_x = phi_b_q3.X * inv_z_q;
        let q_y = phi_b_q3.Y * inv_z_q;

        let expected_a_str = "550f146c857d327df9749ef9092ef9310d4d793c4de236e88a6f432c7642bbd0f92b5f64349cd0a9816064c4cdc2dc32f12514ca01e77f982cafb615";
        let expected_p3_x_str = "f388c4a644b9a217edf4a233d4dbd58b5f017bcccfb72d3c233ee6993911bf01563542b27e2d982e91c347be7f47e8d58f717f33a7b6f45386257848";
        let expected_p3_y_str = "d3545defd7a6fdd6a3396f8fafd7dddcf0aabe6c96711fb4f3b8deeab8393858e2614104bfdb780c97bc0a57de92058bb54f530f4a655eef9cf59713";
        let expected_q3_x_str = "0bb4ebc4850f7215b871ef885f1f895500504803405ee3b4059bc7142a31c5ccbb70e3de616360a950ae6aa9ccc70da8bd6926a8bd14bc186c931311";
        let expected_q3_y_str = "fbae76347b91cd8f45921560c3b9339ae5d88935b0ee335c86411dfd4f23f59d6541dc2423a51241eff2fabe8f53ddbca3a9945e9b458c3f6ee89840";

        let exp_a = decode_fp2_flat(expected_a_str);
        let exp_p3_x = decode_fp2_flat(expected_p3_x_str);
        let exp_p3_y = decode_fp2_flat(expected_p3_y_str);
        let exp_q3_x = decode_fp2_flat(expected_q3_x_str);
        let exp_q3_y = decode_fp2_flat(expected_q3_y_str);

        assert_eq!(eb_prime.A.equals(&exp_a), u32::MAX, "Codomain A parameter mismatch");
        assert_eq!(p_x.equals(&exp_p3_x), u32::MAX, "Evaluated P3 X-coordinate mismatch");
        assert_eq!(p_y.equals(&exp_p3_y), u32::MAX, "Evaluated P3 Y-coordinate mismatch");
        assert_eq!(q_x.equals(&exp_q3_x), u32::MAX, "Evaluated Q3 X-coordinate mismatch");
        assert_eq!(q_y.equals(&exp_q3_y), u32::MAX, "Evaluated Q3 Y-coordinate mismatch");
    }
}
