#![allow(non_snake_case)]

mod benchmark_scholten {
    use criterion::{Criterion, black_box, criterion_group};
    use isogeny::{
        elliptic::{curve::Curve, point::PointX, projective_point::Point},
        theta::elliptic_product::{EllipticProduct, ProductPoint},
    };
    use std::time::Duration;

    fn benchmark_dim_one_isogeny(c: &mut Criterion) {
        use fp2::traits::Fp as _;
        use isogeny::fields::sqisign::SqiField248 as Fp2;

        // Curve coefficients
        static A1_STR: &str = "b6fb75410e12fa25edb72702cf8bbdc307dfaf9f95c267d5bea6fbe7199dcb040dcfec3f4049bcb00d88513ba18c0dd305c6fb02df2da80e4044107b2c96eb00";
        static A2_STR: &str = "9048175252d6505ac864c361add9849120bbe707683ab4183c9d754f24465304f325ea7fe90c5507149769217ac2dbbdcb8869b2ae48d9ea652c001928804704";

        // Kernel Data
        static KER_X_STR: &str = "8760ae8f71373123b3f98b11f0fe332f39f58e45b5e0a7783245e0d9abed0a03651ad0d3213842667e2a416b249ca37b4ff62d0560122e9d8ff518e3b018dc03";

        let (A1, c1) = Fp2::decode(&hex::decode(A1_STR).unwrap());
        assert_eq!(c1, u32::MAX);

        let (A2, c2) = Fp2::decode(&hex::decode(A2_STR).unwrap());
        assert_eq!(c2, u32::MAX);

        let (KER_X, c3) = Fp2::decode(&hex::decode(KER_X_STR).unwrap());
        assert_eq!(c3, u32::MAX);

        // Domain curve
        let E1 = Curve::new(&A1);

        // Codomain curve
        let E2 = Curve::new(&A2);

        // Isogeny kernel
        let ker = PointX::new(&KER_X, &Fp2::ONE);

        // Length of isogeny chain
        let n = 125;

        // Compute chain
        let (E3, ok) = E1.two_isogeny_chain(&ker, n, &mut []);

        // ensure the output matches what we expect
        assert_eq!(ok, u32::MAX);
        assert_eq!(E2.j_invariant().equals(&E3.j_invariant()), u32::MAX);

        let bench_id = format!("Benchmarking 2^{n} isogeny over Fp2",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| black_box(E1).two_isogeny_chain(&black_box(ker), black_box(n), &mut []))
        });
    }

    fn benchmark_dim_two_isogeny(c: &mut Criterion) {
        use fp2::traits::Fp as _;
        use isogeny::fields::sqisign::SqiField248 as Fp2;

        // Domain Curve coefficients
        const A1_STR: &str = "baa79b0dc07508bb6fea4685db4b48f237686ad1e12964c985814261bcee97015c28136967c8faa77df1d28ffe4f81c68369091bc503d27ab21c459ce88eb101";
        const A2_STR: &str = "84228651f271b0f39f2f19f2e8718f31ed3365ac9e5cb303afe663d0cfc11f0455d891b0ca6c7e653f9ba2667730bb77befe1b1a31828404284af8fd7baacc01";

        // Codomain Curve coefficients
        const A3_STR: &str = "c9ae315618b89d82c3ba2c2f6ba312b3fd7e67caf4ed47ff22f826f775b7af0178d3a7399d5e086b1913535fe9931d94511c79715e0f9dffbba43a01059a2d00";
        const A4_STR: &str = "2b14b4831c280c523790d20739de9502f35d4d44b4d5abad52d6d15b4bc21a0342a8c6abdf9a7f2fee3be93ca5f71f425051664645eb54af4391682d650bd100";

        // Kernel Data
        const P1_X_STR: &str = "4c60e805dfacd69fc9c55768dafbd9a43081b7dfe0ffc95a80b1e6e5a19e5400d745f4962fe258b14597ed5b0b712ba86bee4e7f92b44998ef8003e6f8316d02";
        const P1_Y_STR: &str = "5edad9e9ab20b0aba7604ede57da0d70b57e03dfd5740f428dda76df58083d035f524c3e854d2be6aae9e8f788b682876920d8b813e42df2ed47eba9551ac902";
        const P2_X_STR: &str = "30c04fdde6b6c8afdb9175601ecb11d6e8d5084f2f6476f6ca7e3461bb44a30201bb5a8d83010c7b190227ecaecc6eb99fa078116543e0157d6f4d87e783a601";
        const P2_Y_STR: &str = "ae52fd0667c7557eea8c0e9810bd5137e3d031653471569b365cf9d5641d84032ae06689365f6e0e5c9a91afc145afcb1a091dc936677ddcc31eab36b1860a00";
        const Q1_X_STR: &str = "2b3d75daaa85396899b42f4cd55dc6006a0e1a32632a5ab1a3de29d9fa859501151736c5825e2471755db281d1289f8358704af1d064c2d0496ece8de6de3b03";
        const Q1_Y_STR: &str = "061cfae6a523376d4c71ce34ff5ab3eda127b9e6cf6c247c4bf12d9eb68382009025b8a5db5912a82096e6982d4c3d5b84916cec91cdfbbee074d023287b3a01";
        const Q2_X_STR: &str = "6baeb9014010f9b7b6afc0f7371c1d91dd6ce55ae87da3049ef76b9fd35e0c04e0c4358edf3b97943eabf4a17556e18437cfdcb92903c2cc9a9325899ee17002";
        const Q2_Y_STR: &str = "122065fdb3e13962329cb8a6405b29e5e3f1c0f8f72c53941705bb655d124500177864fe72d6a648f8af1f434dad1ac00fd884cf8da592a45ba7bc7fbad07604";

        let (A1, _) = Fp2::decode(&hex::decode(A1_STR).unwrap());
        let (A2, _) = Fp2::decode(&hex::decode(A2_STR).unwrap());

        let (A3, _) = Fp2::decode(&hex::decode(A3_STR).unwrap());
        let (A4, _) = Fp2::decode(&hex::decode(A4_STR).unwrap());

        let (P1_X, _) = Fp2::decode(&hex::decode(P1_X_STR).unwrap());
        let (P1_Y, _) = Fp2::decode(&hex::decode(P1_Y_STR).unwrap());
        let (P2_X, _) = Fp2::decode(&hex::decode(P2_X_STR).unwrap());
        let (P2_Y, _) = Fp2::decode(&hex::decode(P2_Y_STR).unwrap());
        let (Q1_X, _) = Fp2::decode(&hex::decode(Q1_X_STR).unwrap());
        let (Q1_Y, _) = Fp2::decode(&hex::decode(Q1_Y_STR).unwrap());
        let (Q2_X, _) = Fp2::decode(&hex::decode(Q2_X_STR).unwrap());
        let (Q2_Y, _) = Fp2::decode(&hex::decode(Q2_Y_STR).unwrap());

        // Curves which define elliptic product
        let E1 = Curve::new(&A1);
        let E2 = Curve::new(&A2);
        let E1E2 = EllipticProduct::new(&E1, &E2);

        // Codomain Curves
        let E3 = Curve::new(&A3);
        let E4 = Curve::new(&A4);

        // Points on E1 x E2
        let P1 = Point::new_xy(&P1_X, &P1_Y);
        let P2 = Point::new_xy(&P2_X, &P2_Y);
        let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
        let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
        let P1P2 = ProductPoint::new(&P1, &P2);
        let Q1Q2 = ProductPoint::new(&Q1, &Q2);

        // Length of isogeny chain
        let n = 125;

        // Compute chain
        let (F3F4, _, ok) = E1E2.elliptic_product_isogeny(&P1P2, &Q1Q2, n, &[], true);
        assert!(ok == u32::MAX);

        let (F3, F4) = F3F4.curves();
        assert_eq!(E3.j_invariant().equals(&F4.j_invariant()), u32::MAX);
        assert_eq!(E4.j_invariant().equals(&F3.j_invariant()), u32::MAX);

        let bench_id = format!("Benchmarking (2^{n}, 2^{n}) isogeny over Fp2.",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| {
                black_box(E1E2).elliptic_product_isogeny(
                    &black_box(P1P2),
                    &black_box(Q1Q2),
                    black_box(n),
                    &[],
                    false,
                )
            })
        });
    }

    criterion_group! {
        name = benchmark_product;
        config = Criterion::default().measurement_time(Duration::from_secs(10));
        targets = benchmark_dim_one_isogeny, benchmark_dim_two_isogeny
    }
}

fn main() {
    benchmark_scholten::benchmark_product();
}
