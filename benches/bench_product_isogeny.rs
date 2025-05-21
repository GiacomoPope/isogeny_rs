#![allow(non_snake_case)]

mod benchmark_product {
    use criterion::{Criterion, black_box, criterion_group};
    use isogeny::{
        elliptic::{curve::Curve, projective_point::Point},
        theta::elliptic_product::{EllipticProduct, ProductPoint},
    };
    use std::time::Duration;

    fn benchmark_product_isogeny(c: &mut Criterion) {
        // Modulus used in test
        // https://github.com/ThetaIsogenies/two-isogenies/
        const MODULUS: [u64; 4] = [
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xB52F88A2BBB638F2,
            0x300C882522D1C193,
        ];
        fp2::define_fp2_from_modulus!(typename = Fp2, base_typename = Fp, modulus = MODULUS,);

        // Curve coefficients
        const A2_STR: &str = "2687f041b47d1ea9c00ae938b2a761aac6bad9ea80dd9dbe1c24d9ef697d7d0475898661998dd3a7b186e2558d1cf0dd771fb49483988c2ff578547815e8f00e";

        // Kernel Data
        const P1_X_STR: &str = "63385515142ea2f7a0c3747c22668aef99b31f43987354116da1915c24ff2f0a279051962feace976834986fc955b11bb8e1ffea47d8ce994ad22ee86c7f7a00";
        const P1_Y_STR: &str = "38243631804d307b72ecd037da591d1f06ca606bf1bb71d77ce10467b00f7b082472068bfea9eb9d80b68e04bb194a23e5214ba41625915d8e590024e5dcf611";
        const P2_X_STR: &str = "3e218d8b18cf09ce29cefaa35467225134910411fe33625136e50f7b9c59e51c517629d8786a98603cc06470a10dea83f7eca03b9b378297c21755bf0aee1324";
        const P2_Y_STR: &str = "06d56830cea82c91cf4059078566145d4b90d992177916ffe1380060057d75277610f27fc5558e5d028699493a300d84521f0c077c6e52c6adc1820c8f53f11a";
        const Q1_X_STR: &str = "3b569c320301eb5aa1aa078b7399d31e2e0e6c70e91223d1d3346be7145c100b16e4c042048452157133174122c04b1c9a17f38c28e959828933a95eebf6a305";
        const Q1_Y_STR: &str = "04a2a994555e67a3ddbb5f87dbef7903a9fc9724cb36e51924d28522222c3d2bd07a4d44b06b176b63d78733a1a4839606276b8c523bd6dc8f23de01e2817e17";
        const Q2_X_STR: &str = "4076f17bf841c5d6acc194e75a4cd020fe3ea03b0914cf3f3db9cc882f8b4724a18ce4bb13c2b3c46abd8e6cdb502dd7f48e58ee49d6c1d632532f6ed995e12b";
        const Q2_Y_STR: &str = "5b88c588b1f27496800acbef34817c5fa5cbcd728de00e31a46fc7aa6ef9af0e173d1ba96e1b2c2ebc5bc3dd3f980344b508b9df1863fb624855dc1a8cc17b12";

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

        // Length of isogeny chain
        let n = 126;

        // Compute chain
        let (E3E4, _, ok) = E1E2.elliptic_product_isogeny(&P1P2, &Q1Q2, n, &[]);
        assert!(ok == u32::MAX);

        let (_, E4) = E3E4.curves();
        assert!(E4.A.is_zero() == u32::MAX);

        let bench_id = format!("Benchmarking (2^n, 2^n) isogeny from two, two repo.",);
        c.bench_function(&bench_id, |b| {
            b.iter(|| {
                black_box(E1E2).elliptic_product_isogeny(
                    &black_box(P1P2),
                    &black_box(Q1Q2),
                    black_box(n),
                    &[],
                )
            })
        });
    }

    criterion_group! {
        name = benchmark_product;
        config = Criterion::default().measurement_time(Duration::from_secs(10));
        targets = benchmark_product_isogeny,
    }
}

fn main() {
    benchmark_product::benchmark_product();
}
