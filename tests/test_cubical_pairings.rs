#![allow(non_snake_case)]

#[cfg(test)]
mod test_pairings {
    use fp2::fq::Fq;
    use isogeny::elliptic::curve::Curve;
    use isogeny::elliptic::point::PointX;
    use isogeny::fields::sike::SikeOne as Fp2;

    // Characteristic: p = 2^216 * 3^137 - 1
    static EA: usize = 216;

    // Supersingular elliptic curve with Montgomery coeff A0
    static A0_STR: &str = "9d9ab33208c854c7c5d70b7429e8ed3c8101cf22464be069336d00de3adfd028826eaa63ef47728a9272d2fe76215171f0fa84bb5e22007d6199a0c10dfae6185579efe044ce22560c3bcedbf594753a52d7568d0857ea804273d176e6e60f1d5a14749e03051105ddd8ccb99701";

    // Below are x(P), x(Q) and x(P - Q) for various bases <P, Q>
    //
    // <P, Q> are the generators of the curve of order N = 2^216 * 3^137
    static PX_STR: &str = "6a372af6de8ee92fe4f1d52990ec0c84648e05254cf87a961fd1d390f4638d48660c2e4ca17b8d5a93fcb59a9efc0028b39db5e9ad2601c2ce85c2eff902a44e2dd4c35a3f8d004967ff20e4c2f976e91373cc648a27fd6bd6ea4ae235b93af5d3968e125948255415d7362e3b01";
    static QX_STR: &str = "3a0faf61f347560686560e692187bfcaff6e7889e442b3ded053c9349fc2e4af4d31cfa7252e023985af73d1064fb80006b3087e490a027feba57e3dc27884b92b0ed7cf9919ab1e02a506d6983edf06433959eb56af2752393aab619f9db4d18ff4f89546befa619d8a297a1300";
    static PMQX_STR: &str = "f2a676ca9e354357c6647e4928ccbf53c9189361bf4c77daea2d9bb42ccaeb475ad4d7f9542d54de90d3df099d59d26d05d471ca754b001ca3948b8f1b226e9a09e9eca4181d079ecf003488105d29db60d509567c40de455fea0297f70c6b5d3064e8158c1dcee308797c2ae301";

    // <P2, Q2> are the generators of E0[A]
    static P2X_STR: &str = "3412971c46b2e48d9fc16de5a7111ca678228033058378a68525cacca39f9041e8a9dba3b7457ce0072e4ad30761bc9baf2e6b474b6600f30e7495bbe4ae5c336d22615250311380208dee97b12f22904ded0c322bab96d65b1543713ba1a1b7cd78f0d44293c98b12ad5b4c4b00";
    static Q2X_STR: &str = "5a79e05fd69c957a6bacc4459574de4400e7c366644a2c11b4dbda6b479728aff2280061864a3beed2856a260115ed8f08e99f55fe5b0039e8522f8b892484a9bb86e063c6316261e962be80dcf13988a749f007f399a84c8a98031f758da061f764c14a38ac99e98638815f0200";
    static P2MQ2X_STR: &str = "640b80005988fd480611b498ddee1e34b5d8b9b27ec6c04ddab0cd3c8c183f596f2ef98bd6dcba97e7dff02dc9c4d63816fd05a74820019e9274d996c135e8bb12f46177fc48fd0246dee712e696c0cf52aacb7f6dcd04b8ca3c5080a02e39c8813a050956d4c659c4fb476cb001";

    // <P3, Q3> are the generators of E0[B]
    static P3X_STR: &str = "9811b1f089f1c6a24361fbc7e0ae603581bc635a97ca89ec314df1318beaea2ea9d651d79f7bcc6136b829ba32bd7f32b2c3443a5412000c0427996562ed0b1792d366605e5e8ff49d06587131d0455a5e37e6dc26c46146e80388ced69fd69a258eed48e90df9a438efec9a0502";
    static Q3X_STR: &str = "ff868d58f50206db4405ddec9a254c7ab11d8fec33342a42c7ac792a283aff045313379361d407e69fb668382604123b863e1edbecce011001f855b9f581eab459ad4426588dfb4952a2925d68483d26f173b925037a9b200dd2a1acfe2a47e6ed7a9a25b9bbce78e9f429057300";
    static P3MQ3X_STR: &str = "044a12aa7153ed1af74775442eb8a54ae7405f9f465d810494895b421ff69df622ce8d134a65441d6acf818693367857358f5b88dfa201e972a0c223b307552ee19574d5e35041bba32a1631a8a10a9bd19bbd8625030b40eeedc7a7691fcfc4398065a5f53eef947d3b6a58c901";

    // [A]P + Q
    static AP_QX_STR: &str = "139165d5f73c15240bc3b8a70d25b1de1337af6ac0e4b9cfafbcb67cf2df6c08ca2641ebd3a251674fffa0ded8c931d882cd0c18313a0161dd34892ca9ab14cf34ba504ee7fae92dcb68477482ef7c271921176e517c3f7cc853f043ca1e1b32fecc6bcfcd127b9cc0c000dbde00";

    // [B]P + Q
    static BP_QX_STR: &str = "ef83e5079133d3c3fc30b6a0bbbb29cdf4dda1d23a61c564ce5c45c71a369a6ab947e15ecf419f22c415729c7073464b77bf7585b95100fc60a327bdd4a7a821f9c6a8fcf7a2ed34a8bbf385d916854779130c34a798ed3e075ecf21db57f411d11f6b929581206f3150f3ab4000";

    // B = 3^137
    static B: [u8; 28] = [
        227, 122, 118, 193, 253, 163, 174, 88, 49, 120, 92, 198, 123, 86, 32, 197, 129, 214, 95,
        252, 108, 68, 115, 23, 39, 31, 52, 2,
    ];
    static B_BITLEN: usize = 218;

    // N = 2^216 * 3^137
    static N: [u8; 55] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 227, 122,
        118, 193, 253, 163, 174, 88, 49, 120, 92, 198, 123, 86, 32, 197, 129, 214, 95, 252, 108,
        68, 115, 23, 39, 31, 52, 2,
    ];
    static N_BITLEN: usize = 434;

    // dX = (p + 1) // X
    static DA: [u8; 28] = [
        227, 122, 118, 193, 253, 163, 174, 88, 49, 120, 92, 198, 123, 86, 32, 197, 129, 214, 95,
        252, 108, 68, 115, 23, 39, 31, 52, 2,
    ];
    static DA_BITLEN: usize = 218;
    static DB: [u8; 28] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    ];
    static DB_BITLEN: usize = 217;
    static DN: [u8; 1] = [1];
    static DN_BITLEN: usize = 1;

    // Tate pairings for order N, A, B: e_X(Q, P)
    static TATE_N_STR: &str = "03a5fc5baf13feee0158f691bbdf53a904f23ae2097d06f376f8ff89dbec74beb144c9884b6283dcdd0f37d0b5f07b1a8bfa1089b3b60197282228e36be9b1f8b94e38f11ce7a2c49d1e15776bc1d81eed893486a9fb540fd1e06bbd59720856d1ee535090ec25ebfca7d4a36f01";
    static TATE_A_STR: &str = "337355edb6fdcc2fcde4bb95cdc383b9cacf126b40887372ed5e84f3303e91b610946413d14c42007824583ca4282c62678010fa971e00db7220bdf3fec4dd8493d5bc6d9121ae92d139a44b20147c2c9a75017609c0165064e9494b5ec36c96c1e8207f6c50ad583b376bb2e200";
    static TATE_B_STR: &str = "665ed36f2b5f5c83ac7b9b0bc3d33d82812a6da55486efca7f76b74b6ddc4718e8d82cf11ad6fa863f12727f8879736b52cc4ba6f932015925f6e6260b3789c6f5146228156cc4ee96b8d88251cb135bbee94f79465932b0692604ccf8690edafb1c65285b9a2585e03f6454d300";

    // Weil pairings for order N, A, B: e_X(Q, P)
    static WEIL_N_STR: &str = "03a5fc5baf13feee0158f691bbdf53a904f23ae2097d06f376f8ff89dbec74beb144c9884b6283dcdd0f37d0b5f07b1a8bfa1089b3b60168d7ddd71c94164e0746b1c70ee3185d3b62e1ea88943e27e11276aef4ccc5a894dd77c5ba025473004fd62d86cf0f4759766f527bc400";
    static WEIL_A_STR: &str = "f2f160e0590892f8e64107b1ca763b0e995a778762683ade3ff618f38161e76ff85dfcf1ca0271239ed8baec56c49bcdbb85bdd84bc301484c3355325156eb01e282fa47f24e0aedc0ea6c9a87a761499112c37cfeb3a2a9e96550760f1f039918d3e0974cdbe8e5433361fc3100";
    static WEIL_B_STR: &str = "13125102956fd55f99301ee8743c8139a6d361d893236464bc56f8292aa2653bbeb01fac2daaaa6a7344f5836c3d873bd76796d7a3750038230b4047ead67767ae158126ace065c2d3239c5eb3b3454ca62b318d55b8f975faa6ffa9426187d5efed64529b95ad0abdce17e3c800";

    #[test]
    fn test_cubical_ladder() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, _) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        let E0 = Curve::new(&A0);

        // Generator x(P)
        let (xP, r) = Fp2::decode(&hex::decode(PX_STR).unwrap());
        assert!(r == u32::MAX);

        // Generator x(Q)
        let (xQ, r) = Fp2::decode(&hex::decode(QX_STR).unwrap());
        assert!(r == u32::MAX);

        // Difference point x(P - Q)
        let (xPQ, r) = Fp2::decode(&hex::decode(PMQX_STR).unwrap());
        assert!(r == u32::MAX);

        // P2 = [B] P (Even torsion part)
        let (xP2, r) = Fp2::decode(&hex::decode(P2X_STR).unwrap());
        assert!(r == u32::MAX);
        let P2 = PointX::new(&xP2, &Fp2::ONE);

        // [B]P + Q
        let (BP_Q_X, r) = Fp2::decode(&hex::decode(BP_QX_STR).unwrap());
        assert!(r == u32::MAX);
        let BP_Q = PointX::new(&BP_Q_X, &Fp2::ONE);

        // Compute [n]P and [n]P + Q
        let (nP_test, nPQ_test) = E0.cubical_ladder(&xP, &xQ, &xPQ, &B, B_BITLEN, false);

        assert!(P2.equals(&nP_test) == u32::MAX);
        assert!(BP_Q.equals(&nPQ_test) == u32::MAX);
    }

    #[test]
    fn test_cubical_ladder_2exp() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // Generator x(P)
        let (xP, r) = Fp2::decode(&hex::decode(PX_STR).unwrap());
        assert!(r == u32::MAX);

        // Generator x(Q)
        let (xQ, r) = Fp2::decode(&hex::decode(QX_STR).unwrap());
        assert!(r == u32::MAX);

        // Difference point x(P - Q)
        let (xPQ, r) = Fp2::decode(&hex::decode(PMQX_STR).unwrap());
        assert!(r == u32::MAX);

        // P3 = [2^ea] P (Odd torsion part)
        let (xP3, r) = Fp2::decode(&hex::decode(P3X_STR).unwrap());
        assert!(r == u32::MAX);
        let P3 = PointX::new(&xP3, &Fp2::ONE);

        // [2^ea] P + Q
        let (AP_Q_X, r) = Fp2::decode(&hex::decode(AP_QX_STR).unwrap());
        assert!(r == u32::MAX);
        let AP_Q = PointX::new(&AP_Q_X, &Fp2::ONE);

        let (nP_test, nPQ_test) = E0.cubical_ladder_2exp(&xP, &xQ, &xPQ, EA);

        assert!(P3.equals(&nP_test) == u32::MAX);
        assert!(AP_Q.equals(&nPQ_test) == u32::MAX);
    }

    #[test]
    fn test_tate_pairing_2exp() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // x-coordinate x(P) of order 2^ea
        let (xP, r) = Fp2::decode(&hex::decode(P2X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(Q) of order 2^ea
        let (xQ, r) = Fp2::decode(&hex::decode(Q2X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(P - Q) of order 2^ea
        let (xPQ, r) = Fp2::decode(&hex::decode(P2MQ2X_STR).unwrap());
        assert!(r == u32::MAX);

        // Tate Pairing computed via SageMath
        let (eQ2P2, r) = Fp2::decode(&hex::decode(TATE_A_STR).unwrap());
        assert!(r == u32::MAX);

        // Compute Tate pairing
        let eQ2P2_test = E0.tate_pairing_2exp(&xP, &xQ, &xPQ, EA, &DA, DA_BITLEN);

        assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);
    }

    #[test]
    fn test_weil_pairing_2exp() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // x-coordinate x(P) of order 2^ea
        let (xP, r) = Fp2::decode(&hex::decode(P2X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(Q) of order 2^ea
        let (xQ, r) = Fp2::decode(&hex::decode(Q2X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(P - Q) of order 2^ea
        let (xPQ, r) = Fp2::decode(&hex::decode(P2MQ2X_STR).unwrap());
        assert!(r == u32::MAX);

        // Weil Pairing computed via SageMath
        let (eQ2P2, r) = Fp2::decode(&hex::decode(WEIL_A_STR).unwrap());
        assert!(r == u32::MAX);

        // Compute Weil pairing
        let eQ2P2_test = E0.weil_pairing_2exp(&xP, &xQ, &xPQ, EA);

        assert!(eQ2P2.equals(&eQ2P2_test) == u32::MAX);
    }

    #[test]
    fn test_tate_pairing_odd() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // x-coordinate x(P) of order 3^eb
        let (xP, r) = Fp2::decode(&hex::decode(P3X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(Q) of order 3^eb
        let (xQ, r) = Fp2::decode(&hex::decode(Q3X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(P-Q) of order 3^eb
        let (xPQ, r) = Fp2::decode(&hex::decode(P3MQ3X_STR).unwrap());
        assert!(r == u32::MAX);

        // Tate Pairing computed via SageMath
        let (eQ3P3, r) = Fp2::decode(&hex::decode(TATE_B_STR).unwrap());
        assert!(r == u32::MAX);

        // Compute Tate pairing
        let eQ3P3_test = E0.tate_pairing(&xP, &xQ, &xPQ, &B, B_BITLEN, &DB, DB_BITLEN);

        assert!(eQ3P3.equals(&eQ3P3_test) == u32::MAX);
    }

    #[test]
    fn test_weil_pairing_odd() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // x-coordinate x(P) of order 3^eb
        let (xP, r) = Fp2::decode(&hex::decode(P3X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(Q) of order 3^eb
        let (xQ, r) = Fp2::decode(&hex::decode(Q3X_STR).unwrap());
        assert!(r == u32::MAX);

        // x-coordinate x(P-Q) of order 3^eb
        let (xPQ, r) = Fp2::decode(&hex::decode(P3MQ3X_STR).unwrap());
        assert!(r == u32::MAX);

        // Weil Pairing computed via SageMath
        let (eQ3P3, r) = Fp2::decode(&hex::decode(WEIL_B_STR).unwrap());
        assert!(r == u32::MAX);

        // Compute Weil pairing
        let eQ3P3_test = E0.weil_pairing(&xP, &xQ, &xPQ, &B, B_BITLEN);

        assert!(eQ3P3.equals(&eQ3P3_test) == u32::MAX);
    }

    #[test]
    fn test_tate_pairing_even() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // Generator x(P)
        let (xP, r) = Fp2::decode(&hex::decode(PX_STR).unwrap());
        assert!(r == u32::MAX);

        // Generator x(Q)
        let (xQ, r) = Fp2::decode(&hex::decode(QX_STR).unwrap());
        assert!(r == u32::MAX);

        // Difference point x(P - Q)
        let (xPQ, r) = Fp2::decode(&hex::decode(PMQX_STR).unwrap());
        assert!(r == u32::MAX);

        // Tate Pairing computed via SageMath
        let (eQP, r) = Fp2::decode(&hex::decode(TATE_N_STR).unwrap());
        assert!(r == u32::MAX);

        // Compute Tate pairing
        let eQP_test = E0.tate_pairing(&xP, &xQ, &xPQ, &N, N_BITLEN, &DN, DN_BITLEN);

        assert!(eQP.equals(&eQP_test) == u32::MAX);
    }

    #[test]
    fn test_weil_pairing_even() {
        // Create E0
        // A supersingular elliptic curve of order 2^ea * 3^eb
        let (A0, r) = Fp2::decode(&hex::decode(A0_STR).unwrap());
        assert!(r == u32::MAX);
        let E0 = Curve::new(&A0);

        // Generator x(P)
        let (xP, r) = Fp2::decode(&hex::decode(PX_STR).unwrap());
        assert!(r == u32::MAX);

        // Generator x(Q)
        let (xQ, r) = Fp2::decode(&hex::decode(QX_STR).unwrap());
        assert!(r == u32::MAX);

        // Difference point x(P - Q)
        let (xPQ, r) = Fp2::decode(&hex::decode(PMQX_STR).unwrap());
        assert!(r == u32::MAX);

        // Weil Pairing computed via SageMath
        let (eQP, r) = Fp2::decode(&hex::decode(WEIL_N_STR).unwrap());
        assert!(r == u32::MAX);

        // Compute Tate pairing
        let eQP_test = E0.weil_pairing(&xP, &xQ, &xPQ, &N, N_BITLEN);

        assert!(eQP.equals(&eQP_test) == u32::MAX);
    }
}
