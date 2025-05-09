#![allow(non_snake_case)]

// TODO: further randomize this test
#[cfg(test)]
mod test_point_compression {

    use isogeny::elliptic::curve::Curve;

    static MODULUS: [u64; 4] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x04FFFFFFFFFFFFFF,
    ];
    fp2::define_fp2_from_modulus!(
        typename = Fp248Ext,
        base_typename = Fp248,
        modulus = MODULUS,
    );

    // Characteristic: p = 5*2^248 - 1
    static F: usize = 248;
    static E: usize = 124;

    // Supersingular elliptic curve with Montgomery coeff A0
    static A0_STR: &str = "5afe2b2360b0f02d30eb4c7e9180f5fe1b9049a0a56c82bf3819b10fb1e4f801940705ebcd0d4bc103dc2a456cd7b12ccff444e30cafff6c54963cc19159e804";

    // <P, Q> are the generators of the curve of order N = 2^f
    static PX_STR: &str = "54e7eead54a6966230bd4f2cbb1f21e98bfc0c92e29873d7692a32a3082003030f7a39fda4a2c5fdcfd603811e849899ac2a8fa1546ec59e5dcab3ce81efbd02";
    static QX_STR: &str = "454eead047b869fe64009c88b4b07d41ac8e36669eb561336258f8f95925df039622947853adb064403b42f487de1b9f8284f144b135cd2849d4c7a5e73a9b01";
    static P_QX_STR: &str = "9e869f279df7d63583bc37328cf538e1b8818f34b9693715cf88620b8ebe8102a32d186ba5227c8276da60c2f952ced4623f6526b299396eab0226b5645e0201";

    // <R, S> are the generators of the curve of order N = 2^e
    static RX_STR: &str = "c361fc21e6a92d7f7d10d9f373a5ca7efb158881408645cd7ca605c0f1c17404ff4a7be793aecea93cb7b0b76df3c95dce8e9af09bf0ab29b3e38d2b1d378402";
    static SX_STR: &str = "3a4d0e7834f5ebbaceb49cd83a0747b226eea16252a0066c4db9383eb11b76022a6b25cb60d7507334e51a305b7bac12b2feec3dd65231905a08138d6f73fc02";
    static R_SX_STR: &str = "881733702997e9667a5263611949fd88779d917b09a5bc20756f7d3931b30d0368488845c39c79bac41c2fca46e92cc39b7ea9a6a8a3f9aebca406032c55eb00";

    // dA = (p + 1) // 2^f
    static DA: [u8; 1] = [5];
    static DA_BITLEN: usize = 3;

    // Randomized variables to find with dlogs
    static A: [u8; 16] = [
        231, 236, 12, 226, 96, 63, 8, 148, 235, 149, 78, 41, 234, 78, 223, 10,
    ];
    static B: [u8; 16] = [
        29, 36, 137, 134, 164, 32, 173, 196, 36, 87, 201, 117, 70, 176, 216, 10,
    ];
    static C: [u8; 16] = [
        13, 1, 11, 179, 59, 21, 131, 227, 185, 72, 84, 22, 159, 90, 10, 14,
    ];
    static D: [u8; 16] = [
        237, 22, 11, 123, 195, 142, 41, 178, 254, 241, 222, 40, 213, 125, 223, 5,
    ];

    #[test]
    fn test_point_compression() {
        // Create E0
        // A supersingular elliptic curve of order 5*2^248
        let (A0, check) = <Fp248Ext>::decode(&hex::decode(A0_STR).unwrap());
        assert!(check == u32::MAX);
        let E0 = Curve::new(&A0);

        // Generator x(P), x(Q), x(P - Q)
        let (xP, check) = <Fp248Ext>::decode(&hex::decode(PX_STR).unwrap());
        assert!(check == u32::MAX);
        let (xQ, check) = <Fp248Ext>::decode(&hex::decode(QX_STR).unwrap());
        assert!(check == u32::MAX);
        let (xPQ, check) = <Fp248Ext>::decode(&hex::decode(P_QX_STR).unwrap());
        assert!(check == u32::MAX);

        // Generator x(R), x(S), x(R - S)
        let (xR, check) = <Fp248Ext>::decode(&hex::decode(RX_STR).unwrap());
        assert!(check == u32::MAX);
        let (xS, check) = <Fp248Ext>::decode(&hex::decode(SX_STR).unwrap());
        assert!(check == u32::MAX);
        let (xRS, check) = <Fp248Ext>::decode(&hex::decode(R_SX_STR).unwrap());
        assert!(check == u32::MAX);

        let (r1, r2, s1, s2, check) =
            E0.point_compression(&xP, &xQ, &xPQ, &xR, &xS, &xRS, F, E, &DA, DA_BITLEN);

        assert!(r1 == A);
        assert!(r2 == B);
        assert!(s1 == C);
        assert!(s2 == D);
        assert!(check == u32::MAX);
    }
}
