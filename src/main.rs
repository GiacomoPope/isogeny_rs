#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use fp2;
use isogeny::curve::Curve;

fn main() {
    // Characteristic: p = 5*2^248 - 1
    pub static ea: usize = 248;

    // Supersingular elliptic curve with Montgomery coeff A0
    pub static A0_str: &str = "5afe2b2360b0f02d30eb4c7e9180f5fe1b9049a0a56c82bf3819b10fb1e4f801940705ebcd0d4bc103dc2a456cd7b12ccff444e30cafff6c54963cc19159e804";

    // <P, Q> are the generators of the curve of order N = 2^ea
    pub static PX_str: &str = "618d5b89152bdbf31078931105c5c6b3c449519ad498f273b058c54d611f2904df88997980684f4a4f8db581c6cce1aca5e17a2d1e114e8346e0a5881b90bf01";
    pub static QX_str: &str = "094276e052db1ce67a71192db1fb0e1c4b41d9bb224fc4e77e7c18e04a584b001f7b31453b2ca17ad36b309b91b7a7c380b442a7e28448446b22ec3513889703";
    // P - Q
    pub static PmQX_str: &str = "7b621910b3f2610e4de0b4fd07c9c444f417ea5b9b7ea91f93b9faa6bc8d3503bfdbc3174c32855427e7c675e16ef89b9ba80b911d582823ee83359becdd3800";

    // dA = (p + 1) // A
    pub static dA: [u8; 1] = [5];
    pub static dA_BITLEN: usize = 3;

    // Tate pairings for order: e_2^ea(P, Q)
    pub static tate_pairing_str: &str = "7558fbdd7877502f500af6cb1a0aaf240545ca8ee0514fb0111a6b097d2fa60131292a781c0418806da296196dfa1e294bbc5a82f0a7a4efe26e20837fe81f03";

    // Weil pairings for order: e_2^ea(P, Q)
    pub static weil_pairing_str: &str = "83b9cd40d9f7cf693b171ff34ce87cbcbb32fb194615790b2a35c252979f170470e6468c7a1f06963a33f70375fb3a4900c9df8ba90abda41322d82baadf2003";

    // Fp251: a finite field element GF(p) with p = 3 mod 4.
    fp2::define_fp_core!(
        typename = Fp251,
        modulus = [
            0xFFFFFFFFFFFFFFFFu64,
            0xFFFFFFFFFFFFFFFFu64,
            0xFFFFFFFFFFFFFFFFu64
        ],
    );

    // Fp251Ext: a finite field element GF(p^2) with modulus x^2 + 1.
    fp2::define_fp2_core!(typename = Fp251Ext, base_field = Fp251,);

    let (A0, _) = <Fp251Ext>::decode(&hex::decode(A0_str).unwrap());
    let E0 = Curve::new(&A0);

    // Point of order 2^ea x(P)
    let (xP, _) = <Fp251Ext>::decode(&hex::decode(PX_str).unwrap());

    // Point of order 2^ea x(Q)
    let (xQ, _) = <Fp251Ext>::decode(&hex::decode(QX_str).unwrap());

    // Difference point x(P - Q)
    let (xPQ, _) = <Fp251Ext>::decode(&hex::decode(PmQX_str).unwrap());

    // Weil and Tate Pairing computed via SageMath
    let (eQ2P2_weil, _) = <Fp251Ext>::decode(&hex::decode(weil_pairing_str).unwrap());
    let (eQ2P2_tate, _) = <Fp251Ext>::decode(&hex::decode(tate_pairing_str).unwrap());

    // Compute Weil pairing
    let eQ2P2_weil_test = E0.weil_pairing_2exp(&xP, &xQ, &xPQ, ea);
    println!("{}", eQ2P2_weil.equals(&eQ2P2_weil_test) == u32::MAX);

    // Compute Tate pairing
    let eQ2P2_tate_test = E0.tate_pairing_2exp(&xP, &xQ, &xPQ, ea, &dA, dA_BITLEN);

    println!("{}", eQ2P2_tate.equals(&eQ2P2_tate_test) == u32::MAX);
}
