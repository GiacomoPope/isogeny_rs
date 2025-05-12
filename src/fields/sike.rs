// NIST lvl 1 SIKE prime: p = 2^216 * 3^137 - 1
// Fp434Ext: a finite field element GF(p^2) with modulus x^2 + 1.
static SIKE_I_MODULUS: [u64; 7] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFDC1767AE2FFFFFF,
    0x7BC65C783158AEA3,
    0x6CFC5FD681C52056,
    0x0002341F27177344,
];
fp2::define_fp2_from_modulus!(
    typename = SikeOne,
    base_typename = SikeOneBase,
    modulus = SIKE_I_MODULUS,
);
