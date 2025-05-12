// SQISign level 1
static SQISIGN_ONE_MODULUS: [u64; 4] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0x04FFFFFFFFFFFFFF,
];

fp2::define_fp2_from_modulus!(
    typename = SqiSqignI,
    base_typename = SqiSignIBase,
    modulus = SQISIGN_ONE_MODULUS,
);
