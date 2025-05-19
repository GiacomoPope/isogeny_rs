// TODO: something in the macro has too many arguements... Catch this in fp2 and fix it.
#![allow(clippy::too_many_arguments)]

// SQISign level 1
static SQISIGN_I_MODULUS: [u64; 4] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0x04FFFFFFFFFFFFFF,
];

fp2::define_fp2_from_modulus!(
    typename = SqiField248,
    base_typename = SqiField248Base,
    modulus = SQISIGN_I_MODULUS,
);

#[cfg(test)]
mod test_sqisign_arithmetic {
    use super::{SQISIGN_I_MODULUS, SqiField248, SqiField248Base};

    fp2::define_fp_tests!(SqiField248Base);
    fp2::define_fp2_tests!(SqiField248, SQISIGN_I_MODULUS, 5);
}
