// SQISign level 1: p = 5 * 2^248 - 1
const SQISIGN_I_MODULUS: [u64; 4] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0x04FFFFFFFFFFFFFF,
];

// SQISign level 3: p = 65 * 2^376 - 1
const SQISIGN_III_MODULUS: [u64; 6] = [
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0x40FFFFFFFFFFFFFF,
];

// SQISign level 5: p = 27 * 2^500 - 1
const SQISIGN_V_MODULUS: [u64; 8] = [
    0xFFFFFFFFFFFFFFFFu64,
    0xFFFFFFFFFFFFFFFFu64,
    0xFFFFFFFFFFFFFFFFu64,
    0xFFFFFFFFFFFFFFFFu64,
    0xFFFFFFFFFFFFFFFFu64,
    0xFFFFFFFFFFFFFFFFu64,
    0xFFFFFFFFFFFFFFFFu64,
    0x01AFFFFFFFFFFFFFu64,
];

fp2::define_fp2_from_modulus!(
    typename = SqiField248,
    base_typename = SqiField248Base,
    modulus = SQISIGN_I_MODULUS,
);

fp2::define_fp2_from_modulus!(
    typename = SqiField376,
    base_typename = SqiField376Base,
    modulus = SQISIGN_III_MODULUS,
);

fp2::define_fp2_from_modulus!(
    typename = SqiField500,
    base_typename = SqiField500Base,
    modulus = SQISIGN_V_MODULUS,
);

use crate::fields::gf_248::Fp248;
fp2::define_fp2_from_type!(typename = Fp248Ext, base_field = Fp248,);

#[cfg(test)]
mod test_sqisign_i_arithmetic {
    use super::{SQISIGN_I_MODULUS, SqiField248, SqiField248Base};

    fp2::define_fp_tests!(SqiField248Base);
    fp2::define_fp2_tests!(SqiField248, SQISIGN_I_MODULUS, 5);
}

#[cfg(test)]
mod test_sqisign_iii_arithmetic {
    use super::{SQISIGN_III_MODULUS, SqiField376, SqiField376Base};

    fp2::define_fp_tests!(SqiField376Base);
    fp2::define_fp2_tests!(SqiField376, SQISIGN_III_MODULUS, 6);
}

#[cfg(test)]
mod test_sqisign_v_arithmetic {
    use super::{SQISIGN_V_MODULUS, SqiField500, SqiField500Base};

    fp2::define_fp_tests!(SqiField500Base);
    fp2::define_fp2_tests!(SqiField500, SQISIGN_V_MODULUS, 4);
}

#[cfg(test)]
mod test_sqisign_i_alt_arithmetic {
    use super::super::gf_248::Fp248;
    use super::Fp248Ext;
    fp2::define_fp_tests!(Fp248);
    fp2::define_fp2_tests!(Fp248Ext, Fp248::MODULUS, 5);
}
