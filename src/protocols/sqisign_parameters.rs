use super::sqisign::{Sqisign, SqisignParameters};
use crate::fields::sqisign::SqiField248 as Fp2;

const SECURITY_BITS: usize = 128;
const COFACTOR: u8 = 5;
const COFACTOR_BITSIZE: usize = 3; // TODO: could do a mul small for u8 type?
const F: usize = 248;
const RESPONSE_LENGTH: usize = 126;
const HASH_ITERATIONS: usize = 64;
const PK_LENGTH: usize = 65;
const SK_LENGTH: usize = 353;
const SIG_LENGTH: usize = 148;

const SQISIGN_PARAMS_I: SqisignParameters = SqisignParameters {
    security_bits: SECURITY_BITS,
    cofactor: COFACTOR,
    cofactor_bitsize: COFACTOR_BITSIZE,
    f: F,
    response_length: RESPONSE_LENGTH,
    hash_iterations: HASH_ITERATIONS,
    pk_len: PK_LENGTH,
    sk_len: SK_LENGTH,
    sig_len: SIG_LENGTH,
};

pub const SQISIGN_I: Sqisign<Fp2> = Sqisign::new(&SQISIGN_PARAMS_I);
