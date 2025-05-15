use super::sqisign::Sqisign;
use crate::fields::sqisign::SqiSignI as Fp2;

static SECURITY_BITS: usize = 128;
static COFACTOR: u8 = 5;
static COFACTOR_BITSIZE: usize = 3; // TODO: could do a mul small for u8 type?
static F: usize = 248;
static RESPONSE_LENGTH: usize = 126;
static HASH_ITERATIONS: usize = 64;
static PK_LENGTH: usize = 65;
static SK_LENGTH: usize = 353;
static SIG_LENGTH: usize = 148;

pub static SQISIGN_I: Sqisign<Fp2> = Sqisign::new(
    SECURITY_BITS,
    COFACTOR,
    COFACTOR_BITSIZE,
    F,
    RESPONSE_LENGTH,
    HASH_ITERATIONS,
    PK_LENGTH,
    SK_LENGTH,
    SIG_LENGTH,
);
