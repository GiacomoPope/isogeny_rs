use super::sqisign::Sqisign;
use crate::fields::sqisign::SqiSignI as Fp2;

static SECURITY_BITS: usize = 128;
static COFACTOR: u8 = 5;
static F: usize = 248;
static RESPONSE_LENGTH: usize = 126;
static HASH_ITERATIONS: usize = 64;
static PK_LENGTH: usize = 65;
static SK_LENGTH: usize = 353;
static SIG_LENGTH: usize = 148;

pub static SQISIGN_I: Sqisign<Fp2> = Sqisign::new(
    SECURITY_BITS,
    COFACTOR,
    F,
    RESPONSE_LENGTH,
    HASH_ITERATIONS,
    PK_LENGTH,
    SK_LENGTH,
    SIG_LENGTH,
);
