use super::sqisign::Sqisign;
use crate::fields::sqisign::{SqiField248, SqiField376, SqiField500};

mod level_one {
    use crate::protocols::sqisign::SqisignParameters;

    const SECURITY_BITS: usize = 128;
    const COFACTOR: u8 = 5;
    const F: usize = 248;
    const RESPONSE_LENGTH: usize = 126;
    const HASH_ITERATIONS: usize = 64;
    const PK_LENGTH: usize = 65;
    const SK_LENGTH: usize = 353;
    const SIG_LENGTH: usize = 148;

    pub const SQISIGN_PARAMS_I: SqisignParameters = SqisignParameters {
        security_bits: SECURITY_BITS,
        cofactor: COFACTOR,
        f: F,
        response_length: RESPONSE_LENGTH,
        hash_iterations: HASH_ITERATIONS,
        pk_len: PK_LENGTH,
        sk_len: SK_LENGTH,
        sig_len: SIG_LENGTH,
    };
}

mod level_three {
    use crate::protocols::sqisign::SqisignParameters;

    const SECURITY_BITS: usize = 192;
    const COFACTOR: u8 = 65;
    const F: usize = 376;
    const RESPONSE_LENGTH: usize = 192;
    const HASH_ITERATIONS: usize = 256;
    const PK_LENGTH: usize = 97;
    const SK_LENGTH: usize = 529;
    const SIG_LENGTH: usize = 224;

    pub const SQISIGN_PARAMS_III: SqisignParameters = SqisignParameters {
        security_bits: SECURITY_BITS,
        cofactor: COFACTOR,
        f: F,
        response_length: RESPONSE_LENGTH,
        hash_iterations: HASH_ITERATIONS,
        pk_len: PK_LENGTH,
        sk_len: SK_LENGTH,
        sig_len: SIG_LENGTH,
    };
}

mod level_five {
    use crate::protocols::sqisign::SqisignParameters;

    const SECURITY_BITS: usize = 256;
    const COFACTOR: u8 = 27;
    const F: usize = 500;
    const RESPONSE_LENGTH: usize = 253;
    const HASH_ITERATIONS: usize = 512;
    const PK_LENGTH: usize = 129;
    const SK_LENGTH: usize = 701;
    const SIG_LENGTH: usize = 292;

    pub const SQISIGN_PARAMS_V: SqisignParameters = SqisignParameters {
        security_bits: SECURITY_BITS,
        cofactor: COFACTOR,
        f: F,
        response_length: RESPONSE_LENGTH,
        hash_iterations: HASH_ITERATIONS,
        pk_len: PK_LENGTH,
        sk_len: SK_LENGTH,
        sig_len: SIG_LENGTH,
    };
}

pub const SQISIGN_I: Sqisign<SqiField248> = Sqisign::new(&level_one::SQISIGN_PARAMS_I);
pub const SQISIGN_III: Sqisign<SqiField376> = Sqisign::new(&level_three::SQISIGN_PARAMS_III);
pub const SQISIGN_V: Sqisign<SqiField500> = Sqisign::new(&level_five::SQISIGN_PARAMS_V);
