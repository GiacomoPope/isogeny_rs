use super::csidh::Csidh;
use crate::fields::csidh::Csidh512;

mod csidh_512 {
    use crate::protocols::csidh::CsidhParameters;

    const NUM_PRIMES: usize = 74;
    const MAX_EXPONENT: usize = 5;
    const COFACTOR: u64 = 4;
    const PRIMES: [u64; NUM_PRIMES] = [
         3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
        61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
       139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
       229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
       317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
    ];

    pub const CSIDH_PARAMS: CsidhParameters = CsidhParameters {
        num_primes: NUM_PRIMES,
        max_exponent: MAX_EXPONENT,
        two_cofactor: COFACTOR,
        primes: PRIMES,
    };
}


pub const CSIDH_512: Csidh<Csidh512> = Csidh::new(&csidh_512::CSIDH_PARAMS);
