
// CSIDH - 512
const CSIDH_512_MODULUS: [u64; 8] = [
    0x1b81b90533c6c87b, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
    0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf,
];



fp2::define_fp2_from_modulus!(
    typename = Csidh512,
    base_typename = Csidh512Base,
    modulus = CSIDH_512_MODULUS,
);

#[cfg(test)]
mod test_csidh_arithmetic {
    use super::{Csidh512Base};

    fp2::define_fp_tests!(Csidh512Base);
    //fp2::define_fp2_tests!(Csidh512, CSIDH_512_MODULUS, 2);
}
