
// CSIDH - 512
const CSIDH_512_MODULUS: [u64; 8] = [
    0x1B81B90533C6C87B,
    0xC2721BF457ACA835,
    0x516730CC1F0B4F25,
    0xA7AAC6C567F35507,
    0x5AFBFCC69322C9CD,
    0xB42D083AEDC88C42,
    0xFC8AB0D15E3E4C4A,
    0x65B48E8F740F89BF,
];



fp2::define_fp_core!(
    typename = Csidh512,
    modulus = CSIDH_512_MODULUS,
);

#[cfg(test)]
mod test_csidh_arithmetic {
    use super::{Csidh512};

    fp2::define_fp_tests!(Csidh512);
}
