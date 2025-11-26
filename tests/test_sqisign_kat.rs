#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

#[cfg(feature = "kat-tests")]
mod test_helpers {
    use std::path::Path;

    use fp2::traits::Fp2 as FqTrait;

    use isogeny::protocols::sqisign::Sqisign;

    #[derive(Debug, Clone, Default)]
    pub struct KatTest {
        pub count: usize,
        pub seed: Vec<u8>,
        pub mlen: usize,
        pub msg: Vec<u8>,
        pub pk_bytes: Vec<u8>,
        pub sk_bytes: Vec<u8>,
        pub smlen: usize,
        pub sm: Vec<u8>,
    }

    fn parse_kat_test(block: &str, test: &mut KatTest) {
        for line in block.lines() {
            let (name, value) = line.split_once(" = ").unwrap();
            match name {
                "count" => test.count = value.parse().unwrap(),
                "seed" => test.seed = hex::decode(value).unwrap(),
                "mlen" => test.mlen = value.parse().unwrap(),
                "msg" => test.msg = hex::decode(value).unwrap(),
                "pk" => test.pk_bytes = hex::decode(value).unwrap(),
                "sk" => test.sk_bytes = hex::decode(value).unwrap(),
                "smlen" => test.smlen = value.parse().unwrap(),
                "sm" => test.sm = hex::decode(value).unwrap(),
                _ => panic!("Unexpected key: {name:?}"),
            }
        }
    }

    pub fn parse_kat_tests(name: &str) -> Vec<KatTest> {
        // Load the KAT file as a string
        let response_file = Path::new("assets")
            .join("sqisign")
            .join(format!("{}.rsp", name));
        let contents = std::fs::read_to_string(response_file).unwrap();

        // We will parse each text block into a KatTest
        let mut kat_tests: Vec<KatTest> = Vec::new();

        // Iterate through each KAT block
        for block in contents.split("\n\n").skip(1) {
            // Ignore empty blocks
            if block.is_empty() {
                continue;
            }

            // Parse the block and create a KatTest from it
            let mut kat_data: KatTest = KatTest::default();
            parse_kat_test(block, &mut kat_data);
            kat_tests.push(kat_data);
        }

        kat_tests
    }

    pub fn test_kat<Fq: FqTrait>(sqisign: &Sqisign<Fq>, test: &KatTest)
    where
        [(); Fq::ENCODED_LENGTH]: Sized,
    {
        // Parse sm into the signature and message
        let (sig, msg) = test.sm.split_at(test.sm.len() - test.mlen);

        // Ensure the parsed message matches the one from sm
        assert!(msg == test.msg);

        // Ensure the signature passes!
        println!("Failed on KAT {}", test.count);
        assert!(sqisign.verify(msg, sig, &test.pk_bytes));
    }
}

#[cfg(test)]
mod test_sqisign_kat_one {
    #[cfg(feature = "kat-tests")]
    use isogeny::protocols::sqisign_parameters::SQISIGN_I;
    #[cfg(feature = "kat-tests")]
    use rayon::prelude::*;

    #[test]
    #[cfg(feature = "kat-tests")]
    fn test_kat_values() {
        let kat_tests = super::test_helpers::parse_kat_tests("PQCsignKAT_353_SQIsign_lvl1");
        kat_tests
            .par_iter()
            .for_each(|kat_test| super::test_helpers::test_kat(&SQISIGN_I, &kat_test));
    }
}

#[cfg(test)]
mod test_sqisign_kat_three {
    #[cfg(feature = "kat-tests")]
    use isogeny::protocols::sqisign_parameters::SQISIGN_III;
    #[cfg(feature = "kat-tests")]
    use rayon::prelude::*;

    #[test]
    #[cfg(feature = "kat-tests")]
    fn test_kat_values() {
        let kat_tests = super::test_helpers::parse_kat_tests("PQCsignKAT_529_SQIsign_lvl3");
        kat_tests
            .par_iter()
            .for_each(|kat_test| super::test_helpers::test_kat(&SQISIGN_III, &kat_test));
    }
}

#[cfg(test)]
mod test_sqisign_kat_five {
    #[cfg(feature = "kat-tests")]
    use isogeny::protocols::sqisign_parameters::SQISIGN_V;
    #[cfg(feature = "kat-tests")]
    use rayon::prelude::*;

    #[test]
    #[cfg(feature = "kat-tests")]
    fn test_kat_values() {
        let kat_tests = super::test_helpers::parse_kat_tests("PQCsignKAT_701_SQIsign_lvl5");
        kat_tests
            .par_iter()
            .for_each(|kat_test| super::test_helpers::test_kat(&SQISIGN_V, &kat_test));
    }
}
