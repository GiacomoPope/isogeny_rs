#![allow(non_snake_case)]

#[cfg(test)]
mod test_csidh {
    use isogeny::protocols::csidh_parameters::CSIDH_512;
    use isogeny::utilities::test_utils::drng::DRNG;


    #[test]
    fn test_csidh_key_exchange() {
        let mut rng = DRNG::from_seed("csidh_key_exchange".as_bytes());
        
        let (alice_sk, alice_pk) = CSIDH_512.keygen(&mut rng);
        let (bob_sk, bob_pk) = CSIDH_512.keygen(&mut rng);

        let alice_shared = CSIDH_512.derive_shared_key(&bob_pk, &alice_sk, &mut rng);
        let bob_shared = CSIDH_512.derive_shared_key(&alice_pk, &bob_sk, &mut rng);

        println!("{:#?}", alice_shared);
        println!("{:#?}", bob_shared);
        assert!(false);
    }
}