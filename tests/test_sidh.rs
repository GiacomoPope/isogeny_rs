#[cfg(test)]
mod test_sidh {
    use fp2::fq::Fq;
    use isogeny::protocols::sidh_parameters::SIDH_434;
    use isogeny::utilities::drng::DRNG;

    #[test]
    fn test_key_exchange() {
        let mut rng = DRNG::from_seed("sidh_key_exchange".as_bytes());

        // Test 5 different public key pairs.
        for _ in 0..5 {
            let (alice_pub, alice_priv) = SIDH_434.keygen_alice(&mut rng);
            let (bob_pub, bob_priv) = SIDH_434.keygen_bob(&mut rng);

            let alice_secret = alice_priv.shared_secret(&bob_pub);
            let bob_secret = bob_priv.shared_secret(&alice_pub);
            assert!(alice_secret.equals(&bob_secret) == u32::MAX);
        }
    }
}
