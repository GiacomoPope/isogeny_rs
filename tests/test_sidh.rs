#[cfg(test)]
mod test_sidh {
    use fp2::traits::Fp as _;
    use isogeny::protocols::sidh_parameters::SIDH_434;
    use isogeny::utilities::test_utils::drng::DRNG;

    #[test]
    fn test_key_exchange() {
        let mut rng = DRNG::from_seed("sidh_key_exchange".as_bytes());

        let (alice_pub, alice_priv) = SIDH_434.keygen_alice(&mut rng);
        let (bob_pub, bob_priv) = SIDH_434.keygen_bob(&mut rng);

        let alice_secret = alice_priv.shared_secret(&bob_pub).unwrap();
        let bob_secret = bob_priv.shared_secret(&alice_pub).unwrap();
        assert!(alice_secret.equals(&bob_secret) == u32::MAX);
    }
}
