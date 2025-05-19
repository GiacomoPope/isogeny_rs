#[cfg(test)]
mod test_sqisign_one {
    use isogeny::protocols::sqisign_parameters::SQISIGN_I;

    static MSG: &str = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
    static PK: &str = "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B";
    static SIG: &str = "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F0455D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC010001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471D86EFF3447018ADB0A6551EE8322AB30010202";

    #[test]
    fn test_verification() {
        // For now, just convert the hex strings to bytes within the test, will be handled later
        // by proper KAT stuff I suppose.
        let msg: &[u8] = &hex::decode(MSG).unwrap();
        let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
        let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
        assert!(SQISIGN_I.verify(msg, sig_bytes, pk_bytes))
    }

    #[test]
    fn test_verification_fail() {
        let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
        let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
        assert!(!SQISIGN_I.verify(b"", sig_bytes, pk_bytes))
    }
}

#[cfg(test)]
mod test_sqisign_three {
    use isogeny::protocols::sqisign_parameters::SQISIGN_III;

    static MSG: &str = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
    static PK: &str = "C32377D6F6D70729884A7F6877EF4791E35D21F751A3E96DE23F9A7A3C01BCD8A5F146DC19E4E2AC63007457F97D8A40EE84AEE7564CA9A7FBE6200FD3E5E55901BFC60EB25C50D39F5C91C96510556BAA22028DF76360841721A601D65E8D0F06";
    static SIG: &str = "0868CFBF275B8E7B19BF597D658D62CC913B9B2933E30A297288FBE687F6F6B8AC8AF7AA007F191386BB1A203CDDBC2BDB42792D05DA69A4507073D12B0BDC47E2B36BC4BA45C68791918281E578F2DC14294504726DCD4CA4C4565FBB89A12800048C7B84746A2CBD8247248E248B70B51AE91994957857692A028D8F5CABABFC91E4BF1C5D350219A0189C57DE4A7710D29E0364C79B2188449EC0397359430D594C7B5980CC67551933A902D3C11F0FBD6DC39711D3E1F501159EE7FB85CE81B4CE24E1016006567DF469315D513E73F69F6301664E6449AF9DCEB4000D15";

    #[test]
    fn test_verification() {
        // For now, just convert the hex strings to bytes within the test, will be handled later
        // by proper KAT stuff I suppose.
        let msg: &[u8] = &hex::decode(MSG).unwrap();
        let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
        let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
        assert!(SQISIGN_III.verify(msg, sig_bytes, pk_bytes))
    }

    #[test]
    fn test_verification_fail() {
        let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
        let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
        assert!(!SQISIGN_III.verify(b"", sig_bytes, pk_bytes))
    }
}

#[cfg(test)]
mod test_sqisign_five {
    use isogeny::protocols::sqisign_parameters::SQISIGN_V;

    static MSG: &str = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
    static PK: &str = "86FFA3B0F73D55A64D13C6F89F28D75FD17C5E2368E1D451127C16D1A97CDB440E20333A233AD2F8E4D70187C8AE31602049ADE949A87F95E79DA4C456F5D400B2485A96D04708A2F30046812B8D65A3BFBFDED0DD6563462F9E2BCE760CD753CAE8471BEC7049EF28FFEFE859C15DAC49DB959AEE99842D97A380A70DD7330106";
    static SIG: &str = "6B8EF5D7689A1EA1CFCE9C6F7495E309E9D1D1B03E61CD97088E679C4901D0B6B6D38217F4AED6C44949B41F9AF80B43E84D0C91BDB1D00E06957BEBF30A58012AD01E52CF7906CE197AD06696F7FCF756908EA980549E7C215D089BDE7117799F628817A1B9C8FB7FEBFF7E9D9B776142460CFAAFC97D48A57E09E0DA378401000229CC8E1B94E1F2F8AFDC42066BEACE076E3E70DD01F90C4D01DAC17BEC58743532848D438A87A574D9DB940C17236AE3566281E27A99EFE5EE26E05B88A1D610A80B3AF38267D845C7FE330F199B43794A9B2E14846924127366B8F6A1F0F24D3C4B54D79DBB61B098BF32D98EA8819F7BE4A5FFBA29E88B1A996C6CDFD32B048BC2ACFFA28870181447FCC8B6F97B63C47CB013C6F3D84CBD07619A5C355B000911";

    #[test]
    fn test_verification() {
        // For now, just convert the hex strings to bytes within the test, will be handled later
        // by proper KAT stuff I suppose.
        let msg: &[u8] = &hex::decode(MSG).unwrap();
        let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
        let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
        assert!(SQISIGN_V.verify(msg, sig_bytes, pk_bytes))
    }

    #[test]
    fn test_verification_fail() {
        let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
        let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
        assert!(!SQISIGN_V.verify(b"", sig_bytes, pk_bytes))
    }
}
