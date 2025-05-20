use ::std::time::Instant;
use isogeny::protocols::sqisign_parameters::SQISIGN_I;

const MSG: &str = "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";
const PK: &str = "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B2029550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A010B";
const SIG: &str = "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F0455D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC010001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471D86EFF3447018ADB0A6551EE8322AB30010202";

fn main() {
    // For now, just convert the hex strings to bytes within the test, will be handled later
    // by proper KAT stuff I suppose.
    let msg: &[u8] = &hex::decode(MSG).unwrap();
    let pk_bytes: &[u8] = &hex::decode(PK).unwrap();
    let sig_bytes: &[u8] = &hex::decode(SIG).unwrap();
    let timer = Instant::now();
    let success = SQISIGN_I.verify(msg, sig_bytes, pk_bytes) as u32;
    let t1 = timer.elapsed();
    println!("verify: {}, time: {:?}", success, t1);
}
