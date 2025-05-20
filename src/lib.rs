// We favour using terms like E for an elliptic curve, or A for its
// Montgomery coefficient, as it is standard in the literature.
#![allow(non_snake_case)]
// We include these so we can have things like
// fn encode(self) -> [u8; Self::ENCODED_LENGTH];
// defined within the Fq trait
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

pub mod elliptic;
pub mod fields;
pub mod protocols;
pub mod theta;
pub mod utilities;
