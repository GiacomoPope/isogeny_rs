use std::marker::PhantomData;

use fp2::fq::Fq as FqTrait;

use crate::elliptic::{basis::BasisX, curve::Curve, point::PointX};

use rand_core::{CryptoRng, RngCore};

/// Public parameters used for a SIDH key exchange.
/// `two_torsion` is the x-only basis for the \[2^ea\] torsion <P2, Q2>
/// `three_torsion` is the x-only basis for the \[3^eb\] torsion <P3, Q3>
/// `N` is the number of bytes needed to represent the secret scalar.
#[derive(Clone, Copy, Debug)]
pub struct SidhParameters<Fq: FqTrait, const N: usize> {
    ea: usize,
    eb: usize,
    two_torsion: BasisX<Fq>,
    three_torsion: BasisX<Fq>,
}

/// Alice's SIDH public key contains the codomain of phi_2 : E0 -> E0/<K2> with
/// K2 = P2 + \[s\] * Q2, with P2, Q2 as public parameters and `s` as the scalar
/// stored in `SidhAlicePrivateKey`. Additionally, it contains the images of
/// phi_2(P3), phi_2(Q3) and phi_2(P3 - Q3)
#[derive(Clone, Copy, Debug)]
pub struct SidhAlicePublicKey<Fq: FqTrait> {
    E: Curve<Fq>,
    basis_img: BasisX<Fq>,
}

/// Bob's SIDH public key contains the codomain of phi_3 : E0 -> E0/<K3> with
/// K3 = P3 + \[s\] * Q3, with P3, Q3 as public parameters and `s` as the scalar
/// stored in `SidhBobPrivateKey`. Additionally, it contains the images of
/// phi_3(P2), phi_3(Q2) and phi_3(P2 - Q2)
#[derive(Clone, Copy, Debug)]
pub struct SidhBobPublicKey<Fq: FqTrait> {
    E: Curve<Fq>,
    basis_img: BasisX<Fq>,
}

/// SIDH pivate key for Alice, simply a scalar represented as an array of bytes
/// of length N. We additionally keep track of the exponent ea for the final isogeny.
#[derive(Clone, Copy, Debug)]
pub struct SidhAlicePrivateKey<Fq: FqTrait, const N: usize> {
    exp: usize,
    scalar: [u8; N],
    _phantom: PhantomData<Fq>,
}

/// SIDH pivate key for Alice, simply a scalar represented as an array of bytes
/// of length N. We additionally keep track of the exponent eb for the final isogeny.
#[derive(Clone, Copy, Debug)]
pub struct SidhBobPrivateKey<Fq: FqTrait, const N: usize> {
    exp: usize,
    scalar: [u8; N],
    _phantom: PhantomData<Fq>,
}

impl<Fq: FqTrait, const N: usize> SidhParameters<Fq, N> {
    pub const fn new(
        ea: usize,
        eb: usize,
        two_torsion: &BasisX<Fq>,
        three_torsion: &BasisX<Fq>,
    ) -> Self {
        Self {
            ea,
            eb,
            two_torsion: *two_torsion,
            three_torsion: *three_torsion,
        }
    }

    // TODO: find a better way to work with randomness.
    /// Sample a secret scalar represented as N bytes.
    fn sample_secret_key<R: CryptoRng + RngCore>(rng: &mut R) -> [u8; N] {
        let mut scalar: [u8; N] = [0; N];
        rng.fill_bytes(&mut scalar);

        scalar
    }

    /// Return the domain E0 : y^2 = x^3 + 6x^2 + x
    pub fn starting_curve() -> Curve<Fq> {
        let A = Fq::from_i32(6);
        Curve::new(&A)
    }

    pub fn keygen_alice<R: CryptoRng + RngCore>(
        self,
        rng: &mut R,
    ) -> (SidhAlicePublicKey<Fq>, SidhAlicePrivateKey<Fq, N>) {
        // Sample a secret key, which is an array of bytes used as a scalar to
        // generate a kernel
        let scalar = Self::sample_secret_key(rng);

        // The domain E0 : y^2 = x^3 + 6x^2 + x
        let E = Self::starting_curve();

        // Compute the kernel K2 = P2 + [s]Q2
        let kernel = E.three_point_ladder(&self.two_torsion, &scalar, N << 3);

        // Compute phi_2 : E0 -> E/<K2> and phi_2(P3), phi_2(Q3), phi_2(P3 - Q3)
        let mut three_torsion_img = self.three_torsion.to_array();
        let codomain = E.two_isogeny_chain(&kernel, self.ea, &mut three_torsion_img);

        // Package the data above into public and private keys
        let public_key = SidhAlicePublicKey::new(&codomain, &three_torsion_img);
        let secret_key = SidhAlicePrivateKey::new(self.ea, scalar);
        (public_key, secret_key)
    }

    pub fn keygen_bob<R: CryptoRng + RngCore>(
        self,
        rng: &mut R,
    ) -> (SidhBobPublicKey<Fq>, SidhBobPrivateKey<Fq, N>) {
        // Sample a secret key, which is an array of bytes used as a scalar to
        // generate a kernel
        let scalar = Self::sample_secret_key(rng);

        // The domain E0 : y^2 = x^3 + 6x^2 + x
        let E = Self::starting_curve();

        // Compute the kernel K3 = P3 + [s]Q3
        let kernel = E.three_point_ladder(&self.three_torsion, &scalar, N << 3);

        // Compute phi_3 : E0 -> E/<K3> and phi_3(P2), phi_3(Q2), phi_3(P2 - Q2)
        let mut two_torsion_img = self.two_torsion.to_array();
        let codomain = E.three_isogeny_chain(&kernel, self.eb, &mut two_torsion_img);

        // Package the data above into public and private keys
        let public_key = SidhBobPublicKey::new(&codomain, &two_torsion_img);
        let secret_key = SidhBobPrivateKey::new(self.eb, scalar);
        (public_key, secret_key)
    }
}

impl<Fq: FqTrait, const N: usize> SidhAlicePrivateKey<Fq, N> {
    pub fn new(exp: usize, scalar: [u8; N]) -> Self {
        Self {
            exp,
            scalar,
            _phantom: PhantomData,
        }
    }

    pub fn shared_secret(self, public_key: &SidhBobPublicKey<Fq>) -> Fq {
        // Extract out the codomain of phi_3 : E -> E/<K3>
        let E = public_key.E;

        // Compute the kernel K = phi_3(P2) + [s] phi_3(Q2)
        let kernel = E.three_point_ladder(&public_key.basis_img, &self.scalar, N << 3);

        // Compute the codomain of E0 -> E3 -> ES = (E / <K3>) / <K>
        let codomain = E.two_isogeny_chain(&kernel, self.exp, &mut []);
        codomain.j_invariant()
    }
}

impl<Fq: FqTrait, const N: usize> SidhBobPrivateKey<Fq, N> {
    pub fn new(exp: usize, scalar: [u8; N]) -> Self {
        Self {
            exp,
            scalar,
            _phantom: PhantomData,
        }
    }

    pub fn shared_secret(self, public_key: &SidhAlicePublicKey<Fq>) -> Fq {
        // Extract out the codomain of phi_2 : E -> E/<K2>
        let E = public_key.E;

        // Compute the kernel K = phi_2(P3) + [s] phi_2(Q3)
        let kernel = E.three_point_ladder(&public_key.basis_img, &self.scalar, N << 3);

        // Compute the codomain of E0 -> E2 -> ES = (E / <K2>) / <K>
        let codomain = E.three_isogeny_chain(&kernel, self.exp, &mut []);
        codomain.j_invariant()
    }
}

impl<Fq: FqTrait> SidhAlicePublicKey<Fq> {
    pub fn new(E: &Curve<Fq>, basis_img: &[PointX<Fq>]) -> Self {
        Self {
            E: *E,
            basis_img: BasisX::from_slice(basis_img),
        }
    }
}

impl<Fq: FqTrait> SidhBobPublicKey<Fq> {
    pub fn new(E: &Curve<Fq>, basis_img: &[PointX<Fq>]) -> Self {
        Self {
            E: *E,
            basis_img: BasisX::from_slice(basis_img),
        }
    }
}
