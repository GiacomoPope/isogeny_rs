use std::marker::PhantomData;

use fp2::fq::Fq as FqTrait;

use crate::elliptic::{
    basis::BasisX, curve::Curve, point::PointX, three_isogeny_chain::three_isogeny_chain,
    two_isogeny_chain::two_isogeny_chain,
};

use rand::TryRngCore;

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

/// A SIDH public key contains the codomain of phi_i : E0 -> E0/<Ki> as well
/// as the images phi_i(Pj), phi_i(!j), phi_i(Pj - Qj) for (i, j) = (2, 3) and (3, 2)
/// for Alice and Bob respectively.
#[derive(Clone, Copy, Debug)]
pub struct SidhPublicKey<Fq: FqTrait> {
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
/// of length N. We additionally keep track of the exponent ea for the final isogeny.
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
    fn sample_secret_key<R: TryRngCore>(rng: &mut R) -> [u8; N] {
        let mut scalar: [u8; N] = [0; N];
        rng.try_fill_bytes(&mut scalar).unwrap();

        scalar
    }

    /// Return the domain E0 : y^2 = x^3 + 6x^2 + x
    pub fn starting_curve() -> Curve<Fq> {
        let A = Fq::from_i32(6);
        Curve::new(&A)
    }

    pub fn keygen_alice<R: TryRngCore>(
        self,
        rng: &mut R,
    ) -> (SidhPublicKey<Fq>, SidhAlicePrivateKey<Fq, N>) {
        // Sample a secret key, which is an array of bytes used as a scalar to
        // generate a kernel
        let scalar = Self::sample_secret_key(rng);

        // The domain E0 : y^2 = x^3 + 6x^2 + x
        let E = Self::starting_curve();

        // Compute the kernel K2 = P2 + [s]Q2
        let kernel = E.three_point_ladder(&self.two_torsion, &scalar, N << 3);

        // Compute phi_2 : E0 -> E/<K2> and phi(P3), phi(Q3), phi(P3 - Q3)
        let mut three_torsion_img = self.three_torsion.to_array();
        let codomain = two_isogeny_chain(&E, &kernel, self.ea, &mut three_torsion_img);

        // Package the data above into public and private keys
        let public_key = SidhPublicKey::new(&codomain, three_torsion_img);
        let secret_key = SidhAlicePrivateKey::new(self.ea, scalar);
        (public_key, secret_key)
    }

    pub fn keygen_bob<R: TryRngCore>(
        self,
        rng: &mut R,
    ) -> (SidhPublicKey<Fq>, SidhBobPrivateKey<Fq, N>) {
        // Sample a secret key, which is an array of bytes used as a scalar to
        // generate a kernel
        let scalar = Self::sample_secret_key(rng);

        // The domain E0 : y^2 = x^3 + 6x^2 + x
        let E = Self::starting_curve();

        // Compute the kernel K3 = P3 + [s]Q3
        let kernel = E.three_point_ladder(&self.three_torsion, &scalar, N << 3);

        // Compute phi_3 : E0 -> E/<K3> and phi(P2), phi(Q2), phi(P2 - Q2)
        let mut two_torsion_img = self.two_torsion.to_array();
        let codomain = three_isogeny_chain(&E, &kernel, self.eb, &mut two_torsion_img);

        // Package the data above into public and private keys
        let public_key = SidhPublicKey::new(&codomain, two_torsion_img);
        let secret_key = SidhBobPrivateKey::new(self.eb, scalar);
        (public_key, secret_key)
    }
}

impl<Fq: FqTrait, const N: usize> SidhAlicePrivateKey<Fq, N> {
    pub fn new(exp: usize, scalar: [u8; N]) -> Self {
        Self {
            exp,
            scalar,
            _phantom: PhantomData::default(),
        }
    }

    pub fn shared_secret(self, public_key: &SidhPublicKey<Fq>) -> Fq {
        // Extract out the codomain of phi_3 : E -> E/<K3>
        let E = public_key.E;

        // Compute the kernel K = phi_3(P2) + [s] phi_3(Q2)
        let kernel = E.three_point_ladder(&public_key.basis_img, &self.scalar, N << 3);

        // Compute the codomain of E0 -> E3 -> ES = (E / <K3>) / <K>
        let codomain = two_isogeny_chain(&E, &kernel, self.exp, &mut []);
        codomain.j_invariant()
    }
}

impl<Fq: FqTrait, const N: usize> SidhBobPrivateKey<Fq, N> {
    pub fn new(exp: usize, scalar: [u8; N]) -> Self {
        Self {
            exp,
            scalar,
            _phantom: PhantomData::default(),
        }
    }

    pub fn shared_secret(self, public_key: &SidhPublicKey<Fq>) -> Fq {
        // Extract out the codomain of phi_2 : E -> E/<K2>
        let E = public_key.E;

        // Compute the kernel K = phi_2(P3) + [s] phi_2(Q3)
        let kernel = E.three_point_ladder(&public_key.basis_img, &self.scalar, N << 3);

        // Compute the codomain of E0 -> E2 -> ES = (E / <K2>) / <K>
        let codomain = three_isogeny_chain(&E, &kernel, self.exp, &mut []);
        codomain.j_invariant()
    }
}

impl<Fq: FqTrait> SidhPublicKey<Fq> {
    pub fn new(E: &Curve<Fq>, basis_img: [PointX<Fq>; 3]) -> Self {
        Self {
            E: *E,
            basis_img: BasisX::from_array(basis_img),
        }
    }
}
