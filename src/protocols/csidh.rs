
use fp2::traits::Fp as FqTrait;

use crate::elliptic::{curve::Curve, point::PointX};

use rand_core::{CryptoRng, RngCore};

#[derive(Clone, Copy, Debug)]
pub struct CsidhParameters<const COUNT: usize> {
    pub max_exponent: usize,
    pub two_cofactor: usize,
    pub primes: [u64; COUNT],
}

pub struct Csidh<Fp: FqTrait, const COUNT: usize> {
    max_exponent: usize,
    two_cofactor: usize,
    base: Fp,
    primes: [u64; COUNT],
}

pub struct CsidhPrivateKey<const COUNT: usize> {
    e: [u64; COUNT],  // secret degree
    d: [bool; COUNT], // secret direction
}

#[derive(Clone, Copy, Debug)]
pub struct CsidhPublicKey<Fp: FqTrait> {
    pub A: Fp,
}

impl<Fp: FqTrait> PartialEq for CsidhPublicKey<Fp> {
    fn eq(&self, other: &Self) -> bool {
        self.A.equals(&other.A) == u32::MAX
    }
}

impl<Fp: FqTrait, const COUNT: usize> Csidh<Fp, COUNT> {

    pub const fn new(params: &CsidhParameters<COUNT>) -> Self {
        Self {
            max_exponent: params.max_exponent,
            two_cofactor: params.two_cofactor,
            base: Fp::ZERO,
            primes: params.primes,
        }
    }

    // Returns the starting curve
    pub const fn starting_curve(&self) -> CsidhPublicKey<Fp> {
        CsidhPublicKey { A: self.base }
    }

    /// get a random point that either on the curve or on its twist
    fn rand_point<R: CryptoRng + RngCore>(
        &self,
        A24: &Fp,
        C24: &Fp,
        rng: &mut R,
    ) -> (PointX<Fp>, bool) {
        let curve = Curve::<Fp>::curve_from_A24_proj(&A24, &C24);
        let P = PointX::rand_point(rng);

        // to check if the point is on the curve,
        // or on the twist
        let (_, twist) = curve.lift_pointx(&P);

        (P, twist == 0)
    }

    fn sample_secret_key<R: CryptoRng + RngCore>(&self, rng: &mut R) -> CsidhPrivateKey<COUNT> {
        let mut e = [0u64; COUNT];
        let mut d = [false; COUNT];

        for i in 0..COUNT {
            // sample exponent between 0 and max_exponent
            e[i] = rng.next_u64();
            e[i] >>= 59; // shift to increase probabily to hit the range
            while e[i] > self.max_exponent as u64 {
                e[i] = rng.next_u64();
                e[i] >>= 59;
            }

            // sample direction
            d[i] = (rng.next_u32() % 2) == 0;
        }

        CsidhPrivateKey { e, d }
    }

    fn action<R: CryptoRng + RngCore>(
        &self,
        public_key: &CsidhPublicKey<Fp>,
        private_key: &CsidhPrivateKey<COUNT>,
        rng: &mut R,
    ) -> CsidhPublicKey<Fp> {
        let mut A24 = public_key.A + Fp::TWO;
        let mut C24 = Fp::FOUR;

        let mut sk_e = private_key.e;
        let sk_d = private_key.d;

        let mut done: u64 = sk_e.iter().sum();
        while done != 0 {
            let (mut P, direction) = self.rand_point(&A24, &C24, rng);

            // clear 2^e cofactor
            Curve::<Fp>::xdbl_proj_iter(&A24, &C24, &mut P, self.two_cofactor);

            for i in (0..COUNT).rev() {
                let secret_e = sk_e[i];
                let secret_d = sk_d[i];
                let degree = self.primes[i];

                // check if the current degree is part of the secrete key
                if secret_e == 0 {
                    continue;
                }

                // check if the sampled point is for
                // the correct direction
                if secret_d != direction {
                    continue;
                }

                // get kernel from point
                // we do this by removing every ell, but the current degree
                // (this can be optimized way more)
                let mut K = P;
                for ell in self.primes.iter() {
                    if *ell == degree {
                        continue;
                    }
                    K = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &K, *ell);
                }

                // check if the kernel is of the correct degree
                // if not we skip it and try again on an new round
                if K.is_zero() == u32::MAX {
                    continue;
                }

                // Finally compute the isogeny and push P
                let mut image_points = [P];
                Curve::<Fp>::velu_odd_isogeny_proj(
                    &mut A24,
                    &mut C24,
                    &K,
                    degree as usize,
                    &mut image_points,
                );
                P = image_points[0];

                // mark step as done
                sk_e[i] -= 1;

                // did we "exhaust" the point?
                if P.is_zero() == u32::MAX {
                    break;
                }
            }

            done = sk_e.iter().sum();
        }

        CsidhPublicKey {
            A: Curve::<Fp>::curve_from_A24_proj(&A24, &C24).A,
        }
    }

    pub fn keygen<R: CryptoRng + RngCore>(
        self,
        rng: &mut R,
    ) -> (CsidhPrivateKey<COUNT>, CsidhPublicKey<Fp>) {
        let sk = self.sample_secret_key(rng);

        let pk = self.action(&self.starting_curve(), &sk, rng);

        (sk, pk)
    }

    pub fn derive_shared_key<R: CryptoRng + RngCore>(
        self,
        public_key: &CsidhPublicKey<Fp>,
        private_key: &CsidhPrivateKey<COUNT>,
        rng: &mut R,
    ) -> CsidhPublicKey<Fp> {
        // todo: verify

        self.action(public_key, private_key, rng)
    }
}
