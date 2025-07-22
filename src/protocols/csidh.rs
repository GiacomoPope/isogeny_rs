use std::{error::Error, usize};
use std::fmt::Display;


use fp2::traits::Fp as FqTrait;

use crate::{
    elliptic::{curve::Curve, point::PointX},
    polynomial_ring::poly::Polynomial,
    utilities::bn::{Bn, mul_bn_by_u64_vartime, bn_compare_vartime},
};

use rand_core::{CryptoRng, RngCore};

#[derive(Clone, Copy, Debug)]
pub struct CsidhParameters<const NUM_ELLS: usize> {
    pub max_exponent: usize,
    pub two_cofactor: usize,
    pub primes: [u64; NUM_ELLS],
    pub four_sqrt_p: Bn,
}

pub struct Csidh<Fp: FqTrait, const NUM_ELLS: usize> {
    max_exponent: usize,
    two_cofactor: usize,
    base: Fp,
    primes: [u64; NUM_ELLS],
    four_sqrt_p: Bn,
}

pub struct CsidhPrivateKey<const NUM_ELLS: usize> {
    e: [u32; NUM_ELLS],  // secret degree
    d: [bool; NUM_ELLS], // secret direction
}

#[derive(Clone, Copy, Debug)]
pub struct CsidhPublicKey<Fp: FqTrait> {
    pub A: Fp,
}

#[derive(Debug)]
pub enum CsidhError {
    PublicKeyVerificationError,
}

impl Display for CsidhError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CsidhError::PublicKeyVerificationError => f.write_str(
                "CSIDH PublicKey verfication faild (EllipticCurve is not supersingular.)",
            ),
        }
    }
}

impl Error for CsidhError {}

impl<Fp: FqTrait> PartialEq for CsidhPublicKey<Fp> {
    fn eq(&self, other: &Self) -> bool {
        self.A.equals(&other.A) == u32::MAX
    }
}

impl<Fp: FqTrait + std::fmt::Debug, const NUM_ELLS: usize> Csidh<Fp, NUM_ELLS> {
    pub const fn new(params: &CsidhParameters<NUM_ELLS>) -> Self {
        Self {
            max_exponent: params.max_exponent,
            two_cofactor: params.two_cofactor,
            base: Fp::ZERO,
            primes: params.primes,
            four_sqrt_p: params.four_sqrt_p,
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
        // get a random point
        let x = Fp::rand(rng);
        let P = PointX::from_x_coord(&x);

        // projective curve: x^3C+4Ax^2−2x^2C+XC​
        let tmp = *A24*x;
        let four_Ax = tmp+tmp+tmp+tmp;
        let tmp = *C24*x;
        let Cxx = tmp*x;
        let two_Cx = tmp+tmp;
        
        // if RHS is a sqrt, it is on the curve
        let RHS = (Cxx+four_Ax-two_Cx+(*C24))*x;
        (P, RHS.is_square() == u32::MAX)
    }

    fn sample_secret_key<R: CryptoRng + RngCore>(&self, rng: &mut R) -> CsidhPrivateKey<NUM_ELLS> {
        let mut e = [0u32; NUM_ELLS];
        let mut d = [false; NUM_ELLS];

        for i in 0..NUM_ELLS {
            // sample exponent between 0 and max_exponent
            e[i] = rng.next_u32();
            e[i] >>= 27; // shift to increase probabily to hit the range
            while e[i] > self.max_exponent as u32 {
                e[i] = rng.next_u32() ;
                e[i] >>= 27;
            }

            // sample direction
            d[i] = (rng.next_u32() % 2) == 0;
        }

        CsidhPrivateKey { e, d }
    }

    //
    //      ACTION
    // 
    pub fn action<R: CryptoRng + RngCore>(
        &self,
        public_key: &CsidhPublicKey<Fp>,
        private_key: &CsidhPrivateKey<NUM_ELLS>,
        rng: &mut R,
    ) -> CsidhPublicKey<Fp> {
        let mut A24 = public_key.A + Fp::TWO;
        let mut C24 = Fp::FOUR;

        let mut sk_e = private_key.e;
        let sk_d = private_key.d;

        let mut done: u32 = sk_e.iter().sum();
        while done != 0 {
            let (mut P, direction) = self.rand_point(&A24, &C24, rng);

            // clear 2^e cofactor
            Curve::<Fp>::xdbl_proj_iter(&A24, &C24, &mut P, self.two_cofactor);

            for i in (0..NUM_ELLS).rev() {
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
                let mut img_points = [P];
                Curve::<Fp>::velu_prime_isogeny_proj::<Polynomial<Fp>>(
                    &mut A24,
                    &mut C24,
                    &K,
                    degree as usize,
                    &mut img_points,
                );

                P = img_points[0];

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



    //  
    //      VERIFICATION
    //

    /// recursively computes [p+1/l] for all l (primes)
    fn cofactor_multiples(&self, P : &mut [PointX<Fp>;NUM_ELLS], A24: &Fp , C24: &Fp, lower: usize, upper: usize) {
        if upper - lower == 1 {
            return;
        }

        let mid = lower + (upper - lower + 1) / 2;

        let mut cl = vec![1 as u64];
        let mut cu = vec![1 as u64];

        for i in lower..mid{
            cu = mul_bn_by_u64_vartime(&cu, self.primes[i]);
        }for i in mid..upper{
            cl = mul_bn_by_u64_vartime(&cl, self.primes[i]);
        }

        
        P[mid] = Curve::<Fp>::xmul_proj_bn_vartime(A24, C24, &mut P[lower], &cu[..]);
        P[lower] = Curve::<Fp>::xmul_proj_bn_vartime(A24, C24, &mut P[lower], &cl[..]);

        self.cofactor_multiples(P, A24, C24, lower, mid);
        self.cofactor_multiples(P, A24, C24, mid, upper);
    }


    pub fn verify<R: CryptoRng + RngCore>(&self, public_key: &CsidhPublicKey<Fp>, rng: &mut R) -> bool {
        let A24 = public_key.A + Fp::TWO;
        let C24 = Fp::FOUR;

        let fsqrtp = &self.four_sqrt_p.limbs[..self.four_sqrt_p.len];

        loop {
            let mut P : [PointX<Fp>; NUM_ELLS] = [PointX::INFINITY ; NUM_ELLS];

            let (tmp, _) = self.rand_point(&A24, &C24, rng);
            P[0] = tmp;
            
            Curve::<Fp>::xdbl_proj_iter(&A24, &C24, &mut P[0], self.two_cofactor);

            let mut order = vec![1];

            self.cofactor_multiples(&mut P, &A24, &C24, 0, NUM_ELLS);

            for i in (0..NUM_ELLS).rev() {
                if P[i].is_zero() == u32::MAX {
                    continue;
                }

                // P should now have order l, so we verify that by checking [l]P = 0
                P[i] = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &P[i], self.primes[i]);

                if P[i].is_zero() == 0 {
                    return false;
                }

                // if the order of out starting Point is > 4sqrt(p), the curve must be supersingular
                order = mul_bn_by_u64_vartime(&order[..], self.primes[i]);
                match bn_compare_vartime(&order[..], fsqrtp) {
                    std::cmp::Ordering::Greater => return true,
                    _ => {},
                }
            }
        }
    }

    pub fn keygen<R: CryptoRng + RngCore>(
        self,
        rng: &mut R,
    ) -> (CsidhPrivateKey<NUM_ELLS>, CsidhPublicKey<Fp>) {
        let sk = self.sample_secret_key(rng);

        let pk = self.action(&self.starting_curve(), &sk, rng);

        (sk, pk)
    }

    pub fn derive_shared_key<R: CryptoRng + RngCore>(
        self,
        public_key: &CsidhPublicKey<Fp>,
        private_key: &CsidhPrivateKey<NUM_ELLS>,
        rng: &mut R,
    ) -> Result<CsidhPublicKey<Fp>, CsidhError> {
        match self.verify(public_key, rng) {
            true => Ok(self.action(public_key, private_key, rng)),
            false => Err(CsidhError::PublicKeyVerificationError),
        }
    }
}
