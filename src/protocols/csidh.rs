use std::{iter::OnceWith, os::unix::raw::uid_t};

use fp2::traits::Fp as FqTrait;

use crate::elliptic::{curve::Curve, point::PointX};

use rand_core::{CryptoRng, RngCore};

#[derive(Clone, Copy, Debug)]
pub struct CsidhParameters<const COUNT: usize> {
    pub num_primes: usize,
    pub max_exponent: usize,
    pub two_cofactor: usize,
    pub primes: [u64; COUNT], 
}


pub struct Csidh<Fp: FqTrait, const COUNT: usize> {
    num_primes: usize,
    max_exponent: usize,
    two_cofactor: usize,
    base: Fp,
    primes: [u64; COUNT],  
}

pub struct CsidhPrivateKey<const COUNT: usize> {
    e : [u64; COUNT]
}

#[derive(Clone, Copy, Debug)]
pub struct CsidhPublicKey<Fp: FqTrait> {
    pub A : Fp
}

impl<Fp: FqTrait> PartialEq for CsidhPublicKey<Fp> {
    fn eq(&self, other : &Self) -> bool {
        self.A.equals(&other.A) == u32::MAX
    }
}

impl<Fp: FqTrait, const COUNT: usize> Csidh<Fp, COUNT> {
    pub const fn new(params: &CsidhParameters<COUNT>) -> Self {
        Self { 
            num_primes: params.num_primes,
            max_exponent: params.max_exponent,
            two_cofactor: params.two_cofactor,
            base: Fp::ZERO,             // Todo: could also be a parameter
            primes: params.primes,
        }
    }

    pub const fn starting_curve(&self) ->  CsidhPublicKey<Fp> {
        CsidhPublicKey{ A: self.base }
    }

    fn sample_secret_key<R: CryptoRng + RngCore>(&self, rng: &mut R) -> CsidhPrivateKey<COUNT> {

        let mut e = [0u64; COUNT];

        for i in 0..COUNT {
            e[i] = rng.next_u64();
            e[i] >>= 59;
            while e[i] > self.max_exponent as u64 {
                e[i] = rng.next_u64();
                e[i] >>= 59;
            }
        }

        CsidhPrivateKey { e }
    }

    fn action<R: CryptoRng + RngCore>(
        &self, 
        public_key: &CsidhPublicKey<Fp>, 
        private_key: &CsidhPrivateKey<COUNT>,
        rng: &mut R) 
        -> CsidhPublicKey<Fp> {
        
        let mut A24 = public_key.A + Fp::TWO;
        let mut C24 = Fp::FOUR;
        let mut sk = private_key.e;

        let mut done : u64 = sk.iter().sum();
        while done != 0{
            // get point
            // can be optimized using projective curves
            // to sample the point
            let curve = Curve::curve_from_A24_proj(&A24, &C24);
            let (mut P, mut ok) = curve.rand_pointX(rng);
            while ok == 0 {
                (P, ok) = curve.rand_pointX(rng);
            }

            // clear cofactor
            Curve::<Fp>::xdbl_proj_iter(&A24, &C24, &mut P, self.two_cofactor);


            for i in (0..self.num_primes).rev() {
                let secret_e = sk[i];
                let degree = self.primes[i];

                if secret_e == 0 {
                    continue;
                }
                   
                // get kernel from point
                let mut K = P;
                for ell in self.primes.iter() {
                    if *ell == degree { continue;}
                    K = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &K, *ell);
                }


                if K.is_zero() == u32::MAX {
                    continue;
                }

                // let K_ = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &K, degree);
                // assert!(K_.is_zero() == u32::MAX); 

                let mut image_points = [P];
                // compute isogeny
                Curve::<Fp>::velu_odd_isogeny_proj(&mut A24, &mut C24, &K, degree as usize, &mut image_points);

                P = image_points[0];
                // mark step as done
                sk[i] -= 1;
            }


            done = sk.iter().sum();
        }

    
        CsidhPublicKey{ A: Curve::<Fp>::curve_from_A24_proj(&A24, &C24).A }
    }


    pub fn keygen<R: CryptoRng + RngCore>(self, 
        rng: &mut R) -> (CsidhPrivateKey<COUNT>, CsidhPublicKey<Fp>) {
        let sk = self.sample_secret_key(rng);

        let pk = self.action(&self.starting_curve(), &sk, rng);

        (sk, pk)
    }

    pub fn derive_shared_key<R: CryptoRng + RngCore>(
        self, 
        public_key: &CsidhPublicKey<Fp>, 
        private_key: &CsidhPrivateKey<COUNT>,
        rng: &mut R,) -> CsidhPublicKey<Fp> {
        
        // todo: verify

        self.action(public_key, private_key, rng)
    }
}

