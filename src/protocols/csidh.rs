use std::os::unix::raw::uid_t;

use fp2::traits::Fp as FqTrait;

use crate::elliptic::{curve::Curve, point::PointX};

use rand_core::{CryptoRng, RngCore};

#[derive(Clone, Copy, Debug)]
pub struct CsidhParameters<const COUNT: usize> {
    pub num_primes: usize,
    pub max_exponent: usize,
    pub two_cofactor: u64,
    pub primes: [u64; COUNT],  //Todo: we need to somehow load this from params...
}


pub struct Csidh<Fp: FqTrait, const COUNT: usize> {
    num_primes: usize,
    max_exponent: usize,
    two_cofactor: u64,
    base: Fp,
    primes: [u64; COUNT],   //Todo: we need to somehow load this from params...
}

pub struct CsidhPrivateKey<const COUNT: usize> {
    e : [u64; COUNT]
}

#[derive(Clone, Copy, Debug)]
pub struct CsidhPublicKey<Fp: FqTrait> {
    pub A : Fp
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

    fn sample_secret_key(&self) -> CsidhPrivateKey<COUNT> {
        // for testing, always use [1,1,1,...,1] as key
        CsidhPrivateKey { e: [1u64; COUNT] }
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
            let x = Fp::rand(rng);
            let mut P = PointX::from_x_coord(&x);
            
            // clear cofactor
            for _ in 0..self.two_cofactor {
                P = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &P, 2);
            }

            for i in (0..self.num_primes).rev() {
                let secret_e = private_key.e[i];
                let degree = self.primes[i];
                println!("step: {} for {}", secret_e, degree);
                
                if secret_e == 0 {
                    // remove degree from point
                    P = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &P, degree);
                    continue;
                }
                   
                // get kernel from point
                let mut K = P;
                for ell in self.primes[0..i].iter() {
                    K = Curve::<Fp>::xmul_proj_u64_vartime(&A24, &C24, &K, *ell);
                }

                if K.is_zero() == u32::MAX {
                    continue;
                }

                let mut image_points = [P];
                // compute isogeny
                Curve::<Fp>::velu_odd_isogeny_proj(&mut A24, &mut C24, &K, degree as usize, &mut image_points);

                // mark one step as done
                sk[i] -= 1;
            }

            done = sk.iter().sum();

            println!("done: {}", done);
        }

    
        CsidhPublicKey{ A: Curve::<Fp>::curve_from_A24_proj(&A24, &C24).A }
    }


    pub fn keygen<R: CryptoRng + RngCore>(self, 
        rng: &mut R) -> (CsidhPrivateKey<COUNT>, CsidhPublicKey<Fp>) {
        let sk = self.sample_secret_key();

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

