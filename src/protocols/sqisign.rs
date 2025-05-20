use core::{error::Error, fmt::Display};
use std::marker::PhantomData;

use ::std::time::Instant;

use fp2::fq::Fq as FqTrait;
use sha3::{
    Shake256,
    digest::{ExtendableOutputReset, Update},
};

use crate::{
    elliptic::{basis::BasisX, curve::Curve},
    theta::elliptic_product::{EllipticProduct, ProductPoint},
    utilities::le_bytes::{byte_slice_difference_into, le_bytes_bit_length},
};

/// Various Errors for SQIsign, to be modified further.
#[derive(Debug)]
pub enum SqisignError {
    LengthError {
        expected_size: usize,
        given_size: usize,
    },
    InvalidFieldEncoding,
}

impl Display for SqisignError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SqisignError::LengthError {
                expected_size,
                given_size,
            } => {
                write!(
                    f,
                    "Expected size of {expected_size} is not equal to {given_size}"
                )
            }
            SqisignError::InvalidFieldEncoding => {
                f.write_str("Decoded field element is not canonically encoded")
            }
        }
    }
}

impl Error for SqisignError {}

/// Public parameters used for SQIsign (currently verification only)
#[derive(Clone, Copy, Debug)]
pub struct SqisignParameters {
    pub security_bits: usize,
    pub cofactor: u8,
    pub f: usize,
    pub response_length: usize,
    pub hash_iterations: usize,
    pub pk_len: usize,
    pub sk_len: usize,
    pub sig_len: usize,
}

/// SQIsign type which implements main methods (currently verification only).
#[derive(Clone, Copy, Debug)]
pub struct Sqisign<Fq: FqTrait> {
    security_bits: usize,
    cofactor: u8,
    f: usize,
    response_length: usize,
    hash_iterations: usize,
    pk_len: usize,
    _sk_len: usize,
    sig_len: usize,
    _phantom: PhantomData<Fq>,
}

/// SQIsign public key holds a Montgomery Curve and a hint for computing
/// a torsion basis E[2^f].
#[derive(Clone, Copy, Debug)]
pub struct SqisignPublicKey<Fq: FqTrait> {
    curve: Curve<Fq>,
    hint: u8,
}

/// SQIsign signature consists of the curve E_aux, integers tracking whether
/// any backtracking was detected during the signature, the length of the
/// response isogeny, four scalars used for a change in bases computation,
/// the scalar to compute the challenge kernel ker(phi_chl) = P + [s]Q and
/// two hints used for deriving torsion bases on E_aux and E_chl.
#[derive(Clone, Copy, Debug)]
pub struct SqisignSignature<'a, Fq: FqTrait> {
    aux_curve: Curve<Fq>,
    backtracking: usize,
    two_resp_length: usize,
    aij: [&'a [u8]; 4],
    chl_scalar: &'a [u8],
    aux_hint: u8,
    chl_hint: u8,
}

/// SQIsign type which holds parameters for a given security level and implements the
/// methods required. Currently verification only.
impl<Fq: FqTrait> Sqisign<Fq> {
    pub const fn new(params: &SqisignParameters) -> Self {
        Self {
            security_bits: params.security_bits,
            cofactor: params.cofactor,
            f: params.f,
            response_length: params.response_length,
            hash_iterations: params.hash_iterations,
            pk_len: params.pk_len,
            _sk_len: params.sk_len,
            sig_len: params.sig_len,
            _phantom: PhantomData,
        }
    }

    /// Decode a buffer of bytes into an Elliptic Curve in Montgomery form.
    fn decode_curve(buf: &[u8]) -> Result<Curve<Fq>, SqisignError> {
        // Ensure that the value was canonically encoded.
        let (A, check) = Fq::decode(buf);
        if check != u32::MAX {
            return Err(SqisignError::InvalidFieldEncoding);
        }

        Ok(Curve::new(&A))
    }

    /// Decode a buffer of bytes into a `SqisignPublicKey<Fq>`.
    fn decode_public_key(&self, buf: &[u8]) -> Result<SqisignPublicKey<Fq>, SqisignError> {
        assert!(self.pk_len == Fq::ENCODED_LENGTH + 1);

        // Ensure that the byte length matches what is expected for the parameter sets.
        if buf.len() != self.pk_len {
            return Err(SqisignError::LengthError {
                expected_size: self.pk_len,
                given_size: buf.len(),
            });
        }

        // Decode all but the last bytes for the Montgomery coefficient A
        let (pk_curve_bytes, buf) = buf.split_at(Fq::ENCODED_LENGTH);
        let curve = Self::decode_curve(pk_curve_bytes)?;

        // The remaining byte is the hint for the torsion basis generation
        let hint = buf[0];

        Ok(SqisignPublicKey { curve, hint })
    }

    /// Decode a buffer of bytes into a `SqisignPublicKey<Fq>`.
    fn decode_signature<'a>(
        &self,
        buf: &'a [u8],
    ) -> Result<SqisignSignature<'a, Fq>, SqisignError> {
        let chl_n_bytes = self.security_bits >> 3;
        let aij_n_bytes = (self.response_length + 9) >> 3;

        // Signature is of the form:
        // Fq ele || byte || byte || a00 || a01 || a10 || a11 || chl || hint_aux || hint_chl
        assert!(self.sig_len == Fq::ENCODED_LENGTH + 2 + 4 * aij_n_bytes + chl_n_bytes + 2);

        // Ensure that the byte length matches what is expected for the parameter sets.
        if buf.len() != self.sig_len {
            return Err(SqisignError::LengthError {
                expected_size: self.pk_len,
                given_size: buf.len(),
            });
        }

        // Extract the bytes for the auxiliary curve
        let (aux_bytes, buf) = buf.split_at(Fq::ENCODED_LENGTH);
        let aux_curve = Self::decode_curve(aux_bytes)?;

        // Extract the two u8 to track backtracking and r such that the
        // response length is 2^r.
        let backtracking = buf[0] as usize;
        let two_resp_length = buf[1] as usize;
        let (_, buf) = buf.split_at(2);

        // Extract out the four scalars used for the change of basis
        let mut aij: [&[u8]; 4] = Default::default();
        let (mut aij_buf, buf) = buf.split_at(4 * aij_n_bytes);
        for scalar in &mut aij {
            (*scalar, aij_buf) = aij_buf.split_at(aij_n_bytes);
        }

        // Extract out the challenge bytes used to create the chl kernel
        let (chl_scalar, buf) = buf.split_at(chl_n_bytes);

        // Extract out the final two bytes, used for torsion basis on E_aux
        // and E_chl
        let hint_aux = buf[0];
        let hint_chl = buf[1];

        Ok(SqisignSignature {
            aux_curve,
            backtracking,
            two_resp_length,
            aij,
            chl_scalar,
            aux_hint: hint_aux,
            chl_hint: hint_chl,
        })
    }

    /// Compute the challenge curve E_chl from the scalar in the signature and
    /// a deterministic basis of E_pk.
    fn compute_challenge_curve<'a>(
        &self,
        pk: &SqisignPublicKey<Fq>,
        sig: &SqisignSignature<'a, Fq>,
    ) -> (Curve<Fq>, u32) {
        // Create a torsion basis from the supplied hint
        let pk_basis = pk.curve.torsion_basis_2e_from_hint(
            0,
            &[self.cofactor],
            (8 - self.cofactor.leading_zeros()) as usize,
            pk.hint,
        );

        // Compute the challenge kernel 2^bt * (P + [scalar]Q)
        let mut chl_kernel =
            pk.curve
                .three_point_ladder(&pk_basis, sig.chl_scalar, self.security_bits);
        chl_kernel = pk.curve.xdouble_iter(&chl_kernel, sig.backtracking);

        // Compute the isogeny chain, a failure u32 is returned for the case of
        // bad input.
        pk.curve
            .two_isogeny_chain(&chl_kernel, self.f - sig.backtracking, &mut [])
    }

    /// Given the public key curve and challenge curve, compute a challenge scalar
    /// from the message. Used for generating and verifiying signatures.
    fn hash_challenge(&self, E_pk: &Curve<Fq>, E_chl: &Curve<Fq>, msg: &[u8]) -> Vec<u8>
    where
        [(); Fq::ENCODED_LENGTH]: Sized,
    {
        let mut shake_256 = Shake256::default();

        // For all but the last steps, we extract out hash_bytes from
        // Shake256.
        let hash_bytes = ((self.security_bits << 1) + 7) >> 3;
        let mut xof_bytes = vec![0; hash_bytes];

        // The first iteration hashes j(E_pk) || j(E_chl) || msg
        shake_256.update(&E_pk.j_invariant().encode());
        shake_256.update(&E_chl.j_invariant().encode());
        shake_256.update(msg);
        shake_256.finalize_xof_reset_into(&mut xof_bytes);

        // Now iterate the hash many times to reach security goal
        // We compute xof_bytes = Shake256(xof_bytes).read(hash_bytes)
        for _ in 0..(self.hash_iterations - 2) {
            shake_256.update(&xof_bytes);
            shake_256.finalize_xof_reset_into(&mut xof_bytes);
        }

        // For the last iteration we request the number of bytes required
        // for the scalar used to generate the challenge kernel.
        let scalar_hash_bytes = (self.security_bits + 7) >> 3;
        let mut scalar = vec![0; scalar_hash_bytes];
        shake_256.update(&xof_bytes);
        shake_256.finalize_xof_reset_into(&mut scalar);

        // Finally we need to reduce this value modulo 2^(f - response_length)
        let modulus_bit_length = self.f - self.response_length;
        let scalar_len = scalar.len();
        let i = modulus_bit_length >> 3;
        if i < scalar_len {
            // Partial mask of top non-zero element.
            // Note in rust 255 >> 8 != 0, instead it returns an error
            // so we have a bit of a song and dance to do when
            // modulus_bit_length & 7 is 0.
            scalar[i] &= u8::MAX
                .checked_shr(8 - (modulus_bit_length & 7) as u32)
                .unwrap_or(0);
            // All other elements are set to zero after modulus.
            for s in scalar.iter_mut().take(scalar_len).skip(i + 1) {
                *s = 0;
            }
        }

        scalar
    }

    /// Apply the change of basis for the torsion basis <P, Q> given entries
    /// of the 2x2 matrix aij such that R = [a00] P + [a10] Q and S = [a01] P + [a11] Q.
    fn apply_change_of_basis(
        E: &Curve<Fq>,
        B: &BasisX<Fq>,
        aij: &[&[u8]],
        bitlen: usize,
    ) -> BasisX<Fq> {
        // Compute R = [a00] P + [a10] Q and S = [a01] P + [a11] Q
        let R = E.ladder_biscalar(B, aij[0], aij[2], bitlen, bitlen);
        let S = E.ladder_biscalar(B, aij[1], aij[3], bitlen, bitlen);

        // Compute a00 - a01 and a10 - a11 modulo 2^bitlen
        let mut diff_a = aij[0].to_vec();
        byte_slice_difference_into(&mut diff_a, aij[1]);

        let mut diff_b = aij[2].to_vec();
        byte_slice_difference_into(&mut diff_b, aij[3]);

        // Compute R - S = [a00 - a01] P + [a10 - a11] Q
        let RS = E.ladder_biscalar(B, &diff_a, &diff_b, bitlen, bitlen);

        BasisX::from_points(&R, &S, &RS)
    }

    /// Compute the deterministic torsion bases for E_chl and E_aux and then
    /// reduce their order and apply a change of basis computation to prepare
    /// the kernel for the (2^n, 2^n)-isogeny.
    fn compute_torsion_bases<'a>(
        &self,
        E_chl: &Curve<Fq>,
        sig: &SqisignSignature<'a, Fq>,
        e_rsp_prime: usize,
        chl_order: usize,
    ) -> (BasisX<Fq>, BasisX<Fq>) {
        // Compute the deterministic torsion basis on E_aux and E_chl
        let mut aux_basis = sig.aux_curve.torsion_basis_2e_from_hint(
            0,
            &[self.cofactor],
            (8 - self.cofactor.leading_zeros()) as usize,
            sig.aux_hint,
        );
        let mut chl_basis = E_chl.torsion_basis_2e_from_hint(
            0,
            &[self.cofactor],
            (8 - self.cofactor.leading_zeros()) as usize,
            sig.chl_hint,
        );

        // Double the bases to get points of the correct even order.
        aux_basis = sig
            .aux_curve
            .basis_double_iter(&aux_basis, self.f - e_rsp_prime - 2);
        chl_basis =
            E_chl.basis_double_iter(&chl_basis, self.f - e_rsp_prime - sig.two_resp_length - 2);

        // Apply the change of basis dictated by the matrix aij contained in the signature.
        chl_basis = Self::apply_change_of_basis(E_chl, &chl_basis, &sig.aij, chl_order);

        (chl_basis, aux_basis)
    }

    /// Compute the short 2^r isogeny and push through the torsion basis of E_chl.
    /// Returns `true` if no error was detected during the 2-isogeny chain, otherwise
    /// returns `false`.
    fn compute_small_isogeny<'a>(
        E: &mut Curve<Fq>,
        B: &mut BasisX<Fq>,
        sig: &SqisignSignature<'a, Fq>,
        e_rsp_prime: usize,
    ) -> bool {
        // Compute the kernel as [2^(e_rsp_prime + 2)] P or [2^(e_rsp_prime + 2)] Q
        // depending on aij values
        let mut basis_img = B.to_array();
        let kernel = if sig.aij[0][0] & 1 == 0 && sig.aij[2][0] & 1 == 0 {
            E.xdouble_iter(&B.Q, e_rsp_prime + 2)
        } else {
            E.xdouble_iter(&B.P, e_rsp_prime + 2)
        };

        // Compute the two isogeny and push the challenge basis through
        let ok: u32;
        (*E, ok) =
            E.two_isogeny_chain_small_vartime(&kernel, sig.two_resp_length, &mut basis_img, false);
        *B = BasisX::from_slice(&basis_img);

        ok == u32::MAX
    }

    /// Return `0` if any element of aij is larger than expected, otherwise return the maximum
    /// bit-length of the matrix coefficients.
    fn check_aij_bitlen<'a>(sig: &SqisignSignature<'a, Fq>, e_rsp_prime: usize) -> usize {
        let chl_order = e_rsp_prime + sig.two_resp_length + 2;
        for scalar in sig.aij.iter() {
            let scalar_bitlen = le_bytes_bit_length(scalar);
            if scalar_bitlen > chl_order {
                return 0;
            }
        }
        chl_order
    }

    fn prepare_product_isogeny_kernel(
        E1: &Curve<Fq>,
        E2: &Curve<Fq>,
        B1: &BasisX<Fq>,
        B2: &BasisX<Fq>,
    ) -> (EllipticProduct<Fq>, ProductPoint<Fq>, ProductPoint<Fq>) {
        let E1E2 = EllipticProduct::new(E1, E2);
        // TODO: lift_basis requires an inversion, we could write a function
        // which normallises B1 and B2 simultaneously to save one inversion
        // here.
        let (P_chl, Q_chl) = E1.lift_basis(B1);
        let (P_aux, Q_aux) = E2.lift_basis(B2);
        let P1P2 = ProductPoint::new(&P_chl, &P_aux);
        let Q1Q2 = ProductPoint::new(&Q_chl, &Q_aux);

        (E1E2, P1P2, Q1Q2)
    }

    /// SQIsign verification.
    pub fn verify(&self, msg: &[u8], sig_bytes: &[u8], pk_bytes: &[u8]) -> bool
    where
        [(); Fq::ENCODED_LENGTH]: Sized,
    {
        let start = Instant::now();

        // Decode the byte encoded public key and signature.
        let pk = self.decode_public_key(pk_bytes).unwrap();
        let sig = self.decode_signature(sig_bytes).unwrap();

        let decode = start.elapsed();
        println!("Decode time: {:?}", decode);

        // Set the modified response length from the signature data;
        let e_rsp_prime = self.response_length - sig.backtracking - sig.two_resp_length;

        // Ensure that all elements of aij have the expected bit length.
        let chl_order = Self::check_aij_bitlen(&sig, e_rsp_prime);
        if chl_order == 0 {
            return false;
        }

        // Compute the challenge kernel and from this, E_chl from E_pk / <K>
        // If the kernel is found to be maleformed in the 2^n isogeny chain.
        // check = 0 and we must reject the signature.
        let (mut chl_curve, check) = self.compute_challenge_curve(&pk, &sig);
        if check == 0 {
            return false;
        }

        let curve = start.elapsed();
        println!("challenge curve time: {:?}", curve - decode);

        // Compute canonical bases on E_chl and E_aux which will be used in the (2^n, 2^n)-isogeny
        let (mut chl_basis, aux_basis) =
            self.compute_torsion_bases(&chl_curve, &sig, e_rsp_prime, chl_order);

        let basis = start.elapsed();
        println!("basis time: {:?}", basis - curve);

        // Compute the small 2-isogeny conditionally and push through the challenge basis.
        // If the kernel is maleformed, `compute_small_isogeny` returns `false` and we must
        // reject the signature.
        if sig.two_resp_length > 0
            && !Self::compute_small_isogeny(&mut chl_curve, &mut chl_basis, &sig, e_rsp_prime)
        {
            return false;
        };

        let small_isogeny = start.elapsed();
        println!("small_isogeny time: {:?}", small_isogeny - basis);

        // In very exceptional cases, no (2,2)-isogeny is needed and the signature
        // can be verified from E_chl directly.
        if e_rsp_prime == 0 {
            return sig.chl_scalar == self.hash_challenge(&pk.curve, &chl_curve, msg);
        }

        // Create the kernel for the (2, 2) isogeny given the x-only bases on E_chl and E_aux.
        let (E1E2, P1P2, Q1Q2) = Self::prepare_product_isogeny_kernel(
            &chl_curve,
            &sig.aux_curve,
            &chl_basis,
            &aux_basis,
        );

        // Compute the isogeny E1 x E2 -> E3 x E4, target codomain is always E4 currently
        // because of linear alegbra choices.
        // If check is zero then an error was detected at some point in the (2, 2) isogeny and the signature
        // is rejected.
        let (E3E4, _, check) = E1E2.elliptic_product_isogeny(&P1P2, &Q1Q2, e_rsp_prime, &[]);
        if check == 0 {
            return false;
        }
        let (_, E4) = E3E4.curves();

        let two_isogeny = start.elapsed();
        println!("two_isogeny time: {:?}", two_isogeny - small_isogeny);

        // The signature is valid if the derived bytes from the hash match the signature scalar.
        sig.chl_scalar == self.hash_challenge(&pk.curve, &E4, msg)
    }
}
