use core::{error::Error, fmt::Display};
use std::marker::PhantomData;

use fp2::fq::Fq as FqTrait;

use crate::elliptic::{basis::BasisX, curve::Curve, two_isogeny_chain::two_isogeny_chain};

/// Various Errors for SQIsign
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
pub struct Sqisign<Fq: FqTrait> {
    security_bits: usize,
    cofactor: u8,
    cofactor_bitsize: usize,
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

/// SQIsign signature consists of
#[derive(Clone, Copy, Debug)]
pub struct SqisignSignature<'a, Fq: FqTrait> {
    aux_curve: Curve<Fq>,
    backtracking: usize,
    two_resp_length: usize,
    // TODO: should I make these array rather than slices by adding
    // in some consts to the structure? I know that aij will have
    // length (response_length + 9) // 8 and scalar will have length
    // security_bits // 8.
    aij: [&'a [u8]; 4],
    chl_scalar: &'a [u8],
    hint_aux: u8,
    hint_chl: u8,
}

impl<Fq: FqTrait> Sqisign<Fq> {
    pub const fn new(
        security_bits: usize,
        cofactor: u8,
        cofactor_bitsize: usize,
        f: usize,
        response_length: usize,
        hash_iterations: usize,
        pk_len: usize,
        sk_len: usize,
        sig_len: usize,
    ) -> Self {
        Self {
            security_bits,
            cofactor,
            cofactor_bitsize,
            f,
            response_length,
            hash_iterations,
            pk_len,
            _sk_len: sk_len,
            sig_len,
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
    fn decode_public_key(self, buf: &[u8]) -> Result<SqisignPublicKey<Fq>, SqisignError> {
        assert!(self.pk_len == Fq::ENCODED_LENGTH + 1);

        // Ensure that the byte length matches what is expected for the parameter sets.
        if buf.len() != self.pk_len {
            return Err(SqisignError::LengthError {
                expected_size: self.pk_len,
                given_size: buf.len(),
            });
        }

        // Decode all but the last bytes for the Montgomery coefficient A
        let pk_curve_bytes = &buf[..Fq::ENCODED_LENGTH];
        let curve = Self::decode_curve(pk_curve_bytes)?;

        // The remaining byte is the hint for the torsion basis generation
        let hint = *buf.last().unwrap();

        Ok(SqisignPublicKey { curve, hint })
    }

    /// Decode a buffer of bytes into a `SqisignPublicKey<Fq>`.
    fn decode_signature<'a>(self, buf: &'a [u8]) -> Result<SqisignSignature<'a, Fq>, SqisignError> {
        let aij_n_bytes = self.security_bits >> 3;
        let chl_n_bytes = (self.response_length + 9) >> 3;

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

        // Extract out all the buffer bytes into slices.
        let mut read = Fq::ENCODED_LENGTH;

        // Extract the bytes for the auxiliary curve
        let aux_bytes = &buf[..read];
        let aux_curve = Self::decode_curve(aux_bytes)?;

        // Extract the two u8 to track backtracking and r such that the
        // response length is 2^r.
        let backtracking = buf[read] as usize;
        read += 1;
        let two_resp_length = buf[Fq::ENCODED_LENGTH + 1] as usize;
        read += 1;

        // Extract out the four scalars used for the change of basis
        let mut aij: [&[u8]; 4] = Default::default();
        for i in 0..4 {
            aij[i] = &buf[read..read + aij_n_bytes];
            read += aij_n_bytes;
        }

        // Extract out the challenge bytes used to create the chl kernel
        let chl_scalar = &buf[read..read + chl_n_bytes];
        read += chl_n_bytes;

        // Extract out the final two bytes, used for torsion basis on E_aux
        // and E_chl
        let hint_aux = buf[read];
        read += 1;
        let hint_chl = buf[read];
        assert!(read + 1 == buf.len());

        Ok(SqisignSignature {
            aux_curve,
            backtracking,
            two_resp_length,
            aij,
            chl_scalar,
            hint_aux,
            hint_chl,
        })
    }

    fn apply_change_of_basis(
        E: &Curve<Fq>,
        B: &BasisX<Fq>,
        aij: &[&[u8]],
        bitlen: usize,
    ) -> BasisX<Fq> {
        let xR = E.ladder_biscalar(B, aij[0], aij[2], bitlen, bitlen);
        let xS = E.ladder_biscalar(B, aij[1], aij[3], bitlen, bitlen);

        let diff_a = &[0];
        let diff_b = &[0];

        let xRS = E.ladder_biscalar(B, diff_a, diff_b, bitlen, bitlen);

        BasisX::from_array([xR, xS, xRS])
    }

    pub fn verify(self, msg: &[u8], sig_bytes: &[u8], pk_bytes: &[u8]) -> bool {
        // Decode the byte encoded public key and signature.
        let pk = self.decode_public_key(pk_bytes).unwrap();
        let sig = self.decode_signature(sig_bytes).unwrap();

        // Set the modified response length from the signature data;
        let e_rsp_prime = self.response_length - sig.backtracking - sig.two_resp_length;

        // Ensure that all elements of aij have the expected bit length.
        let chl_order = e_rsp_prime + sig.two_resp_length + 2;
        for scalar in sig.aij.iter() {
            let scalar_bitlen =
                (scalar.len() << 3) - (scalar.last().unwrap().leading_zeros() as usize);
            if scalar_bitlen > chl_order {
                return false;
            }
        }

        // Create a torsion basis from the supplied hint
        let pk_basis = pk.curve.torsion_basis_2e_from_hint(
            0,
            &[self.cofactor],
            self.cofactor_bitsize,
            pk.hint,
        );

        // Compute the challenge kernel 2^bt * (P + [scalar]Q)
        let mut chl_kernel =
            pk.curve
                .three_point_ladder(&pk_basis, sig.chl_scalar, self.security_bits);
        chl_kernel = pk.curve.xmul_2e(&chl_kernel, sig.backtracking);

        // TODO: we need this chain to return errors on bad input.
        let curve_chl =
            two_isogeny_chain(&pk.curve, &chl_kernel, self.f - sig.backtracking, &mut []);

        // Compute the deterministic torsion basis on E_aux and E_chl
        let mut aux_basis = sig.aux_curve.torsion_basis_2e_from_hint(
            0,
            &[self.cofactor],
            self.cofactor_bitsize,
            sig.hint_aux,
        );
        let mut chl_basis = curve_chl.torsion_basis_2e_from_hint(
            0,
            &[self.cofactor],
            self.cofactor_bitsize,
            sig.hint_chl,
        );

        aux_basis = sig
            .aux_curve
            .basis_xmul_2e(&aux_basis, self.f - e_rsp_prime - 2);
        chl_basis =
            curve_chl.basis_xmul_2e(&chl_basis, self.f - e_rsp_prime - sig.two_resp_length - 2);

        // Apply the change of basis dictated by the matrix aij contained in the signature.
        // chl_basis = Self::apply_change_of_basis(&curve_chl, &chl_basis, &sig.aij, chl_order);

        // Double the basis down to correct the order

        return false;
    }
}
