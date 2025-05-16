# A Rust Isogeny Library

A Rust library for isogeny-based cryptography

## :construction: Everything will Change Always :construction:

Currently this code is "pre-alpha" in that in development of new features, the code is being constantly refactored. Don't expect any of the current code to maintain API / form.

## Motivation

Over the past few years, I've written some Rust for isogeny-based research papers which have now become spread over several GitHub repositories. The aim with this project is to collect all this work into one library with a consistent API.

The hope is that after enough work, this library makes implementing new protocols in Rust more easily. This is helped thanks to the finite field macro: [`fp2`](https://github.com/GiacomoPope/fp2) which allows the easy creation of any field $GF(p^2)$ with modulus $x^2 + 1$ requiring only `p` encoded as little endian `u64` for creation.

## Protocols

This library currently contains:

- SQIsign verification following the [SQIsign spec](https://sqisign.org)
- A toy implementation of SIDH to demonstrate 2-isogenies and 3-isogenies

## Associated Work

This repository has started as a collection and refactoring of some isogeny-based cryptography research papers.

- "An Algorithmic Approach to (2, 2)-isogenies in the Theta Model and Applications to Isogeny-based Cryptography" by Pierrick Dartois, Luciano Maino, Giacomo Pope, and Damien Robert.
  - https://eprint.iacr.org/2023/1747
  - https://github.com/ThetaIsogenies/two-isogenies
- "Simpler and faster pairings from the Montgomery Ladder" by Giacomo Pope, Krijn Reijnders, Damien Robert, Alessandro Sferlazza and Benjamin Smith.
  - https://eprint.iacr.org/2025/672
  - https://github.com/GiacomoPope/cubical-pairings/

## SIDH Example

As a small example, SIDH key exchange is relatively compact and easy to read

```rs
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
    // We ignore the check during keygen as the parameters are trusted.
    let mut three_torsion_img = self.three_torsion.to_array();
    let (codomain, _) = E.two_isogeny_chain(&kernel, self.ea, &mut three_torsion_img);

    // Package the data above into public and private keys
    let public_key = SidhAlicePublicKey::new(&codomain, &three_torsion_img);
    let secret_key = SidhAlicePrivateKey::new(self.ea, scalar);
    (public_key, secret_key)
}
```

and made to emulate the maybe more familiar SageMath API:

```py
p = 2^216 * 3^137 - 1
Fp2.<i> = GF(p^2, modulus=x^2+1)
E = EllipticCurve(Fp2, [0, 6, 0, 1, 0])
P, Q = E.gens()
P2, Q2 = 3^137 * P, 3^137 * Q
P3, Q3 = 2^216 * P, 2^216 * Q

def alice_keygen(E, P2, Q2, P3, Q3):
    s = randint(0, 2**224)
    k = P2 + s*Q2
    phi = E.isogeny(k, algorithm="factored")
    return phi.codomain(), phi(P3), phi(Q3)
```

with the benefit of being signficiantly faster (run `cargo bench`):

```
Benchmarking Alice Keygen for SIKE434 Parameters
                       time:   [4.6825 ms 4.6870 ms 4.6913 ms]

sage: %timeit alice_keygen(E, P2, Q2, P3, Q3)
1.17 s ± 6.48 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

### Tests

Tests can be run:

```
cargo test
```

[//]: # (badges)

[build-image]: https://github.com/GiacomoPope/isogeny_rs/workflows/Rust/badge.svg
[build-link]: https://github.com/GiacomoPope/isogeny_rs/actions?query=workflow%3ARust


## Collaboration

I am very interested in collaboration to improve both the performance and scope of this project. Additionally, I am a mathematican first and Rust person second, so if any Rust experts have opinions / advice of making this project more idomatic to a Rust developer, please let me know.
