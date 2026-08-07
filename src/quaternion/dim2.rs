use crate::bigint::BigIntAlg;
use std::fmt::Debug;

// ========================================================================
// 2D Vector and Matrix Primitives
// ========================================================================

#[derive(Clone, Debug)]
pub struct Vec2<T: BigIntAlg> {
    pub coords: [T; 2],
}

impl<T: BigIntAlg> Vec2<T> {
    pub fn new(c0: T, c1: T) -> Self {
        Self { coords: [c0, c1] }
    }

    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero())
    }
}

#[derive(Clone, Debug)]
pub struct Mat2x2<T: BigIntAlg> {
    pub m: [[T; 2]; 2],
}

impl<T: BigIntAlg> Mat2x2<T> {
    pub fn zero() -> Self {
        Self {
            m: [[T::zero(), T::zero()], [T::zero(), T::zero()]],
        }
    }

    /// Evaluates the matrix-vector multiplication
    pub fn eval(&self, vec: &Vec2<T>) -> Vec2<T> {
        let c0 = self.m[0][0].clone() * vec.coords[0].clone() + self.m[0][1].clone() * vec.coords[1].clone();
        let c1 = self.m[1][0].clone() * vec.coords[0].clone() + self.m[1][1].clone() * vec.coords[1].clone();
        Vec2::new(c0, c1)
    }

    /// Computes the exact integer adjugate matrix and its determinant
    pub fn inv_with_det_as_denom(&self) -> (Self, T) {
        let det = self.m[0][0].clone() * self.m[1][1].clone() - self.m[0][1].clone() * self.m[1][0].clone();

        let inv = Self {
            m: [
                [self.m[1][1].clone(), -self.m[0][1].clone()],
                [-self.m[1][0].clone(), self.m[0][0].clone()],
            ],
        };

        (inv, det)
    }
}

// ========================================================================
// Qlapoti 2D Geometry and CVP Utilities
// ========================================================================

/// Performs rounded division towards the nearest integer (Babai's rounding step)
pub fn rounded_div<T: BigIntAlg>(a: &T, b: &T) -> T {
    let abs_b = b.abs();
    let sign_q = a.clone() * b.clone();

    let mut q = a.clone() / b.clone();
    let r = a.clone() % b.clone();

    let abs_r = r.abs();
    let twice_abs_r = abs_r.clone() + abs_r; // 2 * |r|

    if twice_abs_r > abs_b {
        if sign_q < T::zero() {
            q = q - T::one();
        } else {
            q = q + T::one();
        }
    }

    q
}

/// Computes the norm of a 2D lattice vector: coord1^2 + norm_q * coord2^2
pub fn dim2_lattice_norm<T: BigIntAlg>(v: &Vec2<T>, norm_q: &T) -> T {
    v.coords[0].clone() * v.coords[0].clone() + v.coords[1].clone() * v.coords[1].clone() * norm_q.clone()
}

/// Computes the bilinear form of two 2D vectors: v11*v21 + norm_q * v12*v22
pub fn dim2_lattice_bilinear<T: BigIntAlg>(v1: &Vec2<T>, v2: &Vec2<T>, norm_q: &T) -> T {
    v1.coords[0].clone() * v2.coords[0].clone() + v1.coords[1].clone() * v2.coords[1].clone() * norm_q.clone()
}

/// Exact solution for the shortest vector in dimension 2 (Algorithm 3.1.14 Cohen)
pub fn dim2_lattice_short_basis<T: BigIntAlg>(basis: &Mat2x2<T>, norm_q: &T) -> Mat2x2<T> {
    let mut a = Vec2::new(basis.m[0][0].clone(), basis.m[1][0].clone());
    let mut b = Vec2::new(basis.m[0][1].clone(), basis.m[1][1].clone());

    let mut norm_a = dim2_lattice_norm(&a, norm_q);
    let mut norm_b = dim2_lattice_norm(&b, norm_q);

    if norm_a < norm_b {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut norm_a, &mut norm_b);
    }

    let mut r;
    let mut norm_t; 

    loop {
        let n = dim2_lattice_bilinear(&a, &b, norm_q);
        r = rounded_div(&n, &norm_b);

        let prod1 = T::from_i32(2) * n.clone() * r.clone();
        let prod2 = r.clone() * r.clone() * norm_b.clone();
        norm_t = norm_a.clone() - prod1 + prod2;

        if norm_b > norm_t {
            norm_a = norm_b.clone();
            norm_b = norm_t.clone();

            let t = Vec2::new(
                a.coords[0].clone() - r.clone() * b.coords[0].clone(),
                a.coords[1].clone() - r.clone() * b.coords[1].clone(),
            );

            a = b;
            b = t;
        } else {
            break;
        }
    }

    if norm_t < norm_a {
        a.coords[0] = a.coords[0].clone() - r.clone() * b.coords[0].clone();
        a.coords[1] = a.coords[1].clone() - r.clone() * b.coords[1].clone();
    }

    Mat2x2 {
        m: [
            [b.coords[0].clone(), a.coords[0].clone()],
            [b.coords[1].clone(), a.coords[1].clone()],
        ],
    }
}
