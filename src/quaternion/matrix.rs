use crate::quaternion::algebra::{IntQuat, QuatConfig};
use crate::bigint::BigIntAlg;

pub struct MatrixUtils;

impl MatrixUtils {
    pub fn det_3x3<T: BigIntAlg>(m: &[[T; 3]; 3]) -> T {
        let term1 = m[0][0].clone() * (m[1][1].clone() * m[2][2].clone() - m[1][2].clone() * m[2][1].clone());
        let term2 = m[0][1].clone() * (m[1][0].clone() * m[2][2].clone() - m[1][2].clone() * m[2][0].clone());
        let term3 = m[0][2].clone() * (m[1][0].clone() * m[2][1].clone() - m[1][1].clone() * m[2][0].clone());
        term1 - term2 + term3
    }

    /// Computes the determinant of a 4x4 matrix represented by 4 IntQuat rows.
    /// Optionally returns the adjoint matrix (scaled inverse) as 4 IntQuat rows.
    pub fn mat_4x4_inv_with_det_as_denom<T: BigIntAlg, P: QuatConfig<T>>(
        inv: Option<&mut [IntQuat<T, P>; 4]>,
        mat: &[IntQuat<T, P>; 4],
    ) -> T {
        // TODO: this is naive way using cofactors, can be largely improved
        let mut adj = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];

        for i in 0..4 {
            for j in 0..4 {
                let mut minor = [
                    [T::zero(), T::zero(), T::zero()],
                    [T::zero(), T::zero(), T::zero()],
                    [T::zero(), T::zero(), T::zero()],
                ];

                let mut mi = 0;
                for row in 0..4 {
                    if row == i { continue; }
                    let mut mj = 0;
                    for col in 0..4 {
                        if col == j { continue; }
                        minor[mi][mj] = mat[row].coords[col].clone();
                        mj += 1;
                    }
                    mi += 1;
                }

                let mut c = Self::det_3x3(&minor);
                if (i + j) % 2 != 0 {
                    c = -c;
                }
                adj[j].coords[i] = c;
            }
        }

        let mut det = T::zero();
        for k in 0..4 {
            det = det + mat[0].coords[k].clone() * adj[k].coords[0].clone();
        }

        if let Some(inv_mat) = inv {
            *inv_mat = adj;
        }

        det
    }
}
