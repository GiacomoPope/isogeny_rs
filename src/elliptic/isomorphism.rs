use fp2::fq::Fq as FqTrait;

use super::{curve::Curve, point::PointX};

/// Structure used to evaluate an isomorphism between two
/// Montgomery curves such that theta(P) = (x*X + z*Z : d*Z)
#[derive(Clone, Copy, Debug)]
pub struct Isomorphism<Fq: FqTrait> {
    x: Fq,
    z: Fq,
    d: Fq,
}

impl<Fq: FqTrait> Isomorphism<Fq> {
    // TODO: explain edge cases
    pub fn new(E1: &Curve<Fq>, E2: &Curve<Fq>) -> Self {
        // lambda_x = (2*A2^3 - 9*A2) * (3 - A1^2)
        // lambda_z = (2*A1^3 - 9*A1) * (3 - A2^2)
        let mut t1 = E1.A.square();
        let mut t2 = E2.A.square();
        let mut x = Fq::THREE - t1; // x = (3 - A1^2)
        let mut z = Fq::THREE - t2; // z = (3 - A1^2)
        t1 *= E1.A;
        t2 *= E2.A;
        x *= t2.mul2() - E2.A.mul3().mul3(); // x = (2 * A2^3 - 9 * A2) * (3 - A1^2)
        z *= t1.mul2() - E1.A.mul3().mul3(); // x = (2 * A1^3 - 9 * A1) * (3 - A2^2)

        // The above is a transformation for Short-Weierstrass, so now we map
        // Mont -> SW -> SW -> Mont again to get the isometery between two
        // Montgomery curves.
        // x = 3*lambda_x
        // z = lambda_x * A1 - lambda_z * A2
        // d = 3*lambda_z
        let d = z.mul3();
        z = x * E1.A - z * E2.A;
        x.set_mul3();

        Self { x, z, d }
    }

    pub fn isomorphism_eval(&self, P: &mut PointX<Fq>) {
        P.X *= self.x;
        let tmp = self.z * P.Z;
        P.X += tmp;
        P.Z *= self.d;
    }
}
