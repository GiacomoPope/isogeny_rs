use crate::bigint::BigIntAlg;

pub trait HnfModExt: BigIntAlg {
    /// Strictly positive modulo in [0, m)
    fn positive_mod(&self, m: &Self) -> Self {
        let mut r = self.clone() % m.clone();
        if r < Self::zero() {
            r = r + m.clone();
        }
        r
    }

    /// Centered modulo in (-m/2, m/2]
    fn centered_mod(&self, m: &Self) -> Self {
        let r = self.positive_mod(m);
        let two = Self::from_i32(2);
        let d = m.clone() / two;
        if r > d {
            r - m.clone()
        } else {
            r
        }
    }
}

impl<T: BigIntAlg> HnfModExt for T {}
