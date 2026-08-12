use crate::bigint::BigIntAlg;

#[derive(Clone, Debug)]
pub struct BigRat<T: BigIntAlg> {
    pub num: T,
    pub den: T,
}

impl<T: BigIntAlg> BigRat<T> {
    pub fn new(n: T, d: T) -> Self {
        let mut res = Self { num: n, den: d };
        res.reduce();
        res
    }

    pub fn zero() -> Self {
        Self::new(T::zero(), T::one())
    }

    pub fn one() -> Self {
        Self::new(T::one(), T::one())
    }

    pub fn from_integer(n: T) -> Self {
        Self::new(n, T::one())
    }

    pub fn reduce(&mut self) {
        if self.num.is_zero() {
            self.den = T::one();
            return;
        }
        let g = self.num.gcd(&self.den);
        self.num = self.num.clone() / g.clone();
        self.den = self.den.clone() / g;
        if self.den < T::zero() {
            self.num = -self.num.clone();
            self.den = -self.den.clone();
        }
    }

    pub fn round(&self) -> T {
        let two = T::from_i32(2);
        let half_den = self.den.clone() / two;
        if self.num >= T::zero() {
            (self.num.clone() + half_den) / self.den.clone()
        } else {
            (self.num.clone() - half_den) / self.den.clone()
        }
    }
}

impl<T: BigIntAlg> core::ops::Add for BigRat<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let n = self.num.clone() * rhs.den.clone() + rhs.num.clone() * self.den.clone();
        let d = self.den * rhs.den;
        Self::new(n, d)
    }
}

impl<T: BigIntAlg> core::ops::Sub for BigRat<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let n = self.num.clone() * rhs.den.clone() - rhs.num.clone() * self.den.clone();
        let d = self.den * rhs.den;
        Self::new(n, d)
    }
}

impl<T: BigIntAlg> core::ops::Mul for BigRat<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.num * rhs.num, self.den * rhs.den)
    }
}

impl<T: BigIntAlg> core::ops::Div for BigRat<T> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        Self::new(self.num * rhs.den, self.den * rhs.num)
    }
}

impl<T: BigIntAlg> PartialEq for BigRat<T> {
    fn eq(&self, other: &Self) -> bool {
        (self.num.clone() * other.den.clone()).ct_eq(&(other.num.clone() * self.den.clone()))
            == u32::MAX
    }
}

impl<T: BigIntAlg> PartialOrd for BigRat<T> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        let lhs = self.num.clone() * other.den.clone();
        let rhs = other.num.clone() * self.den.clone();
        lhs.partial_cmp(&rhs)
    }
}
