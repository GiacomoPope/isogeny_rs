use core::ops::{Add, Mul};
use crate::quaternion::algebra::{QuatConfig, IntQuat, RatQuat};
use crate::quaternion::hnf::quat_hnf_mod_core;
use crate::bigint::BigIntAlg;

#[derive(Debug)]
pub struct IntLattice<T: BigIntAlg, P: QuatConfig<T>> {
    pub generators: [IntQuat<T, P>; 4],
}

impl<T: BigIntAlg, P: QuatConfig<T>> Clone for IntLattice<T, P> {
    fn clone(&self) -> Self {
        Self {
            generators: self.generators.clone(),
        }
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> IntLattice<T, P> {
    pub fn zero() -> Self {
        Self {
            generators: [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()],
        }
    }

    pub fn ct_eq(&self, other: &Self) -> u32 {
        self.generators[0].ct_eq(&other.generators[0])
            & self.generators[1].ct_eq(&other.generators[1])
            & self.generators[2].ct_eq(&other.generators[2])
            & self.generators[3].ct_eq(&other.generators[3])
    }

    pub fn scale(&mut self, scalar: &T) {
        for i in 0..4 {
            self.generators[i] = &self.generators[i] * scalar;
        }
    }

    fn det_3x3(m: &[[T; 3]; 3]) -> T {
        let term1 = m[0][0].clone() * (m[1][1].clone() * m[2][2].clone() - m[1][2].clone() * m[2][1].clone());
        let term2 = m[0][1].clone() * (m[1][0].clone() * m[2][2].clone() - m[1][2].clone() * m[2][0].clone());
        let term3 = m[0][2].clone() * (m[1][0].clone() * m[2][1].clone() - m[1][1].clone() * m[2][0].clone());
        term1 - term2 + term3
    }

    // Method from https://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf 3rd of May 2023, 16h15 CEST
    pub fn inv_with_det(&self, inv: Option<&mut IntLattice<T, P>>) -> T {
        let mut adj = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];
        let mat = &self.generators;

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
            inv_mat.generators = adj;
        }

        det
    }

    pub fn gram(&self) -> [[T; 4]; 4] {
        let mut g_mat = [
            [T::zero(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::zero()],
            [T::zero(), T::zero(), T::zero(), T::zero()],
        ];
        let two = T::from_i32(2);
        let p = P::p();

        for i in 0..4 {
            for j in 0..=i {
                let mut sum = self.generators[i].coords[0].clone() * self.generators[j].coords[0].clone();
                sum = sum + self.generators[i].coords[1].clone() * self.generators[j].coords[1].clone();
                let mut p_sum = self.generators[i].coords[2].clone() * self.generators[j].coords[2].clone();
                p_sum = p_sum + self.generators[i].coords[3].clone() * self.generators[j].coords[3].clone();
                sum = sum + p.clone() * p_sum;
                g_mat[i][j] = sum * two.clone();
            }
        }

        for i in 0..4 {
            for j in (i + 1)..4 {
                g_mat[i][j] = g_mat[j][i].clone();
            }
        }
        g_mat
    }

    pub fn hnf(&mut self) {
        let mod_val = self.inv_with_det(None).abs();
        self.generators = quat_hnf_mod_core(&mut self.generators, &mod_val);
    }

    pub fn add_lazy(lat1: &Self, lat2: &Self) -> Self {
        let mut generators = Vec::with_capacity(8);
        for i in 0..4 { generators.push(lat1.generators[i].clone()); }
        for i in 0..4 { generators.push(lat2.generators[i].clone()); }

        let det1 = lat1.inv_with_det(None);
        let det2 = lat2.inv_with_det(None);

        Self {
            generators: quat_hnf_mod_core(&mut generators, &det1.gcd(&det2)),
        }
    }

    pub fn mul_lazy(lat1: &Self, lat2: &Self) -> Self {
        let mut generators = Vec::with_capacity(16);
        let mut detmat = IntLattice::zero();

        for k in 0..4 {
            for i in 0..4 {
                let elem_res = &lat1.generators[k] * &lat2.generators[i];
                if k == 0 { detmat.generators[i] = elem_res.clone(); }
                generators.push(elem_res);
            }
        }

        let mod_val = detmat.inv_with_det(None).abs();
        Self {
            generators: quat_hnf_mod_core(&mut generators, &mod_val),
        }
    }

    pub fn alg_elem_mul(lat: &Self, elem: &IntQuat<T, P>) -> Self {
        let mut generators = Vec::with_capacity(4);
        for i in 0..4 { generators.push(&lat.generators[i] * elem); }

        let mut mat = IntLattice::zero();
        for i in 0..4 { mat.generators[i] = generators[i].clone(); }

        let mod_val = mat.inv_with_det(None).abs();
        Self {
            generators: quat_hnf_mod_core(&mut generators, &mod_val),
        }
    }
}

#[derive(Debug)]
pub struct RatLattice<T: BigIntAlg, P: QuatConfig<T>> {
    pub basis: IntLattice<T, P>,
    pub denom: T,
}

impl<T: BigIntAlg, P: QuatConfig<T>> Clone for RatLattice<T, P> {
    fn clone(&self) -> Self {
        Self {
            basis: self.basis.clone(),
            denom: self.denom.clone(),
        }
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> RatLattice<T, P> {
    pub fn zero() -> Self {
        Self {
            basis: IntLattice::zero(),
            denom: T::one(),
        }
    }

    pub fn reduce_denom(&mut self) {
        let mut g = self.denom.clone();
        for i in 0..4 {
            g = g.gcd(&self.basis.generators[i].coords_gcd());
        }

        let mut sign = T::one();
        if self.denom < T::zero() {
            sign = T::from_i32(-1);
        }

        let divisor = g * sign;
        for i in 0..4 {
            self.basis.generators[i] = &self.basis.generators[i] / &divisor;
        }
        self.denom = self.denom.clone() / divisor;
    }

    pub fn hnf(&mut self) {
        self.basis.hnf();
        self.reduce_denom();
    }

    pub fn ct_eq(&self, other: &Self) -> u32 {
        self.denom.ct_eq(&other.denom) & self.basis.ct_eq(&other.basis)
    }

    pub fn equal(lat1: &Self, lat2: &Self) -> u32 {
        let mut a = lat1.clone();
        let mut b = lat2.clone();
        a.reduce_denom();
        b.reduce_denom();
        a.denom = a.denom.abs();
        b.denom = b.denom.abs();
        a.hnf();
        b.hnf();

        a.ct_eq(&b)
    }

    pub fn inclusion(sublat: &Self, overlat: &Self) -> u32 {
        let sum = Self::add_lazy(overlat, sublat);
        Self::equal(&sum, overlat)
    }

    pub fn conjugate_without_hnf(&self) -> Self {
        Self {
            basis: IntLattice {
                generators: [
                    self.basis.generators[0].conj(),
                    self.basis.generators[1].conj(),
                    self.basis.generators[2].conj(),
                    self.basis.generators[3].conj(),
                ],
            },
            denom: self.denom.clone(),
        }
    }

    // Method described in https://cseweb.ucsd.edu/classes/sp14/cse206A-a/lec4.pdf consulted on 19 of May 2023, 12h40 CEST
    pub fn dual_without_hnf(&self) -> Self {
        let mut inv = IntLattice::zero();
        let det = self.basis.inv_with_det(Some(&mut inv));

        let mut dual_gens = [IntQuat::zero(), IntQuat::zero(), IntQuat::zero(), IntQuat::zero()];

        for i in 0..4 {
            for j in 0..4 {
                dual_gens[i].coords[j] = inv.generators[j].coords[i].clone() * self.denom.clone();
            }
        }

        Self {
            basis: IntLattice { generators: dual_gens },
            denom: det,
        }
    }

    pub fn add_lazy(lat1: &Self, lat2: &Self) -> Self {
        let mut scaled1 = lat1.basis.clone();
        scaled1.scale(&lat2.denom);

        let mut scaled2 = lat2.basis.clone();
        scaled2.scale(&lat1.denom);

        let mut res = Self {
            basis: IntLattice::add_lazy(&scaled1, &scaled2),
            denom: lat1.denom.clone() * lat2.denom.clone(),
        };
        res.reduce_denom();
        res
    }

    // method described in https://cseweb.ucsd.edu/classes/sp14/cse206A-a/lec4.pdf consulted on 19 of May 2023, 12h40 CEST
    pub fn intersect(lat1: &Self, lat2: &Self) -> Self {
        let dual1 = lat1.dual_without_hnf();
        let dual2 = lat2.dual_without_hnf();
        let dual_res = Self::add_lazy(&dual1, &dual2);
        let mut res = dual_res.dual_without_hnf();
        res.hnf(); // from SQISign repo: « could be removed if we do not expect HNF any more »
        res
    }

    pub fn alg_elem_mul(lat: &Self, elem: &RatQuat<T, P>) -> Self {
        let mut res = Self {
            basis: IntLattice::alg_elem_mul(&lat.basis, &elem.num),
            denom: lat.denom.clone() * elem.denom.clone(),
        };
        res.reduce_denom();
        res
    }

    pub fn mul_lazy(lat1: &Self, lat2: &Self) -> Self {
        let mut res = Self {
            basis: IntLattice::mul_lazy(&lat1.basis, &lat2.basis),
            denom: lat1.denom.clone() * lat2.denom.clone(),
        };
        res.reduce_denom();
        res
    }

    pub fn contains(&self, x: &RatQuat<T, P>) -> Option<[T; 4]> {
        let mut inv = IntLattice::zero();
        let det = self.basis.inv_with_det(Some(&mut inv));

        let mut work_coord = [T::zero(), T::zero(), T::zero(), T::zero()];
        for i in 0..4 {
            let mut sum = T::zero();
            for j in 0..4 {
                sum = sum + inv.generators[j].coords[i].clone() * x.num.coords[j].clone();
            }
            work_coord[i] = sum * self.denom.clone();
        }

        let prod = x.denom.clone() * det;
        let mut divisible = true;

        for i in 0..4 {
            if !(work_coord[i].clone() % prod.clone()).is_zero() {
                divisible = false;
            }
            work_coord[i] = work_coord[i].clone() / prod.clone();
        }

        if divisible {
            Some(work_coord)
        } else {
            None
        }
    }

    pub fn index(sublat: &Self, overlat: &Self) -> T {
        let det_sub = sublat.basis.inv_with_det(None);
        let mut tmp_over = overlat.denom.clone() * overlat.denom.clone();
        tmp_over = tmp_over.clone() * tmp_over.clone();
        let num = det_sub * tmp_over;

        let det_over = overlat.basis.inv_with_det(None);
        let mut tmp_sub = sublat.denom.clone() * sublat.denom.clone();
        tmp_sub = tmp_sub.clone() * tmp_sub.clone();
        let den = det_over * tmp_sub;

        (num / den).abs()
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> Add<RatLattice<T, P>> for RatLattice<T, P> {
    type Output = RatLattice<T, P>;
    fn add(self, rhs: RatLattice<T, P>) -> Self::Output {
        let mut res = RatLattice::add_lazy(&self, &rhs);
        res.hnf();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Add<&'b RatLattice<T, P>> for &'a RatLattice<T, P> {
    type Output = RatLattice<T, P>;
    fn add(self, rhs: &'b RatLattice<T, P>) -> Self::Output {
        let mut res = RatLattice::add_lazy(self, rhs);
        res.hnf();
        res
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> Mul<RatLattice<T, P>> for RatLattice<T, P> {
    type Output = RatLattice<T, P>;
    fn mul(self, rhs: RatLattice<T, P>) -> Self::Output {
        let mut res = RatLattice::mul_lazy(&self, &rhs);
        res.hnf();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Mul<&'b RatLattice<T, P>> for &'a RatLattice<T, P> {
    type Output = RatLattice<T, P>;
    fn mul(self, rhs: &'b RatLattice<T, P>) -> Self::Output {
        let mut res = RatLattice::mul_lazy(self, rhs);
        res.hnf();
        res
    }
}

impl<'a, 'b, T: BigIntAlg, P: QuatConfig<T>> Mul<&'b RatQuat<T, P>> for &'a RatLattice<T, P> {
    type Output = RatLattice<T, P>;
    fn mul(self, rhs: &'b RatQuat<T, P>) -> Self::Output {
        RatLattice::alg_elem_mul(self, rhs)
    }
}
