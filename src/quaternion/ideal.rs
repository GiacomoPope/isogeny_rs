use crate::quaternion::algebra::{IntQuat, QuatConfig, RatQuat};
use crate::quaternion::lattice::{IntLattice, RatLattice};
use crate::bigint::BigIntAlg;

#[derive(Debug)]
pub struct QuatLeftIdeal<T: BigIntAlg, P: QuatConfig<T>> {
    pub lattice: RatLattice<T, P>,
    pub norm: T,
    pub parent_order: RatLattice<T, P>,
}

impl<T: BigIntAlg, P: QuatConfig<T>> Clone for QuatLeftIdeal<T, P> {
    fn clone(&self) -> Self {
        Self {
            lattice: self.lattice.clone(),
            norm: self.norm.clone(),
            parent_order: self.parent_order.clone(),
        }
    }
}

impl<T: BigIntAlg, P: QuatConfig<T>> QuatLeftIdeal<T, P> {
    pub fn update_norm(&mut self) {
        let index = RatLattice::index(&self.lattice, &self.parent_order);
        self.norm = index.sqrt(); 
    }

    pub fn verify_norm(&self) -> u32 {
        let index = RatLattice::index(&self.lattice, &self.parent_order);
        self.norm.ct_eq(&index.sqrt())
    }

    pub fn ct_eq(&self, other: &Self) -> u32 {
        RatLattice::equal(&self.parent_order, &other.parent_order)
            & self.norm.ct_eq(&other.norm)
            & RatLattice::equal(&self.lattice, &other.lattice)
    }

    pub fn create_principal(order: &RatLattice<T, P>, x: &RatQuat<T, P>) -> Self {
        let mut lattice = RatLattice::alg_elem_mul(order, x);
        lattice.reduce_denom();

        let (norm_num, norm_denom) = x.norm();
        debug_assert!(
            norm_denom.ct_eq(&T::one()) == u32::MAX, 
            "Principal ideal generator must have integer norm"
        );

        Self {
            lattice,
            norm: norm_num,
            parent_order: order.clone(),
        }
    }

    pub fn create(order: &RatLattice<T, P>, x: &RatQuat<T, P>, n: &T) -> Self {
        debug_assert!(!x.is_zero(), "Generator x cannot be zero");

        let principal = Self::create_principal(order, x);

        let mut order_n = order.clone();
        for i in 0..4 {
            order_n.basis.generators[i] = &order_n.basis.generators[i] * n;
        }

        let mut lattice = RatLattice::add_lazy(&principal.lattice, &order_n);
        lattice.hnf();

        let mut ideal = Self {
            lattice,
            norm: T::zero(),
            parent_order: order.clone(),
        };
        ideal.update_norm();
        ideal
    }

    pub fn add(i1: &Self, i2: &Self) -> Self {
        debug_assert!(
            RatLattice::equal(&i1.parent_order, &i2.parent_order) == u32::MAX,
            "Ideals must have the same parent order"
        );

        let mut sum_lattice = RatLattice::add_lazy(&i1.lattice, &i2.lattice);
        sum_lattice.hnf();

        let mut sum_ideal = Self {
            lattice: sum_lattice,
            norm: T::zero(),
            parent_order: i1.parent_order.clone(),
        };
        sum_ideal.update_norm();
        sum_ideal
    }

    pub fn intersect(i1: &Self, i2: &Self) -> Self {
        debug_assert!(
            RatLattice::equal(&i1.parent_order, &i2.parent_order) == u32::MAX,
            "Ideals must have the same parent order"
        );

        let mut inter_lattice = RatLattice::intersect(&i1.lattice, &i2.lattice);
        inter_lattice.hnf();

        let mut inter_ideal = Self {
            lattice: inter_lattice,
            norm: T::zero(),
            parent_order: i1.parent_order.clone(),
        };
        inter_ideal.update_norm();
        inter_ideal
    }

    pub fn mul(ideal: &Self, alpha: &RatQuat<T, P>) -> Self {
        let mut prod_lattice = RatLattice::alg_elem_mul(&ideal.lattice, alpha);
        prod_lattice.hnf();

        let (alpha_norm_num, alpha_norm_denom) = alpha.norm();

        let mut prod_norm = ideal.norm.clone() * alpha_norm_num;
        debug_assert!(
            (prod_norm.clone() % alpha_norm_denom.clone()).is_zero(),
            "Norm division failure during ideal multiplication"
        );
        prod_norm = prod_norm / alpha_norm_denom;

        let prod_ideal = Self {
            lattice: prod_lattice,
            norm: prod_norm,
            parent_order: ideal.parent_order.clone(),
        };

        debug_assert!(
            prod_ideal.verify_norm() == u32::MAX, 
            "Ideal norm verification failed after multiplication"
        );

        prod_ideal
    }

    pub fn mul_ideal(i1: &Self, i2: &Self) -> Self {
        debug_assert!(
            RatLattice::equal(&i1.parent_order, &i2.parent_order) == u32::MAX,
            "Ideals must have the same parent order"
        );

        let mut prod_lattice = RatLattice::mul_lazy(&i1.lattice, &i2.lattice);
        prod_lattice.hnf();

        Self {
            lattice: prod_lattice,
            norm: i1.norm.clone() * i2.norm.clone(),
            parent_order: i1.parent_order.clone(),
        }
    }

    pub fn conjugate_without_hnf(&self, new_parent_order: &RatLattice<T, P>) -> Self {
        Self {
            lattice: self.lattice.conjugate_without_hnf(),
            norm: self.norm.clone(),
            parent_order: new_parent_order.clone(),
        }
    }

    pub fn generator(&self) -> Option<RatQuat<T, P>> {
        let mut int_norm = 1_i32;

        loop {
            if int_norm > 150 { return None; }

            for a in -int_norm..=int_norm {
                for b in (-int_norm + a.abs())..=(int_norm - a.abs()) {
                    for c in (-int_norm + a.abs() + b.abs())..=(int_norm - a.abs() - b.abs()) {
                        let d = int_norm - a.abs() - b.abs() - c.abs();

                        let a_t = T::from_i32(a);
                        let b_t = T::from_i32(b);
                        let c_t = T::from_i32(c);
                        let d_t = T::from_i32(d);

                        let gcd_val = a_t.gcd(&b_t).gcd(&c_t).gcd(&d_t);

                        if gcd_val.ct_eq(&T::one()) == u32::MAX {
                            let mut coords = [T::zero(), T::zero(), T::zero(), T::zero()];
                            let scalars = [&a_t, &b_t, &c_t, &d_t];

                            for i in 0..4 {
                                let mut sum = T::zero();
                                for j in 0..4 {
                                    sum = sum + self.lattice.basis.generators[j].coords[i].clone() * scalars[j].clone();
                                }
                                coords[i] = sum;
                            }

                            let gen_int = IntQuat::new(coords[0].clone(), coords[1].clone(), coords[2].clone(), coords[3].clone());
                            let gen_rat = RatQuat::new(gen_int, self.lattice.denom.clone());

                            let (norm_n, norm_d) = gen_rat.norm();
                            if norm_d.ct_eq(&T::one()) == u32::MAX {
                                let q = norm_n.clone() / self.norm.clone();
                                let r = norm_n.clone() % self.norm.clone();

                                if r.is_zero() && self.norm.gcd(&q).ct_eq(&T::one()) == u32::MAX {
                                    return Some(gen_rat);
                                }
                            }
                        }
                    }
                }
            }
            int_norm += 1;
        }
    }
}

// ========================================================================
// Order Operations
// ========================================================================

pub fn quat_order_discriminant<T: BigIntAlg, P: QuatConfig<T>>(order: &RatLattice<T, P>) -> T {
    let p = P::p();

    let gram = [
        [T::from_i32(2), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::from_i32(2), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::from_i32(2) * p.clone(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::from_i32(2) * p.clone()],
    ];

    let mut prod = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];

    for i in 0..4 {
        for j in 0..4 {
            let mut sum = T::zero();
            for k in 0..4 {
                sum = sum + order.basis.generators[i].coords[k].clone() * gram[k][j].clone();
            }
            prod[i][j] = sum;
        }
    }

    let mut final_prod = [
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
        [T::zero(), T::zero(), T::zero(), T::zero()],
    ];

    for i in 0..4 {
        for j in 0..4 {
            let mut sum = T::zero();
            for k in 0..4 {
                sum = sum + prod[i][k].clone() * order.basis.generators[j].coords[k].clone();
            }
            final_prod[i][j] = sum;
        }
    }

    let mut final_mat: IntLattice<T, P> = IntLattice::zero();
    for i in 0..4 {
        final_mat.generators[i] = IntQuat::new(
            final_prod[i][0].clone(), final_prod[i][1].clone(),
            final_prod[i][2].clone(), final_prod[i][3].clone()
        );
    }

    let det = final_mat.inv_with_det(None);

    let denom_sq = order.denom.clone() * order.denom.clone();
    let denom_4 = denom_sq.clone() * denom_sq;
    let denom_8 = denom_4.clone() * denom_4.clone();

    let sqr = det / denom_8;
    sqr.sqrt()
}

pub fn quat_order_is_maximal<T: BigIntAlg, P: QuatConfig<T>>(order: &RatLattice<T, P>) -> u32 {
    let disc = quat_order_discriminant(order);
    disc.ct_eq(P::p())
}
