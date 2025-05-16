#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

#[cfg(test)]
mod test_torsion_basis {
    use fp2::fq::Fq;
    use isogeny::elliptic::basis::BasisX;
    use isogeny::elliptic::curve::Curve;
    use isogeny::elliptic::point::PointX;
    use isogeny::fields::sqisign::SqiSignI as Fp2;

    // Characteristic: p = 5*2^248 - 1
    static f: usize = 248;
    static e: usize = 124;

    // Montgomery coefficient with A a NQR in GF(p^2)
    static A_NQR_RE_BYTES: [u8; 32] = [
        203, 54, 167, 139, 157, 158, 26, 48, 238, 132, 6, 103, 183, 231, 90, 126, 110, 104, 120,
        168, 87, 156, 88, 211, 237, 158, 241, 236, 134, 18, 233, 4,
    ];
    static A_NQR_IM_BYTES: [u8; 32] = [
        190, 27, 176, 40, 178, 186, 232, 103, 26, 247, 113, 221, 233, 255, 3, 61, 39, 202, 84, 231,
        25, 178, 88, 174, 237, 194, 209, 208, 111, 84, 74, 3,
    ];
    static A_NQR: Fp2 = Fp2::const_decode_no_check(&A_NQR_RE_BYTES, &A_NQR_IM_BYTES);

    // Montgomery coefficient with A a QR in GF(p^2)
    static A_QR_RE_BYTES: [u8; 32] = [
        126, 182, 40, 112, 55, 198, 126, 217, 66, 174, 71, 99, 110, 224, 253, 208, 144, 241, 16,
        117, 107, 102, 14, 199, 220, 232, 178, 111, 39, 214, 26, 1,
    ];
    static A_QR_IM_BYTES: [u8; 32] = [
        9, 153, 229, 111, 54, 112, 39, 172, 235, 167, 145, 242, 163, 251, 111, 10, 244, 2, 145, 40,
        144, 18, 127, 234, 236, 94, 76, 59, 152, 60, 191, 4,
    ];
    static A_QR: Fp2 = Fp2::const_decode_no_check(&A_QR_RE_BYTES, &A_QR_IM_BYTES);

    fn point_has_order_n(E: &Curve<Fp2>, P: &PointX<Fp2>, n: usize) -> bool {
        let P2 = E.xdouble_iter(P, n - 1);
        if P2.is_zero() == u32::MAX {
            return false;
        }
        let O = E.xdouble_iter(&P2, 1);
        if O.is_zero() != u32::MAX {
            return false;
        }
        return true;
    }

    fn basis_has_order_n(E: &Curve<Fp2>, basis: &BasisX<Fp2>, n: usize) -> bool {
        let ok1 = point_has_order_n(&E, &basis.P, n);
        let ok2 = point_has_order_n(&E, &basis.Q, n);
        let ok3 = point_has_order_n(&E, &basis.PQ, n);
        return ok1 && ok2 && ok3;
    }

    fn basis_are_equal(b1: &BasisX<Fp2>, b2: &BasisX<Fp2>) -> bool {
        let ok1 = b1.P.equals(&b2.P);
        let ok2 = b1.Q.equals(&b2.Q);
        let ok3 = b1.PQ.equals(&b2.PQ);

        (ok1 & ok2 & ok3) == u32::MAX
    }

    #[test]
    fn test_full_order_qr() {
        assert!(A_QR.is_square() == u32::MAX);
        let E = Curve::new(&A_QR);
        let (basis, _) = E.torsion_basis_2e_with_hint(0, &[5u8], 3);
        assert!(basis_has_order_n(&E, &basis, f));
    }

    #[test]
    fn test_full_order_nqr() {
        assert!(A_NQR.is_square() == 0);
        let E = Curve::new(&A_NQR);
        let (basis, _) = E.torsion_basis_2e_with_hint(0, &[5u8], 3);
        assert!(basis_has_order_n(&E, &basis, f));
    }

    #[test]
    fn test_reduced_order_qr() {
        assert!(A_QR.is_square() == u32::MAX);
        let E = Curve::new(&A_QR);
        let (basis, _) = E.torsion_basis_2e_with_hint(f - e, &[5u8], 3);
        assert!(basis_has_order_n(&E, &basis, e));
    }

    #[test]
    fn test_reduced_order_nqr() {
        assert!(A_NQR.is_square() == 0);
        let E = Curve::new(&A_NQR);
        let (basis, _) = E.torsion_basis_2e_with_hint(f - e, &[5u8], 3);
        assert!(basis_has_order_n(&E, &basis, e));
    }

    #[test]
    fn test_hint_fallback_qr() {
        assert!(A_QR.is_square() == u32::MAX);
        let E = Curve::new(&A_QR);
        let basis = E.torsion_basis_2e_from_hint(0, &[5u8], 3, 1);
        assert!(basis_has_order_n(&E, &basis, f));
    }

    #[test]
    fn test_hint_fallback_nqr() {
        assert!(A_NQR.is_square() == 0);
        let E = Curve::new(&A_NQR);
        let basis = E.torsion_basis_2e_from_hint(0, &[5u8], 3, 0);
        assert!(basis_has_order_n(&E, &basis, f));
    }

    #[test]
    fn test_hint_is_canonical_full_qr() {
        assert!(A_QR.is_square() == u32::MAX);
        let E = Curve::new(&A_QR);
        let (b1, hint) = E.torsion_basis_2e_with_hint(0, &[5u8], 3);
        assert!(hint & 1 == 1);

        let b2 = E.torsion_basis_2e_from_hint(0, &[5u8], 3, hint);
        assert!(basis_are_equal(&b1, &b2));
    }

    #[test]
    fn test_hint_is_canonical_full_nqr() {
        assert!(A_NQR.is_square() == 0);
        let E = Curve::new(&A_NQR);
        let (b1, hint) = E.torsion_basis_2e_with_hint(0, &[5u8], 3);
        assert!(hint & 1 == 0);

        let b2 = E.torsion_basis_2e_from_hint(0, &[5u8], 3, hint);
        assert!(basis_are_equal(&b1, &b2));
    }

    #[test]
    fn test_hint_is_canonical_qr() {
        assert!(A_QR.is_square() == u32::MAX);
        let E = Curve::new(&A_QR);
        let (b1, hint) = E.torsion_basis_2e_with_hint(f - e, &[5u8], 3);
        assert!(hint & 1 == 1);

        let b2 = E.torsion_basis_2e_from_hint(f - e, &[5u8], 3, hint);
        assert!(basis_are_equal(&b1, &b2));
    }

    #[test]
    fn test_hint_is_canonical_nqr() {
        assert!(A_NQR.is_square() == 0);
        let E = Curve::new(&A_NQR);
        let (b1, hint) = E.torsion_basis_2e_with_hint(f - e, &[5u8], 3);
        assert!(hint & 1 == 0);

        let b2 = E.torsion_basis_2e_from_hint(f - e, &[5u8], 3, hint);
        assert!(basis_are_equal(&b1, &b2));
    }
}
