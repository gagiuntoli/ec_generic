use num_bigint::BigUint;

///
/// A struct which implements the bottom layer finite field group needed to
/// operate with the coordinates of the elliptic curve group.
///
pub struct FiniteField {}

impl FiniteField {
    ///
    /// Adds to elements in the set
    ///
    /// `a + b = a mod p`
    ///
    pub fn add(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
        assert!(a < p, "{a} >= {p}");
        assert!(b < p, "{b} >= {p}");

        let r = a + b;
        r.modpow(&BigUint::from(1u32), p)
    }

    ///
    /// Multiplies to elements in the set
    ///
    /// `a * b = a mod p`
    ///
    pub fn mult(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
        assert!(a < p, "{a} >= {p}");
        assert!(b < p, "{b} >= {p}");

        let r = a * b;
        r.modpow(&BigUint::from(1u32), p)
    }

    ///
    /// Finds the additive inverse of an element in the set:
    ///
    /// `a + (-a) = 0 mod p`
    ///
    pub fn inv_add(a: &BigUint, p: &BigUint) -> BigUint {
        assert!(a < p, "a: {a} >= p: {p}");

        if *a == BigUint::from(0u32) {
            return a.clone();
        }

        let res = p - a;

        assert!(res < *p, "res: {res} >= p: {p}");

        res
    }

    ///
    /// Subtract two elements in the set:
    ///
    /// `a - b = a + (-b) = a mod p`
    ///
    pub fn subtract(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
        assert!(a < p, "a: {a} >= {p}");
        assert!(b < p, "b: {b} >= {p}");

        let b_inv = FiniteField::inv_add(b, p);

        assert!(b_inv < p.clone(), "b_inv: {b_inv} >= p: {p}");

        let res = FiniteField::add(a, &b_inv, p);

        assert!(b_inv < p.clone(), "res: {res} >= p: {p}");

        res
    }

    ///
    /// Finds the multiplicative inverse of an element in the set if p is a
    /// prime number using Fermat's Little Theorem:
    ///
    /// `a^(-1) mod p = a^(p-2) mod p`
    ///
    /// Such that:
    /// `a * a^(-1) = 1 mod p`
    ///
    pub fn inv_mult_prime(a: &BigUint, p: &BigUint) -> BigUint {
        assert!(a < p, "{a} >= {p}");

        a.modpow(&(p - BigUint::from(2u32)), p)
    }

    ///
    /// Divides two elements in the set:
    ///
    /// `a / b = a * b^(-1) = c mod p`
    ///
    pub fn divide(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
        assert!(a < p, "{a} >= {p}");
        assert!(b < p, "{b} >= {p}");

        let d_inv = FiniteField::inv_mult_prime(b, p);
        assert!(d_inv < p.clone(), "{a} >= {p}");

        FiniteField::mult(a, &d_inv, p)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_add_1() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::add(&c, &d, &p);

        assert_eq!(r, BigUint::from(3u32));
    }

    #[test]
    fn test_add_result_0() {
        let c = BigUint::from(10u32);
        let d = BigUint::from(1u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::add(&c, &d, &p);

        assert_eq!(r, BigUint::from(0u32));
    }

    #[test]
    fn test_add_2() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(31u32);

        let r = FiniteField::add(&c, &d, &p);

        assert_eq!(r, BigUint::from(14u32));
    }

    #[test]
    fn test_mul_1() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::mult(&c, &d, &p);

        assert_eq!(r, BigUint::from(7u32));
    }

    #[test]
    fn test_mul_2() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(51u32);

        let r = FiniteField::mult(&c, &d, &p);

        assert_eq!(r, BigUint::from(40u32));
    }

    #[test]
    fn test_inv_add_1() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let r = FiniteField::inv_add(&c, &p);

        assert_eq!(r, BigUint::from(47u32));
    }

    #[test]
    fn test_inv_add_zero() {
        let c = BigUint::from(0u32);
        let p = BigUint::from(51u32);

        let r = FiniteField::inv_add(&c, &p);

        assert_eq!(r, BigUint::from(0u32));
    }

    #[test]
    #[should_panic]
    fn test_inv_add_2() {
        let c = BigUint::from(52u32);
        let p = BigUint::from(51u32);

        FiniteField::inv_add(&c, &p);
    }

    #[test]
    fn test_inv_add_identity() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let c_inv = FiniteField::inv_add(&c, &p);

        assert_eq!(c_inv, BigUint::from(47u32));
        assert_eq!(FiniteField::add(&c, &c_inv, &p), BigUint::from(0u32));
    }

    #[test]
    fn test_subtract() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        assert_eq!(FiniteField::subtract(&c, &c, &p), BigUint::from(0u32));
    }

    #[test]
    fn test_inv_mult_prime_identity() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let c_inv = FiniteField::inv_mult_prime(&c, &p);

        // 4 * 3 mod 11 = 12 mod 11 = 1
        assert_eq!(c_inv, BigUint::from(3u32));
        assert_eq!(FiniteField::mult(&c, &c_inv, &p), BigUint::from(1u32));
    }

    #[test]
    fn test_divide() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        assert_eq!(FiniteField::divide(&c, &c, &p), BigUint::from(1u32));
    }
}
