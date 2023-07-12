use num_bigint::BigUint;

///
/// A struct which implements the bottom layer finite field group needed to
/// operate with the coordinates of the elliptic curve group.
///
pub struct FiniteField {}

#[derive(Debug, PartialEq)]
pub enum FiniteFieldError {
    InvalidArgument(String),
    InvalidResult(String),
}

impl FiniteField {
    ///
    /// Adds to elements in the set
    ///
    /// `a + b = a mod p`
    ///
    pub fn add(a: &BigUint, b: &BigUint, p: &BigUint) -> Result<BigUint, FiniteFieldError> {
        FiniteField::check_less_than(a, p)?;
        FiniteField::check_less_than(b, p)?;

        Ok((a + b).modpow(&BigUint::from(1u32), p))
    }

    ///
    /// Multiplies to elements in the set
    ///
    /// `a * b = a mod p`
    ///
    pub fn mult(a: &BigUint, b: &BigUint, p: &BigUint) -> Result<BigUint, FiniteFieldError> {
        FiniteField::check_less_than(a, p)?;
        FiniteField::check_less_than(b, p)?;

        Ok((a * b).modpow(&BigUint::from(1u32), p))
    }

    ///
    /// Finds the additive inverse of an element in the set:
    ///
    /// `a + (-a) = 0 mod p`
    ///
    pub fn inv_add(a: &BigUint, p: &BigUint) -> Result<BigUint, FiniteFieldError> {
        FiniteField::check_less_than(a, p)?;

        if *a == BigUint::from(0u32) {
            return Ok(a.clone());
        }

        Ok(p - a)
    }

    ///
    /// Subtract two elements in the set:
    ///
    /// `a - b = a + (-b) = a mod p`
    ///
    pub fn subtract(a: &BigUint, b: &BigUint, p: &BigUint) -> Result<BigUint, FiniteFieldError> {
        FiniteField::check_less_than(a, p)?;
        FiniteField::check_less_than(b, p)?;

        let b_inv = FiniteField::inv_add(b, p)?;

        FiniteField::add(a, &b_inv, p)
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
    pub fn inv_mult_prime(a: &BigUint, p: &BigUint) -> Result<BigUint, FiniteFieldError> {
        FiniteField::check_less_than(a, p)?;

        Ok(a.modpow(&(p - BigUint::from(2u32)), p))
    }

    ///
    /// Divides two elements in the set:
    ///
    /// `a / b = a * b^(-1) = a mod p`
    ///
    pub fn divide(a: &BigUint, b: &BigUint, p: &BigUint) -> Result<BigUint, FiniteFieldError> {
        FiniteField::check_less_than(a, p)?;
        FiniteField::check_less_than(b, p)?;

        let d_inv = FiniteField::inv_mult_prime(b, p)?;

        FiniteField::mult(a, &d_inv, p)
    }

    pub fn check_less_than(a: &BigUint, b: &BigUint) -> Result<(), FiniteFieldError> {
        if a >= b {
            return Err(FiniteFieldError::InvalidArgument(format!("{} >= {}", a, b)));
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_add() {
        let a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let res = FiniteField::add(&a, &b, &p).unwrap();
        assert_eq!(res, BigUint::from(3u32));

        let a = BigUint::from(10u32);
        let b = BigUint::from(1u32);
        let p = BigUint::from(11u32);

        let res = FiniteField::add(&a, &b, &p).unwrap();
        assert_eq!(res, BigUint::from(0u32));

        let a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let p = BigUint::from(31u32);

        let res = FiniteField::add(&a, &b, &p).unwrap();
        assert_eq!(res, BigUint::from(14u32));
    }

    #[test]
    fn test_mul_1() {
        let a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let res = FiniteField::mult(&a, &b, &p).unwrap();

        assert_eq!(res, BigUint::from(7u32));
    }

    #[test]
    fn test_multiply() {
        let a = BigUint::from(4u32);
        let b = BigUint::from(10u32);
        let p = BigUint::from(51u32);

        let res = FiniteField::mult(&a, &b, &p).unwrap();
        assert_eq!(res, BigUint::from(40u32));
    }

    #[test]
    fn test_inv_add() {
        let a = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let res = FiniteField::inv_add(&a, &p).unwrap();
        assert_eq!(res, BigUint::from(47u32));

        let a = BigUint::from(0u32);
        let p = BigUint::from(51u32);

        let res = FiniteField::inv_add(&a, &p).unwrap();
        assert_eq!(res, BigUint::from(0u32));

        let a = BigUint::from(52u32);
        let p = BigUint::from(51u32);

        assert_eq!(
            FiniteField::inv_add(&a, &p),
            Err(FiniteFieldError::InvalidArgument(format!("{} >= {}", a, p)))
        );

        let a = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let c_inv = FiniteField::inv_add(&a, &p);

        assert_eq!(c_inv, Ok(BigUint::from(47u32)));
        assert_eq!(
            FiniteField::add(&a, &c_inv.unwrap(), &p),
            Ok(BigUint::from(0u32))
        );
    }

    #[test]
    fn test_subtract() {
        // a - a = 0 mod p
        let a = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        assert_eq!(FiniteField::subtract(&a, &a, &p), Ok(BigUint::from(0u32)));
    }

    #[test]
    fn test_inv_mult() {
        // 4 * 3 mod 11 = 12 mod 11 = 1
        let a = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let c_inv = FiniteField::inv_mult_prime(&a, &p);

        assert_eq!(c_inv, Ok(BigUint::from(3u32)));
        assert_eq!(
            FiniteField::mult(&a, &c_inv.unwrap(), &p),
            Ok(BigUint::from(1u32))
        );
    }

    #[test]
    fn test_divide() {
        // a / a = 1 mod p
        let a = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        assert_eq!(FiniteField::divide(&a, &a, &p), Ok(BigUint::from(1u32)));
    }
}
