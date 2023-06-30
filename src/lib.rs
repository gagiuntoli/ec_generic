/*!
This crate is a minimal and simple-to-use elliptic curve library. This
library allows to perform the following operations over an elliptic curve
finite cyclic group:

- Point Addition: `R = P + Q`
- Point Doubling: `R = P + P = 2 * P`
- Scalar Multiplication: `R = d * P`

The library could be use in any cryptographic algorithm that requires elliptic
curve groups, for example:

- Digital Signature Algorithm (DSA)
- Zero-Knowledge Proofs (ZKP)

# Usage

This crate is [on crates.io](https://crates.io/crates/regex) and can be
used by adding `regex` to your dependencies in your project's `Cargo.toml`.

```toml
[dependencies]
ec_generic = "0.1.2"
```

# Example: `y^2 = x^3 + 2x + 2`

```rust
use ec_generic::{EllipticCurve, Point};
use num_bigint::BigUint;

fn main() {
    let ec = EllipticCurve {
        a: BigUint::from(2u32),
        b: BigUint::from(2u32),
        p: BigUint::from(17u32),
    };

    // (6,3) + (5,1) = (10,6)
    let p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
    let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
    let pr = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));

    let res = ec.add(&p1, &p2);
    assert_eq!(res, pr);

    let res = ec.add(&p2, &p1);
    assert_eq!(res, pr);
}
```

# Example: `secp256k1`: `y^2 = x^3 + 7`

```rust
use ec_generic::{EllipticCurve, Point};
use num_bigint::BigUint;

fn main() {
    let p = BigUint::parse_bytes(
        b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
        16,
    )
    .expect("could not convert p");

    let n = BigUint::parse_bytes(
        b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
        16,
    )
    .expect("could not convert n");

    let gx = BigUint::parse_bytes(
        b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
        16,
    )
    .expect("could not convert gx");

    let gy = BigUint::parse_bytes(
        b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
        16,
    )
    .expect("could not convert gy");

    let ec = EllipticCurve {
        a: BigUint::from(0u32),
        b: BigUint::from(7u32),
        p,
    };

    let g = Point::Coor(gx, gy);

    // n * G = I (Identity)
    let res = ec.scalar_mul(&g, &n);

    assert_eq!(res, Point::Identity);
}
```

*/

use num_bigint::BigUint;

#[derive(PartialEq, Clone, Debug)]
pub enum Point {
    Coor(BigUint, BigUint),
    Identity,
}

#[derive(PartialEq, Clone, Debug)]
pub struct EllipticCurve {
    // y^2 = x^2 + a * x + b
    pub a: BigUint,
    pub b: BigUint,
    pub p: BigUint,
}

impl EllipticCurve {
    ///
    /// Perform a point addition: C = A + B where A and B are points which
    /// belong to the curve. Geometrically speaking, the point C is the
    /// x-reflection of the intersection of the lines that passes through A and
    /// B and intersects the curve.
    ///
    pub fn add(&self, c: &Point, d: &Point) -> Point {
        assert!(self.is_on_curve(c), "Point is not in curve");
        assert!(self.is_on_curve(d), "Point is not in curve");
        assert!(*c != *d, "Points should not be the same");

        match (c, d) {
            (Point::Identity, _) => d.clone(),
            (_, Point::Identity) => c.clone(),
            (Point::Coor(x1, y1), Point::Coor(x2, y2)) => {
                // Check that they are not additive inverses
                let y1plusy2 = FiniteField::add(&y1, &y2, &self.p);
                if x1 == x2 && y1plusy2 == BigUint::from(0u32) {
                    return Point::Identity;
                }

                // s = (y2 - y1) / (x2 - x1) mod p
                // x3 = s^2 - x1 - x2  mod p
                // y3 = s(x1 - x3) - y1 mod p
                let numerator = FiniteField::subtract(y2, y1, &self.p);
                let denominator = FiniteField::subtract(x2, x1, &self.p);
                let s = FiniteField::divide(&numerator, &denominator, &self.p);

                let (x3, y3) = self.compute_x3_y3(&x1, &y1, &x2, &s);
                Point::Coor(x3, y3)
            }
        }
    }

    ///
    /// Perform a point doubling: B = A + A = 2 * A where A is a point in the
    /// curve. Geometrically speaking, the point B is the intersection of the
    /// tangent line over A that intersects the curve.
    ///
    pub fn double(&self, c: &Point) -> Point {
        assert!(self.is_on_curve(c), "Point is not in curve");

        match c {
            Point::Identity => Point::Identity,
            Point::Coor(x1, y1) => {
                // s = (3 * x1^2 + a) / (2 * y1) mod p
                // x3 = s^2 - 2 * x1 mod p
                // y3 = s(x1 - x3) - y1 mod p
                let numerator = x1.modpow(&BigUint::from(2u32), &self.p);
                let numerator = FiniteField::mult(&BigUint::from(3u32), &numerator, &self.p);
                let numerator = FiniteField::add(&self.a, &numerator, &self.p);

                let denominator = FiniteField::mult(&BigUint::from(2u32), y1, &self.p);
                let s = FiniteField::divide(&numerator, &denominator, &self.p);

                let (x3, y3) = self.compute_x3_y3(&x1, &y1, &x1, &s);
                Point::Coor(x3, y3)
            }
        }
    }

    fn compute_x3_y3(
        &self,
        x1: &BigUint,
        y1: &BigUint,
        x2: &BigUint,
        s: &BigUint,
    ) -> (BigUint, BigUint) {
        // x3 = s^2 - x1 - x2 mod p
        // y3 = s(x1 - x3) - y1 mod p
        let s2 = s.modpow(&BigUint::from(2u32), &self.p);
        let x3 = FiniteField::subtract(&s2, x1, &self.p);
        let x3 = FiniteField::subtract(&x3, x2, &self.p);

        let y3 = FiniteField::subtract(x1, &x3, &self.p);
        let y3 = FiniteField::mult(&s, &y3, &self.p);
        let y3 = FiniteField::subtract(&y3, &y1, &self.p);

        assert!(x3 < self.p, "{} >= {}", x3, &self.p);
        assert!(y3 < self.p, "{} >= {}", y3, &self.p);

        (x3, y3)
    }

    ///
    /// Perform a scalar multiplication of a point: B = d * A where A is a point
    /// in the curve and d is a positive scalar of any value.
    ///
    /// It uses the addition/doubling algorithm - B = d * A:
    ///
    /// T = A
    /// for i in range(bits of d - 1, 0)
    ///      T = 2 * T
    ///      if bit i of d == 1
    ///          T = T + A
    ///
    pub fn scalar_mul(&self, c: &Point, d: &BigUint) -> Point {
        let mut t = c.clone();
        for i in (0..(d.bits() - 1)).rev() {
            t = self.double(&t);
            if d.bit(i) {
                t = self.add(&t, c);
            }
        }
        t
    }

    pub fn is_on_curve(&self, c: &Point) -> bool {
        match c {
            Point::Coor(x, y) => {
                // y^2 = x^3 + a * x + b
                let y2 = y.modpow(&BigUint::from(2u32), &self.p);
                let x3 = x.modpow(&BigUint::from(3u32), &self.p);
                let ax = FiniteField::mult(&self.a, x, &self.p);
                let x3plusax = FiniteField::add(&x3, &ax, &self.p);
                y2 == FiniteField::add(&x3plusax, &self.b, &self.p)
            }
            Point::Identity => true,
        }
    }
}

struct FiniteField {}

impl FiniteField {
    fn add(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        // c + d = r mod p

        assert!(c < p, "{c} >= {p}");
        assert!(d < p, "{d} >= {p}");

        let r = c + d;
        r.modpow(&BigUint::from(1u32), p)
    }

    fn mult(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        // c * d = r mod p

        assert!(c < p, "{c} >= {p}");
        assert!(d < p, "{d} >= {p}");

        let r = c * d;
        r.modpow(&BigUint::from(1u32), p)
    }

    fn inv_addition(c: &BigUint, p: &BigUint) -> BigUint {
        // -c mod p

        assert!(c < p, "{c} >= {p}");

        p - c
    }

    fn subtract(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        // c - d mod p
        assert!(c < p, "{c} >= {p}");
        assert!(d < p, "{d} >= {p}");

        let d_inv = FiniteField::inv_addition(d, p);
        assert!(d_inv < p.clone(), "{d_inv} >= {p}");

        FiniteField::add(c, &d_inv, p)
    }

    // TODO: this function uses Fermat's Little Theorem and thus it is only
    // valid for a p prime only for p prime
    fn inv_multiplication(c: &BigUint, p: &BigUint) -> BigUint {
        // c^(-1) mod p = c^(p-2) mod p

        assert!(c < p, "{c} >= {p}");

        c.modpow(&(p - BigUint::from(2u32)), p)
    }

    fn divide(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        assert!(c < p, "{c} >= {p}");
        assert!(d < p, "{d} >= {p}");

        let d_inv = FiniteField::inv_multiplication(d, p);
        assert!(d_inv < p.clone(), "{c} >= {p}");

        FiniteField::mult(c, &d_inv, p)
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
    fn test_inv_addition_1() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let r = FiniteField::inv_addition(&c, &p);

        assert_eq!(r, BigUint::from(47u32));
    }

    #[test]
    #[should_panic]
    fn test_inv_addition_2() {
        let c = BigUint::from(52u32);
        let p = BigUint::from(51u32);

        FiniteField::inv_addition(&c, &p);
    }

    #[test]
    fn test_inv_addition_identity() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let c_inv = FiniteField::inv_addition(&c, &p);

        assert_eq!(c_inv, BigUint::from(47u32));
        assert_eq!(FiniteField::add(&c, &c_inv, &p), BigUint::from(0u32));
    }

    #[test]
    fn test_substract() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        assert_eq!(FiniteField::subtract(&c, &c, &p), BigUint::from(0u32));
    }

    #[test]
    fn test_inv_multiplication_identity() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let c_inv = FiniteField::inv_multiplication(&c, &p);

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

    #[test]
    fn test_ec_point_addition() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        // (6,3) + (5,1) = (10,6)
        let p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));

        let res = ec.add(&p1, &p2);
        assert_eq!(res, pr);

        let res = ec.add(&p2, &p1);
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_point_addition_identity() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        // (6,3) + (5,1) = (10,6)
        let p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let p2 = Point::Identity;
        let pr = p1.clone();

        let res = ec.add(&p1, &p2);
        assert_eq!(res, pr);

        let res = ec.add(&p2, &p1);
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_point_addition_reflected_in_x() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        // (5,16) + (5,1) = Point::Identity
        let p1 = Point::Coor(BigUint::from(5u32), BigUint::from(16u32));
        let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Identity;

        let res = ec.add(&p1, &p2);
        assert_eq!(res, pr);

        let res = ec.add(&p2, &p1);
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_point_doubling() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        // (5,1) + (5,1) = 2 (5, 1) = (6,3)
        let p1 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));

        let res = ec.double(&p1);
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_point_doubling_identity() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        // I + I = 2 I = I
        let p1 = Point::Identity;
        let pr = Point::Identity;

        let res = ec.double(&p1);
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_scalar_multiplication() {
        // y^2 = x^3 + 2x + 2 mod 17   |G| = 19  19 * A = I
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        let c = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));

        // 2 (5, 1) = (6,3)
        let pr = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let res = ec.scalar_mul(&c, &BigUint::from(2u32));
        assert_eq!(res, pr);

        // 10 (5, 1) = (7,11)
        let pr = Point::Coor(BigUint::from(7u32), BigUint::from(11u32));
        let res = ec.scalar_mul(&c, &BigUint::from(10u32));
        assert_eq!(res, pr);

        // 16 (5, 1) = (10,11)
        let pr = Point::Coor(BigUint::from(10u32), BigUint::from(11u32));
        let res = ec.scalar_mul(&c, &BigUint::from(16u32));
        assert_eq!(res, pr);

        // 17 (5, 1) = (6,14)
        let pr = Point::Coor(BigUint::from(6u32), BigUint::from(14u32));
        let res = ec.scalar_mul(&c, &BigUint::from(17u32));
        assert_eq!(res, pr);

        // 18 (5, 1) = (5,16)
        let pr = Point::Coor(BigUint::from(5u32), BigUint::from(16u32));
        let res = ec.scalar_mul(&c, &BigUint::from(18u32));
        assert_eq!(res, pr);

        // 19 (5, 1) = I
        let pr = Point::Identity;
        let res = ec.scalar_mul(&c, &BigUint::from(19u32));
        assert_eq!(res, pr);
    }

    #[test]
    fn test_ec_secp256k1() {
        /*
          y^2 = x^3 + 7

          p = FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F
          n = FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141
        G = (
            x = 79BE667E F9DCBBAC 55A06295 CE870B07 029BFCDB 2DCE28D9 59F2815B 16F81798,
            y = 483ADA77 26A3C465 5DA4FBFC 0E1108A8 FD17B448 A6855419 9C47D08F FB10D4B8
        )
        a = 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
        b = 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000007
        */

        // n * G = I
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .expect("could not convert p");

        let n = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .expect("could not convert n");

        let gx = BigUint::parse_bytes(
            b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .expect("could not convert gx");

        let gy = BigUint::parse_bytes(
            b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .expect("could not convert gy");

        let ec = EllipticCurve {
            a: BigUint::from(0u32),
            b: BigUint::from(7u32),
            p,
        };

        let g = Point::Coor(gx, gy);

        let res = ec.scalar_mul(&g, &n); // n * G

        assert_eq!(res, Point::Identity);
    }

    #[test]
    fn test_bits() {
        let a = BigUint::from(2u32);
        assert!(!a.bit(0));
        assert!(a.bit(1));
    }
}
