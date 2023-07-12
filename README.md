# ec_generic

This crate is a minimal and easy-to-use elliptic curve library. This library
allows to perform the following operations over an elliptic curve finite cyclic
group:

- Point Addition: `R = P + Q`
- Point Doubling: `R = P + P = 2 * P`
- Scalar Multiplication: `R = d * P`

A generic elliptic curve is defined as `y^2 = x^3 + ax + b mod p`, and in this
particular library the constrains are:

 - `p` should be **a prime number** bigger than 3
 - `4 a^3 + 27 b^2 != 0`

The library could be use in any cryptographic algorithm that requires elliptic
curve groups, for example:

- Digital Signature Algorithm (DSA)
- Zero-Knowledge Proofs (ZKP)

# Usage

This crate is [on crates.io](https://crates.io/crates/regex) and can be
used by adding `regex` to your dependencies in your project's `Cargo.toml`.

```toml
[dependencies]
ec_generic = "0.1.12"
```

## Example: `y^2 = x^3 + 2x + 2 mod 17`

```rust
use ec_generic::{EllipticCurve, Point};
use num_bigint::BigUint;

let ec = EllipticCurve {
    a: BigUint::from(2u32),
    b: BigUint::from(2u32),
    p: BigUint::from(17u32),
};

// (6,3) + (5,1) = (10,6)
let p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
let pr = Ok(Point::Coor(BigUint::from(10u32), BigUint::from(6u32)));

let res = ec.add(&p1, &p2);
assert_eq!(res, pr);

