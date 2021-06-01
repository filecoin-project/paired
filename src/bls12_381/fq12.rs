use super::fq::FROBENIUS_COEFF_FQ12_C1;
use super::fq::{Fq, FqRepr};
use super::fq2::Fq2;
use super::fq6::Fq6;
use fff::{Field, PrimeField, PrimeFieldRepr};
use rand_core::RngCore;
use std::fmt;

/// An element of Fq12, represented by c0 + c1 * w.
#[derive(Copy, Clone, Debug, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Fq12 {
    pub c0: Fq6,
    pub c1: Fq6,
}

impl ::std::fmt::Display for Fq12 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
        write!(f, "Fq12({} + {} * w)", self.c0, self.c1)
    }
}

impl Fq12 {
    pub fn conjugate(&mut self) {
        self.c1.negate();
    }

    pub fn mul_by_014(&mut self, c0: &Fq2, c1: &Fq2, c4: &Fq2) {
        let mut aa = self.c0;
        aa.mul_by_01(c0, c1);
        let mut bb = self.c1;
        bb.mul_by_1(c4);
        let mut o = *c1;
        o.add_assign(c4);
        self.c1.add_assign(&self.c0);
        self.c1.mul_by_01(c0, &o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = bb;
        self.c0.mul_by_nonresidue();
        self.c0.add_assign(&aa);
    }

    /// Compress this point. Returns `None` if the element is not in the cyclomtomic subgroup.
    pub fn compress(&self) -> Option<Fq12Compressed> {
        if !self.is_cyc() {
            return None;
        }

        // Use torus-based compression from Section 4.1 in
        // "On Compressible Pairings and Their Computation" by Naehrig et al.
        let mut b = self.c0;

        b.c0.add_assign(&Fq2::one());
        b.mul_assign(&self.c1.inverse().unwrap());

        Some(Fq12Compressed(b))
    }

    fn is_cyc(&self) -> bool {
        // Check if a^(p^4 - p^2 + 1) == 1.
        let mut t0 = *self;
        t0.frobenius_map(4);
        t0.mul_assign(self);
        let mut t1 = *self;
        t1.frobenius_map(2);

        t0 == t1
    }
}

impl Field for Fq12 {
    fn random<R: RngCore>(rng: &mut R) -> Self {
        Fq12 {
            c0: Fq6::random(rng),
            c1: Fq6::random(rng),
        }
    }

    fn zero() -> Self {
        Fq12 {
            c0: Fq6::zero(),
            c1: Fq6::zero(),
        }
    }

    fn one() -> Self {
        Fq12 {
            c0: Fq6::one(),
            c1: Fq6::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
    }

    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
    }

    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }

    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }

    fn frobenius_map(&mut self, power: usize) {
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);

        self.c1.c0.mul_assign(&FROBENIUS_COEFF_FQ12_C1[power % 12]);
        self.c1.c1.mul_assign(&FROBENIUS_COEFF_FQ12_C1[power % 12]);
        self.c1.c2.mul_assign(&FROBENIUS_COEFF_FQ12_C1[power % 12]);
    }

    fn square(&mut self) {
        let mut ab = self.c0;
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0;
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1;
        c0.mul_by_nonresidue();
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab;
        self.c1.add_assign(&ab);
        ab.mul_by_nonresidue();
        c0.sub_assign(&ab);
        self.c0 = c0;
    }

    fn mul_assign(&mut self, other: &Self) {
        let mut aa = self.c0;
        aa.mul_assign(&other.c0);
        let mut bb = self.c1;
        bb.mul_assign(&other.c1);
        let mut o = other.c0;
        o.add_assign(&other.c1);
        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = bb;
        self.c0.mul_by_nonresidue();
        self.c0.add_assign(&aa);
    }

    fn inverse(&self) -> Option<Self> {
        let mut c0s = self.c0;
        c0s.square();
        let mut c1s = self.c1;
        c1s.square();
        c1s.mul_by_nonresidue();
        c0s.sub_assign(&c1s);

        c0s.inverse().map(|t| {
            let mut tmp = Fq12 { c0: t, c1: t };
            tmp.c0.mul_assign(&self.c0);
            tmp.c1.mul_assign(&self.c1);
            tmp.c1.negate();

            tmp
        })
    }
}

/// Compressed representation of `Fq12`.
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Fq12Compressed(Fq6);

impl fmt::Debug for Fq12Compressed {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Fq12Compressed")
            .field("c0", &self.0)
            .finish()
    }
}

impl Fq12Compressed {
    /// Uncompress the given Fq12 element, returns `None` if the element is an invalid compression
    /// format.
    pub fn uncompress(self) -> Option<Fq12> {
        // Formula for decompression for the odd q case from Section 2 in
        // "Compression in finite fields and torus-based cryptography" by
        // Rubin-Silverberg.
        let mut fp6_neg_one = Fq6::one();
        fp6_neg_one.negate();
        let t = Fq12 {
            c0: self.0,
            c1: fp6_neg_one,
        }
        .inverse()
        .unwrap();
        let mut c = Fq12 {
            c0: self.0,
            c1: Fq6::one(),
        };
        c.mul_assign(&t);

        if c.is_cyc() {
            return Some(c);
        }

        None
    }
}

impl crate::Compress for Fq12 {
    fn write_compressed<W: std::io::Write>(self, mut out: W) -> std::io::Result<()> {
        let c = self.compress().unwrap();

        c.0.c0.c0.into_repr().write_le(&mut out)?;
        c.0.c0.c1.into_repr().write_le(&mut out)?;

        c.0.c1.c0.into_repr().write_le(&mut out)?;
        c.0.c1.c1.into_repr().write_le(&mut out)?;

        c.0.c2.c0.into_repr().write_le(&mut out)?;
        c.0.c2.c1.into_repr().write_le(&mut out)?;

        Ok(())
    }

    fn read_compressed<R: std::io::Read>(mut source: R) -> std::io::Result<Self> {
        let read_fp = |source: &mut dyn std::io::Read| {
            let mut repr = FqRepr::default();
            repr.read_le(source)?;
            let fp = Fq::from_repr(repr);
            fp.map_err(|err| std::io::Error::new(std::io::ErrorKind::InvalidData, err.to_string()))
        };

        let x0 = read_fp(&mut source)?;
        let x1 = read_fp(&mut source)?;

        let y0 = read_fp(&mut source)?;
        let y1 = read_fp(&mut source)?;

        let z0 = read_fp(&mut source)?;
        let z1 = read_fp(&mut source)?;

        let x = Fq2 { c0: x0, c1: x1 };
        let y = Fq2 { c0: y0, c1: y1 };
        let z = Fq2 { c0: z0, c1: z1 };

        let compressed = Fq12Compressed(Fq6 {
            c0: x,
            c1: y,
            c2: z,
        });
        compressed.uncompress().ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "invalid compression point")
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::bls12_381::{Bls12, Fq, G1, G2};
    use crate::Engine;
    use fff::PrimeField;
    use groupy::CurveProjective;
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_fq12_mul_by_014() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            let c0 = Fq2::random(&mut rng);
            let c1 = Fq2::random(&mut rng);
            let c5 = Fq2::random(&mut rng);
            let mut a = Fq12::random(&mut rng);
            let mut b = a;

            a.mul_by_014(&c0, &c1, &c5);
            b.mul_assign(&Fq12 {
                c0: Fq6 {
                    c0,
                    c1,
                    c2: Fq2::zero(),
                },
                c1: Fq6 {
                    c0: Fq2::zero(),
                    c1: c5,
                    c2: Fq2::zero(),
                },
            });

            assert_eq!(a, b);
        }
    }

    #[test]
    fn fq12_field_tests() {
        crate::tests::field::random_field_tests::<Fq12>();
        crate::tests::field::random_frobenius_tests::<Fq12, _>(Fq::char(), 13);
    }

    #[test]
    fn fp12_compression_test() {
        use crate::Compress;

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for i in 0..100 {
            let a = Fq12::random(&mut rng);
            // usually not cyclomatic, so not compressable
            if let Some(b) = a.compress() {
                let c = b.uncompress().unwrap();
                assert_eq!(a, c, "{}", i);
            } else {
                println!("skipping {}", i);
            }

            // pairing result, should be compressable
            let p = G1::random(&mut rng).into_affine();
            let q = G2::random(&mut rng).into_affine();
            let a = Bls12::pairing(p, q);
            assert!(a.is_cyc());

            let b = a.compress().unwrap();
            let c = b.uncompress().unwrap();
            assert_eq!(a, c, "{}", i);

            let mut buffer = Vec::new();
            a.write_compressed(&mut buffer).unwrap();
            let out = Fq12::read_compressed(std::io::Cursor::new(buffer)).unwrap();
            assert_eq!(a, out);
        }
    }
}
