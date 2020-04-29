use std::fmt;

use fff::{BitIterator, Field, PrimeField, PrimeFieldRepr, SqrtField};
use groupy::{CurveAffine, CurveProjective, EncodedPoint, GroupDecodingError};
use rand_core::RngCore;

use super::super::{Bls12, Fq, Fq12, FqRepr, Fr, FrRepr, IsogenyMap, OsswuMap};
use super::chain::chain_pm3div4;
use super::g2::G2Affine;
use super::util::osswu_helper;
use crate::{Engine, PairingCurveAffine, Signum0};

curve_impl!(
    "G1",
    G1,
    G1Affine,
    G1Prepared,
    Fq,
    Fr,
    G1Uncompressed,
    G1Compressed,
    G2Affine,
    16,
    15
);

#[derive(Copy, Clone)]
pub struct G1Uncompressed([u8; 96]);

encoded_point_delegations!(G1Uncompressed);

impl fmt::Debug for G1Uncompressed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        self.0[..].fmt(formatter)
    }
}

impl EncodedPoint for G1Uncompressed {
    type Affine = G1Affine;

    fn empty() -> Self {
        G1Uncompressed([0; 96])
    }
    fn size() -> usize {
        96
    }
    fn into_affine(&self) -> Result<G1Affine, GroupDecodingError> {
        let affine = self.into_affine_unchecked()?;

        if !affine.is_on_curve() {
            Err(GroupDecodingError::NotOnCurve)
        } else if !affine.is_in_correct_subgroup_assuming_on_curve() {
            Err(GroupDecodingError::NotInSubgroup)
        } else {
            Ok(affine)
        }
    }
    fn into_affine_unchecked(&self) -> Result<G1Affine, GroupDecodingError> {
        // Create a copy of this representation.
        let mut copy = self.0;

        if copy[0] & (1 << 7) != 0 {
            // Distinguisher bit is set, but this should be uncompressed!
            return Err(GroupDecodingError::UnexpectedCompressionMode);
        }

        if copy[0] & (1 << 6) != 0 {
            // This is the point at infinity, which means that if we mask away
            // the first two bits, the entire representation should consist
            // of zeroes.
            copy[0] &= 0x3f;

            if copy.iter().all(|b| *b == 0) {
                Ok(G1Affine::zero())
            } else {
                Err(GroupDecodingError::UnexpectedInformation)
            }
        } else {
            if copy[0] & (1 << 5) != 0 {
                // The bit indicating the y-coordinate should be lexicographically
                // largest is set, but this is an uncompressed element.
                return Err(GroupDecodingError::UnexpectedInformation);
            }

            // Unset the three most significant bits.
            copy[0] &= 0x1f;

            let mut x = FqRepr([0; 6]);
            let mut y = FqRepr([0; 6]);

            {
                let mut reader = &copy[..];

                x.read_be(&mut reader).unwrap();
                y.read_be(&mut reader).unwrap();
            }

            Ok(G1Affine {
                x: Fq::from_repr(x)
                    .map_err(|e| GroupDecodingError::CoordinateDecodingError("x coordinate", e))?,
                y: Fq::from_repr(y)
                    .map_err(|e| GroupDecodingError::CoordinateDecodingError("y coordinate", e))?,
                infinity: false,
            })
        }
    }
    fn from_affine(affine: G1Affine) -> Self {
        let mut res = Self::empty();

        if affine.is_zero() {
            // Set the second-most significant bit to indicate this point
            // is at infinity.
            res.0[0] |= 1 << 6;
        } else {
            let mut writer = &mut res.0[..];

            affine.x.into_repr().write_be(&mut writer).unwrap();
            affine.y.into_repr().write_be(&mut writer).unwrap();
        }

        res
    }
}

#[derive(Copy, Clone)]
pub struct G1Compressed([u8; 48]);

encoded_point_delegations!(G1Compressed);

impl fmt::Debug for G1Compressed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        self.0[..].fmt(formatter)
    }
}

impl EncodedPoint for G1Compressed {
    type Affine = G1Affine;

    fn empty() -> Self {
        G1Compressed([0; 48])
    }
    fn size() -> usize {
        48
    }
    fn into_affine(&self) -> Result<G1Affine, GroupDecodingError> {
        let affine = self.into_affine_unchecked()?;

        // NB: Decompression guarantees that it is on the curve already.

        if !affine.is_in_correct_subgroup_assuming_on_curve() {
            Err(GroupDecodingError::NotInSubgroup)
        } else {
            Ok(affine)
        }
    }
    fn into_affine_unchecked(&self) -> Result<G1Affine, GroupDecodingError> {
        // Create a copy of this representation.
        let mut copy = self.0;

        if copy[0] & (1 << 7) == 0 {
            // Distinguisher bit isn't set.
            return Err(GroupDecodingError::UnexpectedCompressionMode);
        }

        if copy[0] & (1 << 6) != 0 {
            // This is the point at infinity, which means that if we mask away
            // the first two bits, the entire representation should consist
            // of zeroes.
            copy[0] &= 0x3f;

            if copy.iter().all(|b| *b == 0) {
                Ok(G1Affine::zero())
            } else {
                Err(GroupDecodingError::UnexpectedInformation)
            }
        } else {
            // Determine if the intended y coordinate must be greater
            // lexicographically.
            let greatest = copy[0] & (1 << 5) != 0;

            // Unset the three most significant bits.
            copy[0] &= 0x1f;

            let mut x = FqRepr([0; 6]);

            {
                let mut reader = &copy[..];

                x.read_be(&mut reader).unwrap();
            }

            // Interpret as Fq element.
            let x = Fq::from_repr(x)
                .map_err(|e| GroupDecodingError::CoordinateDecodingError("x coordinate", e))?;

            G1Affine::get_point_from_x(x, greatest).ok_or(GroupDecodingError::NotOnCurve)
        }
    }
    fn from_affine(affine: G1Affine) -> Self {
        let mut res = Self::empty();

        if affine.is_zero() {
            // Set the second-most significant bit to indicate this point
            // is at infinity.
            res.0[0] |= 1 << 6;
        } else {
            {
                let mut writer = &mut res.0[..];

                affine.x.into_repr().write_be(&mut writer).unwrap();
            }

            let mut negy = affine.y;
            negy.negate();

            // Set the third most significant bit if the correct y-coordinate
            // is lexicographically largest.
            if affine.y > negy {
                res.0[0] |= 1 << 5;
            }
        }

        // Set highest bit to distinguish this as a compressed element.
        res.0[0] |= 1 << 7;

        res
    }
}

impl G1Affine {
    fn scale_by_cofactor(&self) -> G1 {
        // G1 cofactor = (x - 1)^2 / 3  = 76329603384216526031706109802092473003
        let cofactor = BitIterator::new([0x8c00aaab0000aaab, 0x396c8c005555e156]);
        self.mul_bits(cofactor)
    }

    fn get_generator() -> Self {
        G1Affine {
            x: super::super::fq::G1_GENERATOR_X,
            y: super::super::fq::G1_GENERATOR_Y,
            infinity: false,
        }
    }

    fn get_coeff_b() -> Fq {
        super::super::fq::B_COEFF
    }

    fn perform_pairing(&self, other: &G2Affine) -> Fq12 {
        super::super::Bls12::pairing(*self, *other)
    }
}

impl G1 {
    fn empirical_recommended_wnaf_for_scalar(scalar: FrRepr) -> usize {
        let num_bits = scalar.num_bits() as usize;

        if num_bits >= 130 {
            4
        } else if num_bits >= 34 {
            3
        } else {
            2
        }
    }

    fn empirical_recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
        const RECOMMENDATIONS: [usize; 12] =
            [1, 3, 7, 20, 43, 120, 273, 563, 1630, 3128, 7933, 62569];

        let mut ret = 4;
        for r in &RECOMMENDATIONS {
            if num_scalars > *r {
                ret += 1;
            } else {
                break;
            }
        }

        ret
    }
}

#[derive(Clone, Debug)]
pub struct G1Prepared(pub(crate) G1Affine);

impl G1Prepared {
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    pub fn from_affine(p: G1Affine) -> Self {
        G1Prepared(p)
    }
}

/// Coefficients of the 11-isogeny x map's numerator
const XNUM: [Fq; 12] = [
    Fq(FqRepr([
        0x4d18b6f3af00131c,
        0x19fa219793fee28c,
        0x3f2885f1467f19ae,
        0x23dcea34f2ffb304,
        0xd15b58d2ffc00054,
        0x0913be200a20bef4,
    ])),
    Fq(FqRepr([
        0x898985385cdbbd8b,
        0x3c79e43cc7d966aa,
        0x1597e193f4cd233a,
        0x8637ef1e4d6623ad,
        0x11b22deed20d827b,
        0x07097bc5998784ad,
    ])),
    Fq(FqRepr([
        0xa542583a480b664b,
        0xfc7169c026e568c6,
        0x5ba2ef314ed8b5a6,
        0x5b5491c05102f0e7,
        0xdf6e99707d2a0079,
        0x0784151ed7605524,
    ])),
    Fq(FqRepr([
        0x494e212870f72741,
        0xab9be52fbda43021,
        0x26f5577994e34c3d,
        0x049dfee82aefbd60,
        0x65dadd7828505289,
        0x0e93d431ea011aeb,
    ])),
    Fq(FqRepr([
        0x90ee774bd6a74d45,
        0x7ada1c8a41bfb185,
        0x0f1a8953b325f464,
        0x104c24211be4805c,
        0x169139d319ea7a8f,
        0x09f20ead8e532bf6,
    ])),
    Fq(FqRepr([
        0x6ddd93e2f43626b7,
        0xa5482c9aa1ccd7bd,
        0x143245631883f4bd,
        0x2e0a94ccf77ec0db,
        0xb0282d480e56489f,
        0x18f4bfcbb4368929,
    ])),
    Fq(FqRepr([
        0x23c5f0c953402dfd,
        0x7a43ff6958ce4fe9,
        0x2c390d3d2da5df63,
        0xd0df5c98e1f9d70f,
        0xffd89869a572b297,
        0x1277ffc72f25e8fe,
    ])),
    Fq(FqRepr([
        0x79f4f0490f06a8a6,
        0x85f894a88030fd81,
        0x12da3054b18b6410,
        0xe2a57f6505880d65,
        0xbba074f260e400f1,
        0x08b76279f621d028,
    ])),
    Fq(FqRepr([
        0xe67245ba78d5b00b,
        0x8456ba9a1f186475,
        0x7888bff6e6b33bb4,
        0xe21585b9a30f86cb,
        0x05a69cdcef55feee,
        0x09e699dd9adfa5ac,
    ])),
    Fq(FqRepr([
        0x0de5c357bff57107,
        0x0a0db4ae6b1a10b2,
        0xe256bb67b3b3cd8d,
        0x8ad456574e9db24f,
        0x0443915f50fd4179,
        0x098c4bf7de8b6375,
    ])),
    Fq(FqRepr([
        0xe6b0617e7dd929c7,
        0xfe6e37d442537375,
        0x1dafdeda137a489e,
        0xe4efd1ad3f767ceb,
        0x4a51d8667f0fe1cf,
        0x054fdf4bbf1d821c,
    ])),
    Fq(FqRepr([
        0x72db2a50658d767b,
        0x8abf91faa257b3d5,
        0xe969d6833764ab47,
        0x464170142a1009eb,
        0xb14f01aadb30be2f,
        0x18ae6a856f40715d,
    ])),
];

/// Coefficients of the 11-isogeny x map's denominator
const XDEN: [Fq; 11] = [
    Fq(FqRepr([
        0xb962a077fdb0f945,
        0xa6a9740fefda13a0,
        0xc14d568c3ed6c544,
        0xb43fc37b908b133e,
        0x9c0b3ac929599016,
        0x0165aa6c93ad115f,
    ])),
    Fq(FqRepr([
        0x23279a3ba506c1d9,
        0x92cfca0a9465176a,
        0x3b294ab13755f0ff,
        0x116dda1c5070ae93,
        0xed4530924cec2045,
        0x083383d6ed81f1ce,
    ])),
    Fq(FqRepr([
        0x9885c2a6449fecfc,
        0x4a2b54ccd37733f0,
        0x17da9ffd8738c142,
        0xa0fba72732b3fafd,
        0xff364f36e54b6812,
        0x0f29c13c660523e2,
    ])),
    Fq(FqRepr([
        0xe349cc118278f041,
        0xd487228f2f3204fb,
        0xc9d325849ade5150,
        0x43a92bd69c15c2df,
        0x1c2c7844bc417be4,
        0x12025184f407440c,
    ])),
    Fq(FqRepr([
        0x587f65ae6acb057b,
        0x1444ef325140201f,
        0xfbf995e71270da49,
        0xccda066072436a42,
        0x7408904f0f186bb2,
        0x13b93c63edf6c015,
    ])),
    Fq(FqRepr([
        0xfb918622cd141920,
        0x4a4c64423ecaddb4,
        0x0beb232927f7fb26,
        0x30f94df6f83a3dc2,
        0xaeedd424d780f388,
        0x06cc402dd594bbeb,
    ])),
    Fq(FqRepr([
        0xd41f761151b23f8f,
        0x32a92465435719b3,
        0x64f436e888c62cb9,
        0xdf70a9a1f757c6e4,
        0x6933a38d5b594c81,
        0x0c6f7f7237b46606,
    ])),
    Fq(FqRepr([
        0x693c08747876c8f7,
        0x22c9850bf9cf80f0,
        0x8e9071dab950c124,
        0x89bc62d61c7baf23,
        0xbc6be2d8dad57c23,
        0x17916987aa14a122,
    ])),
    Fq(FqRepr([
        0x1be3ff439c1316fd,
        0x9965243a7571dfa7,
        0xc7f7f62962f5cd81,
        0x32c6aa9af394361c,
        0xbbc2ee18e1c227f4,
        0x0c102cbac531bb34,
    ])),
    Fq(FqRepr([
        0x997614c97bacbf07,
        0x61f86372b99192c0,
        0x5b8c95fc14353fc3,
        0xca2b066c2a87492f,
        0x16178f5bbf698711,
        0x12a6dcd7f0f4e0e8,
    ])),
    Fq(FqRepr([
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493,
    ])),
];

/// Coefficients of the 11-isogeny y map's numerator
const YNUM: [Fq; 16] = [
    Fq(FqRepr([
        0x2b567ff3e2837267,
        0x1d4d9e57b958a767,
        0xce028fea04bd7373,
        0xcc31a30a0b6cd3df,
        0x7d7b18a682692693,
        0x0d300744d42a0310,
    ])),
    Fq(FqRepr([
        0x99c2555fa542493f,
        0xfe7f53cc4874f878,
        0x5df0608b8f97608a,
        0x14e03832052b49c8,
        0x706326a6957dd5a4,
        0x0a8dadd9c2414555,
    ])),
    Fq(FqRepr([
        0x13d942922a5cf63a,
        0x357e33e36e261e7d,
        0xcf05a27c8456088d,
        0x0000bd1de7ba50f0,
        0x83d0c7532f8c1fde,
        0x13f70bf38bbf2905,
    ])),
    Fq(FqRepr([
        0x5c57fd95bfafbdbb,
        0x28a359a65e541707,
        0x3983ceb4f6360b6d,
        0xafe19ff6f97e6d53,
        0xb3468f4550192bf7,
        0x0bb6cde49d8ba257,
    ])),
    Fq(FqRepr([
        0x590b62c7ff8a513f,
        0x314b4ce372cacefd,
        0x6bef32ce94b8a800,
        0x6ddf84a095713d5f,
        0x64eace4cb0982191,
        0x0386213c651b888d,
    ])),
    Fq(FqRepr([
        0xa5310a31111bbcdd,
        0xa14ac0f5da148982,
        0xf9ad9cc95423d2e9,
        0xaa6ec095283ee4a7,
        0xcf5b1f022e1c9107,
        0x01fddf5aed881793,
    ])),
    Fq(FqRepr([
        0x65a572b0d7a7d950,
        0xe25c2d8183473a19,
        0xc2fcebe7cb877dbd,
        0x05b2d36c769a89b0,
        0xba12961be86e9efb,
        0x07eb1b29c1dfde1f,
    ])),
    Fq(FqRepr([
        0x93e09572f7c4cd24,
        0x364e929076795091,
        0x8569467e68af51b5,
        0xa47da89439f5340f,
        0xf4fa918082e44d64,
        0x0ad52ba3e6695a79,
    ])),
    Fq(FqRepr([
        0x911429844e0d5f54,
        0xd03f51a3516bb233,
        0x3d587e5640536e66,
        0xfa86d2a3a9a73482,
        0xa90ed5adf1ed5537,
        0x149c9c326a5e7393,
    ])),
    Fq(FqRepr([
        0x462bbeb03c12921a,
        0xdc9af5fa0a274a17,
        0x9a558ebde836ebed,
        0x649ef8f11a4fae46,
        0x8100e1652b3cdc62,
        0x1862bd62c291dacb,
    ])),
    Fq(FqRepr([
        0x05c9b8ca89f12c26,
        0x0194160fa9b9ac4f,
        0x6a643d5a6879fa2c,
        0x14665bdd8846e19d,
        0xbb1d0d53af3ff6bf,
        0x12c7e1c3b28962e5,
    ])),
    Fq(FqRepr([
        0xb55ebf900b8a3e17,
        0xfedc77ec1a9201c4,
        0x1f07db10ea1a4df4,
        0x0dfbd15dc41a594d,
        0x389547f2334a5391,
        0x02419f98165871a4,
    ])),
    Fq(FqRepr([
        0xb416af000745fc20,
        0x8e563e9d1ea6d0f5,
        0x7c763e17763a0652,
        0x01458ef0159ebbef,
        0x8346fe421f96bb13,
        0x0d2d7b829ce324d2,
    ])),
    Fq(FqRepr([
        0x93096bb538d64615,
        0x6f2a2619951d823a,
        0x8f66b3ea59514fa4,
        0xf563e63704f7092f,
        0x724b136c4cf2d9fa,
        0x046959cfcfd0bf49,
    ])),
    Fq(FqRepr([
        0xea748d4b6e405346,
        0x91e9079c2c02d58f,
        0x41064965946d9b59,
        0xa06731f1d2bbe1ee,
        0x07f897e267a33f1b,
        0x1017290919210e5f,
    ])),
    Fq(FqRepr([
        0x872aa6c17d985097,
        0xeecc53161264562a,
        0x07afe37afff55002,
        0x54759078e5be6838,
        0xc4b92d15db8acca8,
        0x106d87d1b51d13b9,
    ])),
];

/// Coefficients of the 11-isogeny y map's denominator
const YDEN: [Fq; 16] = [
    Fq(FqRepr([
        0xeb6c359d47e52b1c,
        0x18ef5f8a10634d60,
        0xddfa71a0889d5b7e,
        0x723e71dcc5fc1323,
        0x52f45700b70d5c69,
        0x0a8b981ee47691f1,
    ])),
    Fq(FqRepr([
        0x616a3c4f5535b9fb,
        0x6f5f037395dbd911,
        0xf25f4cc5e35c65da,
        0x3e50dffea3c62658,
        0x6a33dca523560776,
        0x0fadeff77b6bfe3e,
    ])),
    Fq(FqRepr([
        0x2be9b66df470059c,
        0x24a2c159a3d36742,
        0x115dbe7ad10c2a37,
        0xb6634a652ee5884d,
        0x04fe8bb2b8d81af4,
        0x01c2a7a256fe9c41,
    ])),
    Fq(FqRepr([
        0xf27bf8ef3b75a386,
        0x898b367476c9073f,
        0x24482e6b8c2f4e5f,
        0xc8e0bbd6fe110806,
        0x59b0c17f7631448a,
        0x11037cd58b3dbfbd,
    ])),
    Fq(FqRepr([
        0x31c7912ea267eec6,
        0x1dbf6f1c5fcdb700,
        0xd30d4fe3ba86fdb1,
        0x3cae528fbee9a2a4,
        0xb1cce69b6aa9ad9a,
        0x044393bb632d94fb,
    ])),
    Fq(FqRepr([
        0xc66ef6efeeb5c7e8,
        0x9824c289dd72bb55,
        0x71b1a4d2f119981d,
        0x104fc1aafb0919cc,
        0x0e49df01d942a628,
        0x096c3a09773272d4,
    ])),
    Fq(FqRepr([
        0x9abc11eb5fadeff4,
        0x32dca50a885728f0,
        0xfb1fa3721569734c,
        0xc4b76271ea6506b3,
        0xd466a75599ce728e,
        0x0c81d4645f4cb6ed,
    ])),
    Fq(FqRepr([
        0x4199f10e5b8be45b,
        0xda64e495b1e87930,
        0xcb353efe9b33e4ff,
        0x9e9efb24aa6424c6,
        0xf08d33680a237465,
        0x0d3378023e4c7406,
    ])),
    Fq(FqRepr([
        0x7eb4ae92ec74d3a5,
        0xc341b4aa9fac3497,
        0x5be603899e907687,
        0x03bfd9cca75cbdeb,
        0x564c2935a96bfa93,
        0x0ef3c33371e2fdb5,
    ])),
    Fq(FqRepr([
        0x7ee91fd449f6ac2e,
        0xe5d5bd5cb9357a30,
        0x773a8ca5196b1380,
        0xd0fda172174ed023,
        0x6cb95e0fa776aead,
        0x0d22d5a40cec7cff,
    ])),
    Fq(FqRepr([
        0xf727e09285fd8519,
        0xdc9d55a83017897b,
        0x7549d8bd057894ae,
        0x178419613d90d8f8,
        0xfce95ebdeb5b490a,
        0x0467ffaef23fc49e,
    ])),
    Fq(FqRepr([
        0xc1769e6a7c385f1b,
        0x79bc930deac01c03,
        0x5461c75a23ede3b5,
        0x6e20829e5c230c45,
        0x828e0f1e772a53cd,
        0x116aefa749127bff,
    ])),
    Fq(FqRepr([
        0x101c10bf2744c10a,
        0xbbf18d053a6a3154,
        0xa0ecf39ef026f602,
        0xfc009d4996dc5153,
        0xb9000209d5bd08d3,
        0x189e5fe4470cd73c,
    ])),
    Fq(FqRepr([
        0x7ebd546ca1575ed2,
        0xe47d5a981d081b55,
        0x57b2b625b6d4ca21,
        0xb0a1ba04228520cc,
        0x98738983c2107ff3,
        0x13dddbc4799d81d6,
    ])),
    Fq(FqRepr([
        0x09319f2e39834935,
        0x039e952cbdb05c21,
        0x55ba77a9a2f76493,
        0xfd04e3dfc6086467,
        0xfb95832e7d78742e,
        0x0ef9c24eccaf5e0e,
    ])),
    Fq(FqRepr([
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493,
    ])),
];

const ELLP_A: Fq = Fq(FqRepr([
    0x2f65aa0e9af5aa51u64,
    0x86464c2d1e8416c3u64,
    0xb85ce591b7bd31e2u64,
    0x27e11c91b5f24e7cu64,
    0x28376eda6bfc1835u64,
    0x155455c3e5071d85u64,
]));

const ELLP_B: Fq = Fq(FqRepr([
    0xfb996971fe22a1e0u64,
    0x9aa93eb35b742d6fu64,
    0x8c476013de99c5c4u64,
    0x873e27c3a221e571u64,
    0xca72b5e45a52d888u64,
    0x06824061418a386bu64,
]));

const XI: Fq = Fq(FqRepr([
    0x886c00000023ffdcu64,
    0xf70008d3090001du64,
    0x77672417ed5828c3u64,
    0x9dac23e943dc1740u64,
    0x50553f1b9c131521u64,
    0x78c712fbe0ab6e8u64,
]));

const SQRT_M_XI_CUBED: Fq = Fq(FqRepr([
    0x43b571cad3215f1fu64,
    0xccb460ef1c702dc2u64,
    0x742d884f4f97100bu64,
    0xdb2c3e3238a3382bu64,
    0xe40f3fa13fce8f88u64,
    0x73a2af9892a2ffu64,
]));

impl IsogenyMap for G1 {
    fn isogeny_map(&mut self) {
        self.eval_iso([&XNUM[..], &XDEN[..], &YNUM[..], &YDEN[..]]);
    }
}

impl OsswuMap for G1 {
    fn osswu_map(u: &Fq) -> G1 {
        // compute x0 and g(x0)
        let [usq, xi_usq, _, x0_num, x0_den, gx0_num, gx0_den] =
            osswu_helper(u, &XI, &ELLP_A, &ELLP_B);

        // compute g(X0(u)) ^ ((p - 3) // 4)
        let sqrt_candidate = {
            let mut tmp1 = gx0_num;
            tmp1.mul_assign(&gx0_den); // u * v
            let mut tmp2 = gx0_den;
            tmp2.square(); // v^2
            tmp2.mul_assign(&tmp1); // u * v^3
            let tmp3 = tmp2;
            chain_pm3div4(&mut tmp2, &tmp3); // (u v^3) ^ ((p - 3) // 4)
            tmp2.mul_assign(&tmp1); // u v (u v^3) ^ ((p - 3) // 4)
            tmp2
        };

        // select correct values for y and for x numerator
        let (mut x_num, mut y) = {
            let mut test_cand = sqrt_candidate;
            test_cand.square();
            test_cand.mul_assign(&gx0_den);
            if test_cand == gx0_num {
                (x0_num, sqrt_candidate) // g(x0) is square
            } else {
                let mut x1_num = x0_num; // g(x1) is square
                x1_num.mul_assign(&xi_usq); // x1 = xi u^2 g(x0)
                let mut y1 = usq; // y1 = sqrt(-xi**3) * u^3 g(x0) ^ ((p - 1) // 4)
                y1.mul_assign(&u);
                y1.mul_assign(&sqrt_candidate);
                y1.mul_assign(&SQRT_M_XI_CUBED);
                (x1_num, y1)
            }
        };

        // make sure sign of y and sign of u agree
        let sgn0_y_xor_u = y.sgn0() ^ u.sgn0();
        y.negate_if(sgn0_y_xor_u);

        // convert to projective
        x_num.mul_assign(&x0_den); // x_num * x_den / x_den^2 = x_num / x_den
        y.mul_assign(&gx0_den); // y * x_den^3 / x_den^3 = y

        G1 {
            x: x_num,
            y,
            z: x0_den,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::util::check_g_prime;

    #[test]
    fn g1_generator() {
        use crate::SqrtField;

        let mut x = Fq::zero();
        let mut i = 0;
        loop {
            // y^2 = x^3 + b
            let mut rhs = x;
            rhs.square();
            rhs.mul_assign(&x);
            rhs.add_assign(&G1Affine::get_coeff_b());

            if let Some(y) = rhs.sqrt() {
                let yrepr = y.into_repr();
                let mut negy = y;
                negy.negate();
                let negyrepr = negy.into_repr();

                let p = G1Affine {
                    x,
                    y: if yrepr < negyrepr { y } else { negy },
                    infinity: false,
                };
                assert!(!p.is_in_correct_subgroup_assuming_on_curve());

                let g1 = p.scale_by_cofactor();
                if !g1.is_zero() {
                    assert_eq!(i, 4);
                    let g1 = G1Affine::from(g1);

                    assert!(g1.is_in_correct_subgroup_assuming_on_curve());

                    assert_eq!(g1, G1Affine::one());
                    break;
                }
            }

            i += 1;
            x.add_assign(&Fq::one());
        }
    }

    #[test]
    fn g1_test_is_valid() {
        // Reject point on isomorphic twist (b = 24)
        {
            let p = G1Affine {
                x: Fq::from_repr(FqRepr([
                    0xc58d887b66c035dc,
                    0x10cbfd301d553822,
                    0xaf23e064f1131ee5,
                    0x9fe83b1b4a5d648d,
                    0xf583cc5a508f6a40,
                    0xc3ad2aefde0bb13,
                ]))
                .unwrap(),
                y: Fq::from_repr(FqRepr([
                    0x60aa6f9552f03aae,
                    0xecd01d5181300d35,
                    0x8af1cdb8aa8ce167,
                    0xe760f57922998c9d,
                    0x953703f5795a39e5,
                    0xfe3ae0922df702c,
                ]))
                .unwrap(),
                infinity: false,
            };
            assert!(!p.is_on_curve());
            assert!(p.is_in_correct_subgroup_assuming_on_curve());
        }

        // Reject point on a twist (b = 3)
        {
            let p = G1Affine {
                x: Fq::from_repr(FqRepr([
                    0xee6adf83511e15f5,
                    0x92ddd328f27a4ba6,
                    0xe305bd1ac65adba7,
                    0xea034ee2928b30a8,
                    0xbd8833dc7c79a7f7,
                    0xe45c9f0c0438675,
                ]))
                .unwrap(),
                y: Fq::from_repr(FqRepr([
                    0x3b450eb1ab7b5dad,
                    0xa65cb81e975e8675,
                    0xaa548682b21726e5,
                    0x753ddf21a2601d20,
                    0x532d0b640bd3ff8b,
                    0x118d2c543f031102,
                ]))
                .unwrap(),
                infinity: false,
            };
            assert!(!p.is_on_curve());
            assert!(!p.is_in_correct_subgroup_assuming_on_curve());
        }

        // Reject point in an invalid subgroup
        // There is only one r-order subgroup, as r does not divide the cofactor.
        {
            let p = G1Affine {
                x: Fq::from_repr(FqRepr([
                    0x76e1c971c6db8fe8,
                    0xe37e1a610eff2f79,
                    0x88ae9c499f46f0c0,
                    0xf35de9ce0d6b4e84,
                    0x265bddd23d1dec54,
                    0x12a8778088458308,
                ]))
                .unwrap(),
                y: Fq::from_repr(FqRepr([
                    0x8a22defa0d526256,
                    0xc57ca55456fcb9ae,
                    0x1ba194e89bab2610,
                    0x921beef89d4f29df,
                    0x5b6fda44ad85fa78,
                    0xed74ab9f302cbe0,
                ]))
                .unwrap(),
                infinity: false,
            };
            assert!(p.is_on_curve());
            assert!(!p.is_in_correct_subgroup_assuming_on_curve());
        }
    }

    #[test]
    fn test_g1_addition_correctness() {
        let mut p = G1 {
            x: Fq::from_repr(FqRepr([
                0x47fd1f891d6e8bbf,
                0x79a3b0448f31a2aa,
                0x81f3339e5f9968f,
                0x485e77d50a5df10d,
                0x4c6fcac4b55fd479,
                0x86ed4d9906fb064,
            ]))
            .unwrap(),
            y: Fq::from_repr(FqRepr([
                0xd25ee6461538c65,
                0x9f3bbb2ecd3719b9,
                0xa06fd3f1e540910d,
                0xcefca68333c35288,
                0x570c8005f8573fa6,
                0x152ca696fe034442,
            ]))
            .unwrap(),
            z: Fq::one(),
        };

        p.add_assign(&G1 {
            x: Fq::from_repr(FqRepr([
                0xeec78f3096213cbf,
                0xa12beb1fea1056e6,
                0xc286c0211c40dd54,
                0x5f44314ec5e3fb03,
                0x24e8538737c6e675,
                0x8abd623a594fba8,
            ]))
            .unwrap(),
            y: Fq::from_repr(FqRepr([
                0x6b0528f088bb7044,
                0x2fdeb5c82917ff9e,
                0x9a5181f2fac226ad,
                0xd65104c6f95a872a,
                0x1f2998a5a9c61253,
                0xe74846154a9e44,
            ]))
            .unwrap(),
            z: Fq::one(),
        });

        let p = G1Affine::from(p);

        assert_eq!(
            p,
            G1Affine {
                x: Fq::from_repr(FqRepr([
                    0x6dd3098f22235df,
                    0xe865d221c8090260,
                    0xeb96bb99fa50779f,
                    0xc4f9a52a428e23bb,
                    0xd178b28dd4f407ef,
                    0x17fb8905e9183c69
                ]))
                .unwrap(),
                y: Fq::from_repr(FqRepr([
                    0xd0de9d65292b7710,
                    0xf6a05f2bcf1d9ca7,
                    0x1040e27012f20b64,
                    0xeec8d1a5b7466c58,
                    0x4bc362649dce6376,
                    0x430cbdc5455b00a
                ]))
                .unwrap(),
                infinity: false,
            }
        );
    }

    #[test]
    fn test_g1_doubling_correctness() {
        let mut p = G1 {
            x: Fq::from_repr(FqRepr([
                0x47fd1f891d6e8bbf,
                0x79a3b0448f31a2aa,
                0x81f3339e5f9968f,
                0x485e77d50a5df10d,
                0x4c6fcac4b55fd479,
                0x86ed4d9906fb064,
            ]))
            .unwrap(),
            y: Fq::from_repr(FqRepr([
                0xd25ee6461538c65,
                0x9f3bbb2ecd3719b9,
                0xa06fd3f1e540910d,
                0xcefca68333c35288,
                0x570c8005f8573fa6,
                0x152ca696fe034442,
            ]))
            .unwrap(),
            z: Fq::one(),
        };

        p.double();

        let p = G1Affine::from(p);

        assert_eq!(
            p,
            G1Affine {
                x: Fq::from_repr(FqRepr([
                    0xf939ddfe0ead7018,
                    0x3b03942e732aecb,
                    0xce0e9c38fdb11851,
                    0x4b914c16687dcde0,
                    0x66c8baf177d20533,
                    0xaf960cff3d83833
                ]))
                .unwrap(),
                y: Fq::from_repr(FqRepr([
                    0x3f0675695f5177a8,
                    0x2b6d82ae178a1ba0,
                    0x9096380dd8e51b11,
                    0x1771a65b60572f4e,
                    0x8b547c1313b27555,
                    0x135075589a687b1e
                ]))
                .unwrap(),
                infinity: false,
            }
        );
    }

    #[test]
    fn test_g1_same_y() {
        // Test the addition of two points with different x coordinates
        // but the same y coordinate.

        // x1 = 128100205326445210408953809171070606737678357140298133325128175840781723996595026100005714405541449960643523234125
        // x2 = 3821408151224848222394078037104966877485040835569514006839342061575586899845797797516352881516922679872117658572470
        // y = 2291134451313223670499022936083127939567618746216464377735567679979105510603740918204953301371880765657042046687078

        let a = G1Affine {
            x: Fq::from_repr(FqRepr([
                0xea431f2cc38fc94d,
                0x3ad2354a07f5472b,
                0xfe669f133f16c26a,
                0x71ffa8021531705,
                0x7418d484386d267,
                0xd5108d8ff1fbd6,
            ]))
            .unwrap(),
            y: Fq::from_repr(FqRepr([
                0xa776ccbfe9981766,
                0x255632964ff40f4a,
                0xc09744e650b00499,
                0x520f74773e74c8c3,
                0x484c8fc982008f0,
                0xee2c3d922008cc6,
            ]))
            .unwrap(),
            infinity: false,
        };

        let b = G1Affine {
            x: Fq::from_repr(FqRepr([
                0xe06cdb156b6356b6,
                0xd9040b2d75448ad9,
                0xe702f14bb0e2aca5,
                0xc6e05201e5f83991,
                0xf7c75910816f207c,
                0x18d4043e78103106,
            ]))
            .unwrap(),
            y: Fq::from_repr(FqRepr([
                0xa776ccbfe9981766,
                0x255632964ff40f4a,
                0xc09744e650b00499,
                0x520f74773e74c8c3,
                0x484c8fc982008f0,
                0xee2c3d922008cc6,
            ]))
            .unwrap(),
            infinity: false,
        };

        // Expected
        // x = 52901198670373960614757979459866672334163627229195745167587898707663026648445040826329033206551534205133090753192
        // y = 1711275103908443722918766889652776216989264073722543507596490456144926139887096946237734327757134898380852225872709
        let c = G1Affine {
            x: Fq::from_repr(FqRepr([
                0xef4f05bdd10c8aa8,
                0xad5bf87341a2df9,
                0x81c7424206b78714,
                0x9676ff02ec39c227,
                0x4c12c15d7e55b9f3,
                0x57fd1e317db9bd,
            ]))
            .unwrap(),
            y: Fq::from_repr(FqRepr([
                0x1288334016679345,
                0xf955cd68615ff0b5,
                0xa6998dbaa600f18a,
                0x1267d70db51049fb,
                0x4696deb9ab2ba3e7,
                0xb1e4e11177f59d4,
            ]))
            .unwrap(),
            infinity: false,
        };

        assert!(a.is_on_curve() && a.is_in_correct_subgroup_assuming_on_curve());
        assert!(b.is_on_curve() && b.is_in_correct_subgroup_assuming_on_curve());
        assert!(c.is_on_curve() && c.is_in_correct_subgroup_assuming_on_curve());

        let mut tmp1 = a.into_projective();
        tmp1.add_assign(&b.into_projective());
        assert_eq!(tmp1.into_affine(), c);
        assert_eq!(tmp1, c.into_projective());

        let mut tmp2 = a.into_projective();
        tmp2.add_assign_mixed(&b);
        assert_eq!(tmp2.into_affine(), c);
        assert_eq!(tmp2, c.into_projective());
    }

    #[test]
    fn test_g1_sw_encode_degenerate() {
        // test the degenerate case t = 0
        let p = G1Affine::sw_encode(Fq::zero());
        assert!(p.is_on_curve());
        assert!(p.is_zero());

        // test the degenerate case t^2 = - b - 1
        let mut t = Fq::one();
        t.add_assign(&G1Affine::get_coeff_b());
        t.negate();
        let mut t = t.sqrt().unwrap();
        t.negate(); // If sqrt impl changes, this test will be affected
        let p = G1Affine::sw_encode(t);
        assert!(p.is_on_curve());
        assert!(!p.is_zero());
        assert_eq!(p.y.parity(), t.parity());
        assert_eq!(p, G1Affine::one());
        t.negate();
        let p = G1Affine::sw_encode(t);
        assert!(p.is_on_curve());
        assert!(!p.is_zero());
        assert_eq!(p.y.parity(), t.parity());
        {
            let mut negone = G1Affine::one();
            negone.negate();
            assert_eq!(p, negone);
        }

        // test that the encoding function is odd for the above t
        t.negate();
        let mut minus_p = G1Affine::sw_encode(t).into_projective();
        minus_p.add_assign_mixed(&p);
        assert!(minus_p.is_zero());
    }

    #[test]
    fn g1_hash_test_vectors() {
        // Obtained via python/sage

        let p = G1::hash(&[]);
        let q = G1 {
            x: Fq::from_str("315124130825307604287835216317628428134609737854237653839182597515996444073032649481416725367158979153513345579672").unwrap(),
            y: Fq::from_str("3093537746211397858160667262592024570071165158580434464756577567510401504168962073691924150397172185836012224315174").unwrap(),
            z: Fq::one()
        };

        assert_eq!(p, q);
    }

    #[test]
    fn g1_curve_tests() {
        use groupy::tests::curve_tests;
        curve_tests::<G1>();
    }

    #[test]
    fn test_iso11_zero() {
        let zero = Fq::zero();
        let mut pt = G1 {
            x: Fq::zero(),
            y: Fq::zero(),
            z: Fq::zero(),
        };
        pt.isogeny_map();
        assert_eq!(pt.x, zero);
        assert_eq!(pt.y, zero);
        assert_eq!(pt.z, zero);
    }

    #[test]
    fn test_iso11_one() {
        let mut pt = G1 {
            x: Fq::one(),
            y: Fq::one(),
            z: Fq::one(),
        };
        pt.isogeny_map();
        assert_eq!(
            pt.x,
            Fq::from_repr(FqRepr([
                0xb129fab9bef88eddu64,
                0x1c5429e2f4b8bc35u64,
                0xcaab8cc9ec4893f2u64,
                0x9e9c31f30a607c8bu64,
                0x9661fcf22bedfddbu64,
                0x10fc4a3ba5f48e07u64,
            ]))
            .unwrap()
        );
        assert_eq!(
            pt.y,
            Fq::from_repr(FqRepr([
                0xaf52c5fbd490f370u64,
                0x1533c0f27b46c02fu64,
                0xc8890dd0987b134fu64,
                0x43e2d5f172257d50u64,
                0x538ebef63fb145beu64,
                0x11eab1145b95cb9fu64,
            ]))
            .unwrap()
        );
        assert_eq!(
            pt.z,
            Fq::from_repr(FqRepr([
                0x7441c43513e11f49u64,
                0x620b0af2483ad30fu64,
                0x678c5bf3ad4090b4u64,
                0xc75152c6f387d070u64,
                0x5f3cc0ed1bd3f0eeu64,
                0x12514e630a486abbu64,
            ]))
            .unwrap()
        );
    }

    #[test]
    fn test_iso11_fixed() {
        let xi = Fq::from_repr(FqRepr([
            0xf6adc4118ae592abu64,
            0xa384a7ab165def35u64,
            0x2365b1fb1c8a73bfu64,
            0xc40dc338ca285231u64,
            0x47ff3364428c59b3u64,
            0x1789051238d025e3u64,
        ]))
        .unwrap();
        let yi = Fq::from_repr(FqRepr([
            0x1a635634e9cced27u64,
            0x03f604e47bc51aa9u64,
            0x06f6ff472fa7276eu64,
            0x0459ed10f1f8abb1u64,
            0x8e76c82bd4a29d21u64,
            0x088cb5712bf81924u64,
        ]))
        .unwrap();
        let zi = Fq::from_repr(FqRepr([
            0x0416411fe2e97d06u64,
            0xaced7fec7a63fe65u64,
            0x683295bcaed54202u64,
            0xbdc3405df9ff0a3bu64,
            0xf9698f57510273fbu64,
            0x064bb4b501466b2au64,
        ]))
        .unwrap();
        let mut pt = G1 {
            x: xi,
            y: yi,
            z: zi,
        };
        pt.isogeny_map();
        assert_eq!(
            pt.x,
            Fq::from_repr(FqRepr([
                0xa51741657e71601du64,
                0x1771cef34519b6f2u64,
                0x2d1aff4e4ae28379u64,
                0x9ddcd540391389adu64,
                0x0db61b8544450f53u64,
                0x0f34c6cea2fc0199u64,
            ]))
            .unwrap()
        );
        assert_eq!(
            pt.y,
            Fq::from_repr(FqRepr([
                0xd1d70b485ea22464u64,
                0xd3a592a3ffc2c77cu64,
                0x72ef2afff097ad4fu64,
                0xf1c66e0e000b5673u64,
                0x1d32499c9f462716u64,
                0x19284e38020f6072u64,
            ]))
            .unwrap()
        );
        assert_eq!(
            pt.z,
            Fq::from_repr(FqRepr([
                0x583946b46d152c9fu64,
                0xb7f34ad188fdc105u64,
                0x47f7edb38429108au64,
                0xb6602e02d0d7ac4du64,
                0xc27121d0eb3d5efcu64,
                0x16f4243bf7230576u64,
            ]))
            .unwrap()
        );
    }

    fn check_g1_prime(x: &Fq, y: &Fq, z: &Fq) {
        check_g_prime(x, y, z, &ELLP_A, &ELLP_B);
    }

    #[test]
    fn test_osswu_g1() {
        // exceptional case: zero
        let p = G1::osswu_map(&Fq::zero());
        let G1 { x, y, z } = &p;
        let xo = Fq::from_repr(FqRepr([
            0x6144f0e146df0250u64,
            0x9e9fd4264a7edcbau64,
            0x519289c2e473a9c7u64,
            0xfc9e9c179c1c484fu64,
            0x1bde5cc11dc20ba5u64,
            0x119d96b86f8b3b8bu64,
        ]))
        .unwrap();
        let yo = Fq::from_repr(FqRepr([
            0x2c26d31ff8057aa2u64,
            0x9f824897b954500eu64,
            0xd6b1bcf4165f3575u64,
            0x8d267d9b89fb2b31u64,
            0x905bde90d4b39d8au64,
            0x8327183f6473933u64,
        ]))
        .unwrap();
        let zo = Fq::from_repr(FqRepr([
            0xfe7db859f2cb453fu64,
            0x8e55cb15e9aab878u64,
            0x51fe89284e4d926au64,
            0x9a148b96ab3e6941u64,
            0xa3857e1ea7b2289du64,
            0xdf088f08f205e3u64,
        ]))
        .unwrap();
        assert_eq!(x, &xo);
        assert_eq!(y, &yo);
        assert_eq!(z, &zo);
        check_g1_prime(x, y, z);

        // exceptional case: sqrt(-1/XI) (positive)
        let excp = Fq::from_repr(FqRepr([
            0x7cc51062bde821b8u64,
            0x88b69520ee5c57fbu64,
            0x46edbdd403fc310u64,
            0x12f01df4948d09ffu64,
            0xdb38f4a9a3d71bdau64,
            0x1f7462c8b6cbf74u64,
        ]))
        .unwrap();
        let p = G1::osswu_map(&excp);
        let G1 { x, y, z } = &p;
        assert_eq!(x, &xo);
        assert_eq!(y, &yo);
        assert_eq!(z, &zo);
        check_g1_prime(x, y, z);

        // exceptional case: sqrt(-1/XI) (negative)
        let excp = Fq::from_repr(FqRepr([
            0x3d39ef9d421788f3u64,
            0x95f56addc2f7a804u64,
            0x62c1f6c3b6713313u64,
            0x51872d905ef808c0u64,
            0x6fe2b30c9f7490fdu64,
            0x1809cbbdae132725u64,
        ]))
        .unwrap();
        let p = G1::osswu_map(&excp);
        let G1 { x, y, z } = &p;
        let myo = {
            let mut tmp = yo;
            tmp.negate();
            tmp
        };
        assert_eq!(x, &xo);
        assert_eq!(y, &myo);
        assert_eq!(z, &zo);
        check_g1_prime(x, y, z);

        let u = Fq::from_repr(FqRepr([
            0xd4e2aa3bbf9a8255u64,
            0xa79f2ece3390978cu64,
            0x48c1a8fdff541ebau64,
            0x2b17303f8af1ec82u64,
            0x86657cd3fc3d08b5u64,
            0x14f05da1ad4eddc8u64,
        ]))
        .unwrap();
        let xo = Fq::from_repr(FqRepr([
            0xb8e5b32b10dd26f7u64,
            0x8a114aa4ef26ad27u64,
            0xad97709b49ae7c62u64,
            0x9bc765ec50b53945u64,
            0xae99d020a70ca4feu64,
            0x1803cbf9bd2e3815u64,
        ]))
        .unwrap();
        let yo = Fq::from_repr(FqRepr([
            0x498ec4b38b052163u64,
            0xdfb4b3c21c64a917u64,
            0xa6ad223eeba44938u64,
            0xa564373b4a3b1d49u64,
            0x4f3ba7671555ba8eu64,
            0x141f3b7a3a3bc9a1u64,
        ]))
        .unwrap();
        let zo = Fq::from_repr(FqRepr([
            0xc75f9dc8b69d09eeu64,
            0x80824ef4608083ceu64,
            0xfcd339725e80194au64,
            0xda50cf8999450757u64,
            0x35da50fd75b53f96u64,
            0xade87be1822999bu64,
        ]))
        .unwrap();
        let p = G1::osswu_map(&u);
        let G1 { x, y, z } = &p;
        assert_eq!(x, &xo);
        assert_eq!(y, &yo);
        assert_eq!(z, &zo);
        check_g1_prime(x, y, z);

        let mut rng = rand_xorshift::XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        for _ in 0..32 {
            let input = Fq::random(&mut rng);
            let p = G1::osswu_map(&input);
            let G1 { x, y, z } = &p;
            check_g1_prime(x, y, z);
        }
    }

    #[test]
    fn test_encode_to_curve_07() {
        use crate::{ExpandMsgXmd, HashToCurve};

        struct TestCase {
            msg: &'static [u8],
            expected: [&'static str; 2],
        }
        impl TestCase {
            fn expected(&self) -> String {
                self.expected[0].to_string() + self.expected[1]
            }
        }

        const DOMAIN: &[u8] = b"BLS12381G1_XMD:SHA-256_SSWU_NU_TESTGEN";

        let cases = vec![
            TestCase {
                msg: b"",
                expected: [
		    "1223effdbb2d38152495a864d78eee14cb0992d89a241707abb03819a91a6d2fd65854ab9a69e9aacb0cbebfd490732c",
		    "0f925d61e0b235ecd945cbf0309291878df0d06e5d80d6b84aa4ff3e00633b26f9a7cb3523ef737d90e6d71e8b98b2d5",
                ],
            },
            TestCase {
                msg: b"abc",
                expected: [
		    "179d3fd0b4fb1da43aad06cea1fb3f828806ddb1b1fa9424b1e3944dfdbab6e763c42636404017da03099af0dcca0fd6",
		    "0d037cb1c6d495c0f5f22b061d23f1be3d7fe64d3c6820cfcd99b6b36fa69f7b4c1f4addba2ae7aa46fb25901ab483e4",

                ],
            },
            TestCase {
                msg: b"abcdef0123456789",
                expected: [
		    "15aa66c77eded1209db694e8b1ba49daf8b686733afaa7b68c683d0b01788dfb0617a2e2d04c0856db4981921d3004af",
		    "0952bb2f61739dd1d201dd0a79d74cda3285403d47655ee886afe860593a8a4e51c5b77a22d2133e3a4280eaaaa8b788",
                ]
            },
            TestCase {
                msg: b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                expected: [
		    "06328ce5106e837935e8da84bd9af473422e62492930aa5f460369baad9545defa468d9399854c23a75495d2a80487ee",
		    "094bfdfe3e552447433b5a00967498a3f1314b86ce7a7164c8a8f4131f99333b30a574607e301d5f774172c627fd0bca",
                ]
            }
        ];

        for case in cases {
            let g = <G1 as HashToCurve<ExpandMsgXmd<sha2ni::Sha256>>>::encode_to_curve(
                &case.msg, DOMAIN,
            );
            let g_uncompressed = g.into_affine().into_uncompressed();

            assert_eq!(case.expected(), hex::encode(&g_uncompressed.0[..]));
        }
    }

    #[test]
    fn test_hash_to_curve_07() {
        use crate::{ExpandMsgXmd, HashToCurve};

        struct TestCase {
            msg: &'static [u8],
            expected: [&'static str; 2],
        }
        impl TestCase {
            fn expected(&self) -> String {
                self.expected[0].to_string() + self.expected[1]
            }
        }

        const DOMAIN: &[u8] = b"BLS12381G1_XMD:SHA-256_SSWU_RO_TESTGEN";

        let cases = vec![
            TestCase {
                msg: b"",
                expected: [
                    "0576730ab036cbac1d95b38dca905586f28d0a59048db4e8778782d89bff856ddef89277ead5a21e2975c4a6e3d8c79e",
                    "1273e568bebf1864393c517f999b87c1eaa1b8432f95aea8160cd981b5b05d8cd4a7cf00103b6ef87f728e4b547dd7ae",
                ],
            },
            TestCase {
                msg: b"abc",
                expected: [
                    "061daf0cc00d8912dac1d4cf5a7c32fca97f8b3bf3f805121888e5eb89f77f9a9f406569027ac6d0e61b1229f42c43d6",
                    "0de1601e5ba02cb637c1d35266f5700acee9850796dc88e860d022d7b9e7e3dce5950952e97861e5bb16d215c87f030d"
                ],
            },
            TestCase {
                msg: b"abcdef0123456789",
                expected: [
                    "0fb3455436843e76079c7cf3dfef75e5a104dfe257a29a850c145568d500ad31ccfe79be9ae0ea31a722548070cf98cd",
                    "177989f7e2c751658df1b26943ee829d3ebcf131d8f805571712f3a7527ee5334ecff8a97fc2a50cea86f5e6212e9a57"
                ]
            },
            TestCase {
                msg: b"a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
                       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                expected: [
                    "0514af2137c1ae1d78d5cb97ee606ea142824c199f0f25ac463a0c78200de57640d34686521d3e9cf6b3721834f8a038",
                    "047a85d6898416a0899e26219bca7c4f0fa682717199de196b02b95eaf9fb55456ac3b810e78571a1b7f5692b7c58ab6"
                ]
            }
        ];

        for case in cases {
            let g =
                <G1 as HashToCurve<ExpandMsgXmd<sha2ni::Sha256>>>::hash_to_curve(&case.msg, DOMAIN);
            let g_uncompressed = g.into_affine().into_uncompressed();

            assert_eq!(case.expected(), hex::encode(&g_uncompressed.0[..]));
        }
    }
}
