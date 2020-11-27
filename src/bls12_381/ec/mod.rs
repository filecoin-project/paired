macro_rules! curve_impl {
    (
        $name:expr,
        $projective:ident,
        $affine:ident,
        $prepared:ident,
        $basefield:ident,
        $scalarfield:ident,
        $uncompressed:ident,
        $compressed:ident,
        $pairing:ident,
        $iso_1:expr,
        $iso_2:expr
    ) => {
        #[derive(Copy, Clone, PartialEq, Eq, Debug)]
        pub struct $affine {
            pub(crate) x: $basefield,
            pub(crate) y: $basefield,
            pub(crate) infinity: bool,
        }

        impl ::std::fmt::Display for $affine {
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                if self.infinity {
                    write!(f, "{}(Infinity)", $name)
                } else {
                    write!(f, "{}(x={}, y={})", $name, self.x, self.y)
                }
            }
        }

        fn y2_from_x(x: $basefield) -> $basefield {
            let mut y2 = x.clone();
            y2.square();
            y2.mul_assign(&x);
            y2.add_assign(&$affine::get_coeff_b());
            y2
        }

        #[derive(Copy, Clone, Debug, Eq)]
        pub struct $projective {
            pub(crate) x: $basefield,
            pub(crate) y: $basefield,
            pub(crate) z: $basefield,
        }

        impl ::std::fmt::Display for $projective {
            fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
                write!(f, "{}", self.into_affine())
            }
        }

        impl PartialEq for $projective {
            fn eq(&self, other: &$projective) -> bool {
                if self.is_zero() {
                    return other.is_zero();
                }
                if other.is_zero() {
                    return false;
                }

                if self.is_normalized() {
                    if other.is_normalized() {
                        return self.into_affine() == other.into_affine();
                    }
                }

                // The points (X, Y, Z) and (X', Y', Z')
                // are equal when (X * Z^2) = (X' * Z'^2)
                // and (Y * Z^3) = (Y' * Z'^3).

                let mut z1 = self.z;
                z1.square();
                let mut z2 = other.z;
                z2.square();

                let mut tmp1 = self.x;
                tmp1.mul_assign(&z2);

                let mut tmp2 = other.x;
                tmp2.mul_assign(&z1);

                if tmp1 != tmp2 {
                    return false;
                }

                z1.mul_assign(&self.z);
                z2.mul_assign(&other.z);
                z2.mul_assign(&self.y);
                z1.mul_assign(&other.y);

                if z1 != z2 {
                    return false;
                }

                true
            }
        }

        impl $projective {
            /// Generic isogeny evaluation function.
            pub(crate) fn eval_iso(&mut self, coeffs: [&[$basefield]; 4]) {
                // Rust (still) can't handle generic array sizes (issue #43408)
                let mut tmp = [$basefield::zero(); $iso_1];
                let mut mapvals = [$basefield::zero(); 4];
                // scope for pt borrow
                {
                    // unpack input point
                    let x = &self.x;
                    let y = &self.y;
                    let z = &self.z;

                    // precompute powers of z
                    let zpows = {
                        let mut zpows = [$basefield::zero(); $iso_2];
                        zpows[0] = *z;
                        zpows[0].square(); // z^2
                        zpows[1] = zpows[0];
                        zpows[1].square(); // z^4
                        {
                            let (z_squared, rest) = zpows.split_at_mut(1);
                            for idx in 1..coeffs[2].len() - 2 {
                                if idx % 2 == 0 {
                                    rest[idx] = rest[idx / 2 - 1];
                                    rest[idx].square();
                                } else {
                                    rest[idx] = rest[idx - 1];
                                    rest[idx].mul_assign(&z_squared[0]);
                                }
                            }
                        }
                        zpows
                    };

                    for idx in 0..4 {
                        let clen = coeffs[idx].len() - 1;
                        // multiply coeffs by powers of Z
                        for jdx in 0..clen {
                            tmp[jdx] = coeffs[idx][clen - 1 - jdx];
                            tmp[jdx].mul_assign(&zpows[jdx]);
                        }
                        // compute map value by Horner's rule
                        mapvals[idx] = coeffs[idx][clen];
                        for tmpval in &tmp[..clen] {
                            mapvals[idx].mul_assign(x);
                            mapvals[idx].add_assign(tmpval);
                        }
                    }

                    // x denominator is order 1 less than x numerator, so we need an extra factor of Z^2
                    mapvals[1].mul_assign(&zpows[0]);

                    // multiply result of Y map by the y-coord, y / z^3
                    mapvals[2].mul_assign(y);
                    mapvals[3].mul_assign(z);
                    mapvals[3].mul_assign(&zpows[0]);
                } // pt is no longer borrowed here

                // hack to simultaneously access elements of tmp
                let (xx, yy, zz) = {
                    let (xx, rest) = tmp.split_at_mut(1);
                    let (yy, rest) = rest.split_at_mut(1);
                    (&mut xx[0], &mut yy[0], &mut rest[0])
                };

                // compute Jacobian coordinates of resulting point
                *zz = mapvals[1];
                zz.mul_assign(&mapvals[3]); // Zout = xden * yden

                *xx = mapvals[0];
                xx.mul_assign(&mapvals[3]); // xnum * yden
                xx.mul_assign(zz); // xnum * xden * yden^2

                *yy = *zz;
                yy.square(); // xden^2 * yden^2
                yy.mul_assign(&mapvals[2]); // ynum * xden^2 * yden^2
                yy.mul_assign(&mapvals[1]); // ynum * xden^3 * yden^2

                self.x = *xx;
                self.y = *yy;
                self.z = *zz;
            }
        }

        impl $affine {
            fn mul_bits<S: AsRef<[u64]>>(&self, bits: BitIterator<S>) -> $projective {
                let mut res = $projective::zero();
                for i in bits {
                    res.double();
                    if i {
                        res.add_assign_mixed(self)
                    }
                }
                res
            }

            /// Attempts to construct an affine point given an x-coordinate. The
            /// point is not guaranteed to be in the prime order subgroup.
            ///
            /// If and only if `greatest` is set will the lexicographically
            /// largest y-coordinate be selected.
            fn get_point_from_x(x: $basefield, greatest: bool) -> Option<$affine> {
                // Compute x^3 + b
                let mut x3b = x;
                x3b.square();
                x3b.mul_assign(&x);
                x3b.add_assign(&$affine::get_coeff_b());

                x3b.sqrt().map(|y| {
                    let mut negy = y;
                    negy.negate();

                    $affine {
                        x,
                        y: if (y < negy) ^ greatest { y } else { negy },
                        infinity: false,
                    }
                })
            }

            fn is_on_curve(&self) -> bool {
                if self.is_zero() {
                    true
                } else {
                    // Check that the point is on the curve
                    let mut y2 = self.y;
                    y2.square();

                    y2 == y2_from_x(self.x)
                }
            }

            fn is_in_correct_subgroup_assuming_on_curve(&self) -> bool {
                self.mul($scalarfield::char()).is_zero()
            }

            /// Implements the Shallue–van de Woestijne encoding described in
            /// Section 3, "Indifferentiable Hashing to Barreto–Naehrig Curves"
            /// from Foque-Tibouchi: <https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf>.
            ///
            /// The encoding is adapted for BLS12-381.
            ///
            /// This encoding produces a point in E/E'. It does not reach every
            /// point. The resulting point may not be in the prime order subgroup,
            /// but it will be on the curve. It could be the point at infinity.
            ///
            /// ## Description
            ///
            /// Lemma 3 gives us three points:
            ///
            /// x_1 = (-1 + sqrt(-3))/2 - (sqrt(-3) * t^2)/(1 + b + t^2)
            /// x_2 = (-1 - sqrt(-3))/2 + (sqrt(-3) * t^2)/(1 + b + t^2)
            /// x_3 = 1 - (1 + b + t^2)^2/(3 * t^2)
            ///
            /// Given t != 0 and t != 1 + b + t^2 != 0, at least one of
            /// these three points (x1, x2, x3) is valid on the curve.
            ///
            /// In the paper, 1 + b + t^2 != 0 has no solutions, but for
            /// E(Fq) in our construction, it does have two solutions.
            /// We follow the convention of the paper by mapping these
            /// to some arbitrary points; in our case, the positive/negative
            /// fixed generator (with the parity of the y-coordinate
            /// corresponding to the t value).
            ///
            /// Unlike the paper, which maps t = 0 to an arbitrary point,
            /// we map it to the point at infinity. This arrangement allows
            /// us to preserve sw_encode(t) = sw_encode(-t) for all t.
            ///
            /// We choose the smallest i such that x_i is on the curve.
            /// We choose the corresponding y-coordinate with the same
            /// parity, defined as the point being lexicographically larger
            /// than its negative.
            fn sw_encode(t: $basefield) -> Self {
                // Handle the case t == 0
                if t.is_zero() {
                    return Self::zero();
                }

                // We choose the corresponding y-coordinate with the same parity as t.
                let parity = t.parity();

                // w = (t^2 + b + 1)^(-1) * sqrt(-3) * t
                let mut w = t;
                w.square();
                w.add_assign(&$affine::get_coeff_b());
                w.add_assign(&$basefield::one());

                // Handle the case t^2 + b + 1 == 0
                if w.is_zero() {
                    let mut ret = Self::one();
                    if parity {
                        ret.negate()
                    }
                    return ret;
                }

                w = w.inverse().unwrap();
                w.mul_assign(&$basefield::get_swenc_sqrt_neg_three());
                w.mul_assign(&t);

                // x1 = - wt  + (sqrt(-3) - 1) / 2
                let mut x1 = w;
                x1.mul_assign(&t);
                x1.negate();
                x1.add_assign(&$basefield::get_swenc_sqrt_neg_three_minus_one_div_two());
                if let Some(p) = Self::get_point_from_x(x1, parity) {
                    return p;
                }

                // x2 = -1 - x1
                let mut x2 = x1;
                x2.negate();
                x2.sub_assign(&$basefield::one());
                if let Some(p) = Self::get_point_from_x(x2, parity) {
                    return p;
                }

                // x3 = 1/w^2 + 1
                let mut x3 = w;
                x3.square();
                x3 = x3.inverse().unwrap();
                x3.add_assign(&$basefield::one());
                Self::get_point_from_x(x3, parity)
                    .expect("this point must be valid if the other two are not")
            }
        }

        impl CurveAffine for $affine {
            type Engine = Bls12;
            type Scalar = $scalarfield;
            type Base = $basefield;
            type Projective = $projective;
            type Uncompressed = $uncompressed;
            type Compressed = $compressed;

            fn zero() -> Self {
                $affine {
                    x: $basefield::zero(),
                    y: $basefield::one(),
                    infinity: true,
                }
            }

            fn one() -> Self {
                Self::get_generator()
            }

            fn is_zero(&self) -> bool {
                self.infinity
            }

            fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, by: S) -> $projective {
                let bits = BitIterator::new(by.into());
                self.mul_bits(bits)
            }

            fn negate(&mut self) {
                if !self.is_zero() {
                    self.y.negate();
                }
            }

            fn into_projective(&self) -> $projective {
                (*self).into()
            }

            fn as_bytes(&self) -> Vec<u8> {
                let mut out = Vec::new();
                out.extend_from_slice(&self.x.as_bytes());
                out.extend_from_slice(&self.y.as_bytes());
                out.push(self.infinity as u8);
                out
            }
        }

        impl PairingCurveAffine for $affine {
            type Prepared = $prepared;
            type Pair = $pairing;
            type PairingResult = Fq12;

            fn prepare(&self) -> Self::Prepared {
                $prepared::from_affine(*self)
            }

            fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult {
                self.perform_pairing(other)
            }
        }

        impl CurveProjective for $projective {
            type Engine = Bls12;
            type Scalar = $scalarfield;
            type Base = $basefield;
            type Affine = $affine;

            fn random<R: RngCore>(rng: &mut R) -> Self {
                loop {
                    let x = $basefield::random(rng);
                    let greatest = rng.next_u32() % 2 != 0;

                    if let Some(p) = $affine::get_point_from_x(x, greatest) {
                        let p = p.scale_by_cofactor();

                        if !p.is_zero() {
                            return p;
                        }
                    }
                }
            }

            // The point at infinity is always represented by
            // Z = 0.
            fn zero() -> Self {
                $projective {
                    x: $basefield::zero(),
                    y: $basefield::one(),
                    z: $basefield::zero(),
                }
            }

            fn one() -> Self {
                $affine::one().into()
            }

            // The point at infinity is always represented by
            // Z = 0.
            fn is_zero(&self) -> bool {
                self.z.is_zero()
            }

            fn is_normalized(&self) -> bool {
                self.is_zero() || self.z == $basefield::one()
            }

            fn batch_normalization<S: std::borrow::BorrowMut<Self>>(v: &mut [S]) {
                // Montgomery’s Trick and Fast Implementation of Masked AES
                // Genelle, Prouff and Quisquater
                // Section 3.2

                // First pass: compute [a, ab, abc, ...]
                let mut prod = Vec::with_capacity(v.len());
                let mut tmp = $basefield::one();
                for g in v
                    .iter_mut()
                    .map(|g| g.borrow_mut())
                    // Ignore normalized elements
                    .filter(|g| !g.is_normalized())
                {
                    tmp.mul_assign(&g.z);
                    prod.push(tmp);
                }

                // Invert `tmp`.
                tmp = tmp.inverse().unwrap(); // Guaranteed to be nonzero.

                // Second pass: iterate backwards to compute inverses
                for (g, s) in v
                    .iter_mut()
                    .map(|g| g.borrow_mut())
                    // Backwards
                    .rev()
                    // Ignore normalized elements
                    .filter(|g| !g.is_normalized())
                    // Backwards, skip last element, fill in one for last term.
                    .zip(
                        prod.into_iter()
                            .rev()
                            .skip(1)
                            .chain(Some($basefield::one())),
                    )
                {
                    // tmp := tmp * g.z; g.z := tmp * s = 1/z
                    let mut newtmp = tmp;
                    newtmp.mul_assign(&g.z);
                    g.z = tmp;
                    g.z.mul_assign(&s);
                    tmp = newtmp;
                }

                // Perform affine transformations
                for g in v
                    .iter_mut()
                    .map(|g| g.borrow_mut())
                    .filter(|g| !g.is_normalized())
                {
                    let mut z = g.z; // 1/z
                    z.square(); // 1/z^2
                    g.x.mul_assign(&z); // x/z^2
                    z.mul_assign(&g.z); // 1/z^3
                    g.y.mul_assign(&z); // y/z^3
                    g.z = $basefield::one(); // z = 1
                }
            }

            fn double(&mut self) {
                if self.is_zero() {
                    return;
                }

                // Other than the point at infinity, no points on E or E'
                // can double to equal the point at infinity, as y=0 is
                // never true for points on the curve. (-4 and -4u-4
                // are not cubic residue in their respective fields.)

                // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

                // A = X1^2
                let mut a = self.x;
                a.square();

                // B = Y1^2
                let mut b = self.y;
                b.square();

                // C = B^2
                let mut c = b;
                c.square();

                // D = 2*((X1+B)2-A-C)
                let mut d = self.x;
                d.add_assign(&b);
                d.square();
                d.sub_assign(&a);
                d.sub_assign(&c);
                d.double();

                // E = 3*A
                let mut e = a;
                e.double();
                e.add_assign(&a);

                // F = E^2
                let mut f = e;
                f.square();

                // Z3 = 2*Y1*Z1
                self.z.mul_assign(&self.y);
                self.z.double();

                // X3 = F-2*D
                self.x = f;
                self.x.sub_assign(&d);
                self.x.sub_assign(&d);

                // Y3 = E*(D-X3)-8*C
                self.y = d;
                self.y.sub_assign(&self.x);
                self.y.mul_assign(&e);
                c.double();
                c.double();
                c.double();
                self.y.sub_assign(&c);
            }

            fn add_assign(&mut self, other: &Self) {
                if self.is_zero() {
                    *self = *other;
                    return;
                }

                if other.is_zero() {
                    return;
                }

                // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl

                // Z1Z1 = Z1^2
                let mut z1z1 = self.z;
                z1z1.square();

                // Z2Z2 = Z2^2
                let mut z2z2 = other.z;
                z2z2.square();

                // U1 = X1*Z2Z2
                let mut u1 = self.x;
                u1.mul_assign(&z2z2);

                // U2 = X2*Z1Z1
                let mut u2 = other.x;
                u2.mul_assign(&z1z1);

                // S1 = Y1*Z2*Z2Z2
                let mut s1 = self.y;
                s1.mul_assign(&other.z);
                s1.mul_assign(&z2z2);

                // S2 = Y2*Z1*Z1Z1
                let mut s2 = other.y;
                s2.mul_assign(&self.z);
                s2.mul_assign(&z1z1);

                if u1 == u2 && s1 == s2 {
                    // The two points are equal, so we double.
                    self.double();
                } else {
                    // If we're adding -a and a together, self.z becomes zero as H becomes zero.

                    // H = U2-U1
                    let mut h = u2;
                    h.sub_assign(&u1);

                    // I = (2*H)^2
                    let mut i = h;
                    i.double();
                    i.square();

                    // J = H*I
                    let mut j = h;
                    j.mul_assign(&i);

                    // r = 2*(S2-S1)
                    let mut r = s2;
                    r.sub_assign(&s1);
                    r.double();

                    // V = U1*I
                    let mut v = u1;
                    v.mul_assign(&i);

                    // X3 = r^2 - J - 2*V
                    self.x = r;
                    self.x.square();
                    self.x.sub_assign(&j);
                    self.x.sub_assign(&v);
                    self.x.sub_assign(&v);

                    // Y3 = r*(V - X3) - 2*S1*J
                    self.y = v;
                    self.y.sub_assign(&self.x);
                    self.y.mul_assign(&r);
                    s1.mul_assign(&j); // S1 = S1 * J * 2
                    s1.double();
                    self.y.sub_assign(&s1);

                    // Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
                    self.z.add_assign(&other.z);
                    self.z.square();
                    self.z.sub_assign(&z1z1);
                    self.z.sub_assign(&z2z2);
                    self.z.mul_assign(&h);
                }
            }

            fn add_assign_mixed(&mut self, other: &Self::Affine) {
                if other.is_zero() {
                    return;
                }

                if self.is_zero() {
                    self.x = other.x;
                    self.y = other.y;
                    self.z = $basefield::one();
                    return;
                }

                // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl

                // Z1Z1 = Z1^2
                let mut z1z1 = self.z;
                z1z1.square();

                // U2 = X2*Z1Z1
                let mut u2 = other.x;
                u2.mul_assign(&z1z1);

                // S2 = Y2*Z1*Z1Z1
                let mut s2 = other.y;
                s2.mul_assign(&self.z);
                s2.mul_assign(&z1z1);

                if self.x == u2 && self.y == s2 {
                    // The two points are equal, so we double.
                    self.double();
                } else {
                    // If we're adding -a and a together, self.z becomes zero as H becomes zero.

                    // H = U2-X1
                    let mut h = u2;
                    h.sub_assign(&self.x);

                    // HH = H^2
                    let mut hh = h;
                    hh.square();

                    // I = 4*HH
                    let mut i = hh;
                    i.double();
                    i.double();

                    // J = H*I
                    let mut j = h;
                    j.mul_assign(&i);

                    // r = 2*(S2-Y1)
                    let mut r = s2;
                    r.sub_assign(&self.y);
                    r.double();

                    // V = X1*I
                    let mut v = self.x;
                    v.mul_assign(&i);

                    // X3 = r^2 - J - 2*V
                    self.x = r;
                    self.x.square();
                    self.x.sub_assign(&j);
                    self.x.sub_assign(&v);
                    self.x.sub_assign(&v);

                    // Y3 = r*(V-X3)-2*Y1*J
                    j.mul_assign(&self.y); // J = 2*Y1*J
                    j.double();
                    self.y = v;
                    self.y.sub_assign(&self.x);
                    self.y.mul_assign(&r);
                    self.y.sub_assign(&j);

                    // Z3 = (Z1+H)^2-Z1Z1-HH
                    self.z.add_assign(&h);
                    self.z.square();
                    self.z.sub_assign(&z1z1);
                    self.z.sub_assign(&hh);
                }
            }

            fn negate(&mut self) {
                if !self.is_zero() {
                    self.y.negate()
                }
            }

            fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S) {
                let mut res = Self::zero();

                let mut found_one = false;

                for i in BitIterator::new(other.into()) {
                    if found_one {
                        res.double();
                    } else {
                        found_one = i;
                    }

                    if i {
                        res.add_assign(self);
                    }
                }

                *self = res;
            }

            fn into_affine(&self) -> $affine {
                (*self).into()
            }

            fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize {
                Self::empirical_recommended_wnaf_for_scalar(scalar)
            }

            fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
                Self::empirical_recommended_wnaf_for_num_scalars(num_scalars)
            }

            /// Implements "Indifferentiable Hashing to Barreto–Naehrig Curves" from Foque-Tibouchi.
            /// <https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf>
            fn hash(msg: &[u8]) -> Self {
                // The construction of Foque et al. requires us to construct two
                // "random oracles" in the field, encode their image with `sw_encode`,
                // and finally add them.
                // We construct them appending to the message the string
                // $name_$oracle
                // For instance, the first oracle in group G1 appends: "G1_0".
                let mut hasher_0 = blake2b_simd::State::new();
                hasher_0.update(msg);
                #[allow(clippy::string_lit_as_bytes)]
                hasher_0.update($name.as_bytes());
                let mut hasher_1 = hasher_0.clone();

                hasher_0.update(b"_0");
                let t0 = Self::Base::hash(hasher_0);
                let t0 = Self::Affine::sw_encode(t0);

                hasher_1.update(b"_1");
                let t1 = Self::Base::hash(hasher_1);
                let t1 = Self::Affine::sw_encode(t1);

                let mut res = t0.into_projective();
                res.add_assign_mixed(&t1);
                res.into_affine().scale_by_cofactor()
            }
            fn as_bytes(&self) -> Vec<u8> {
                let mut out = Vec::new();
                out.extend_from_slice(&self.x.as_bytes());
                out.extend_from_slice(&self.y.as_bytes());
                out.extend_from_slice(&self.z.as_bytes());
                out
            }
        }

        // The affine point X, Y is represented in the jacobian
        // coordinates with Z = 1.
        impl From<$affine> for $projective {
            fn from(p: $affine) -> $projective {
                if p.is_zero() {
                    $projective::zero()
                } else {
                    $projective {
                        x: p.x,
                        y: p.y,
                        z: $basefield::one(),
                    }
                }
            }
        }

        // The projective point X, Y, Z is represented in the affine
        // coordinates as X/Z^2, Y/Z^3.
        impl From<$projective> for $affine {
            fn from(p: $projective) -> $affine {
                if p.is_zero() {
                    $affine::zero()
                } else if p.z == $basefield::one() {
                    // If Z is one, the point is already normalized.
                    $affine {
                        x: p.x,
                        y: p.y,
                        infinity: false,
                    }
                } else {
                    // Z is nonzero, so it must have an inverse in a field.
                    let zinv = p.z.inverse().unwrap();
                    let mut zinv_powered = zinv;
                    zinv_powered.square();

                    // X/Z^2
                    let mut x = p.x;
                    x.mul_assign(&zinv_powered);

                    // Y/Z^3
                    let mut y = p.y;
                    zinv_powered.mul_assign(&zinv);
                    y.mul_assign(&zinv_powered);

                    $affine {
                        x,
                        y,
                        infinity: false,
                    }
                }
            }
        }

        #[cfg(test)]
        use rand_core::SeedableRng;
        #[cfg(test)]
        use rand_xorshift::XorShiftRng;

        #[test]
        fn test_hash() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
                0xbc, 0xe5,
            ]);

            let mut seed: [u8; 32] = [0u8; 32];
            for _ in 0..100 {
                rng.fill_bytes(&mut seed);
                let p = $projective::hash(&seed).into_affine();
                assert!(!p.is_zero());
                assert!(p.is_on_curve());
                assert!(p.is_in_correct_subgroup_assuming_on_curve());
            }
        }

        #[test]
        fn test_sw_encode() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
                0xbc, 0xe5,
            ]);

            for _ in 0..100 {
                let mut t = $basefield::random(&mut rng);
                let p = $affine::sw_encode(t);
                assert!(p.is_on_curve());
                assert!(!p.is_zero());

                t.negate();
                let mut minus_p = $affine::sw_encode(t).into_projective();
                minus_p.add_assign_mixed(&p);
                assert!(minus_p.is_zero());
            }
        }
    };
}

macro_rules! encoded_point_delegations {
    ($t:ident) => {
        impl AsRef<[u8]> for $t {
            fn as_ref(&self) -> &[u8] {
                &self.0
            }
        }
        impl AsMut<[u8]> for $t {
            fn as_mut(&mut self) -> &mut [u8] {
                &mut self.0
            }
        }

        impl PartialEq for $t {
            fn eq(&self, other: &$t) -> bool {
                PartialEq::eq(&self.0[..], &other.0[..])
            }
        }
        impl Eq for $t {}
        impl PartialOrd for $t {
            fn partial_cmp(&self, other: &$t) -> Option<::std::cmp::Ordering> {
                PartialOrd::partial_cmp(&self.0[..], &other.0[..])
            }
        }
        impl Ord for $t {
            fn cmp(&self, other: &Self) -> ::std::cmp::Ordering {
                Ord::cmp(&self.0[..], &other.0[..])
            }
        }

        impl ::std::hash::Hash for $t {
            fn hash<H: ::std::hash::Hasher>(&self, state: &mut H) {
                self.0[..].hash(state);
            }
        }
    };
} // encoded_point_delegations

mod chain;
mod g1;
mod g2;
mod util;

pub use self::g1::*;
pub use self::g2::*;
