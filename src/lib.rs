/*! PolynomialOverP ring over finite prime field $`\mathbb{F}_p[x]`$

```
use num_traits::Zero;
use polynomial_over_finite_prime_field::PolynomialOverP;
let p = PolynomialOverP::<i32>::new(vec![3, 1, 4, 1, 5, 9, 2, 6, 5, 3], 17);
let q = PolynomialOverP::<i32>::new(vec![2, 7, 1, 8, 2, 8], 17);
let mut r = p.clone();
let d = r.division(&q);
assert!((d * q + r - p).is_zero());
```
*/
mod sealed {
    pub trait SizedExt: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    impl<T> SizedExt for T where T: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    #[cfg(not(feature = "__internal_inject_debug"))]
    pub use std::marker::Sized;
    #[cfg(feature = "__internal_inject_debug")]
    pub use SizedExt as Sized;
}
use modulo_n_tools::{add_mod, mul_mod, sub_mod};
use num_traits::{One, Zero};
use ring_algorithm::{modulo_inverse, RingNormalize};
use sealed::Sized;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Rem, Sub, SubAssign};
mod ops;

/** PolynomialOverP ring over finite prime field $`F_p[x]`$

```
use num_traits::Zero;
use polynomial_over_finite_prime_field::PolynomialOverP;
let p = PolynomialOverP::<i32>::new(vec![3, 1, 4, 1, 5], 17);
let q = PolynomialOverP::<i32>::new(vec![2, 7, 1], 17);
let mut r = p.clone();
let d = r.division(&q);
assert!((d * q + r - p).is_zero());
```
*/
#[derive(Clone, Debug)]
pub struct PolynomialOverP<T> {
    coef: Vec<T>,
    prime: T,
}

impl<T: Sized> PartialEq for PolynomialOverP<T>
where
    T: Clone + Eq + Zero,
    for<'x> &'x T: Sub<Output = T> + Rem<Output = T>,
{
    fn eq(&self, other: &Self) -> bool {
        if self.prime != other.prime {
            return false;
        }
        if self.coef.len() != other.coef.len() {
            return false;
        }
        for (v, w) in self.coef.iter().zip(other.coef.iter()) {
            if !(&(v - w) % &self.prime).is_zero() {
                return false;
            }
        }
        true
    }
}
impl<T: Sized> Eq for PolynomialOverP<T>
where
    T: Clone + Eq + Zero,
    for<'x> &'x T: Sub<Output = T> + Rem<Output = T>,
{
}

impl<T: Sized> PolynomialOverP<T> {
    fn len(&self) -> usize {
        self.coef.len()
    }
    /** degree of polynomial

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 2, 1], 5); // 3+2x+x^2
    assert_eq!(p.deg(), Some(2));
    let q = PolynomialOverP::<i32>::new(vec![0], 5); // 0
    assert_eq!(q.deg(), None);
    ```
    */
    pub fn deg(&self) -> Option<usize> {
        if self.coef.is_empty() {
            None
        } else {
            Some(self.len() - 1)
        }
    }
    /** leading coefficent

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 2, 1], 5); // 3+2x+x^2
    assert_eq!(p.lc(), Some(&1));
    let q = PolynomialOverP::<i32>::new(vec![0], 5); // 0
    assert_eq!(q.lc(), None);
    ```
    */
    pub fn lc(&self) -> Option<&T> {
        self.deg().map(|d| &self.coef[d])
    }
    /** get coefficents

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 2, 1], 5); // 3+2x+x^2
    assert_eq!(p.coefs(), vec![3, 2, 1]);
    let q = PolynomialOverP::<i32>::new(vec![0], 5); // 0
    assert_eq!(q.coefs(), Vec::<i32>::new());
    ```
    */
    pub fn coefs(self) -> Vec<T> {
        self.coef
    }
    /** get prime

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 2, 1], 97); // 3+2x+x^2
    assert_eq!(p.prime_ref(), &97);
    let q = PolynomialOverP::<i32>::new(vec![0], 5); // 0
    assert_eq!(q.prime_ref(), &5);
    ```
    */
    pub fn prime_ref(&self) -> &T {
        &self.prime
    }
    fn trim_zero(&mut self)
    where
        T: Zero,
    {
        let len = self
            .coef
            .iter()
            .rposition(|x| !x.is_zero())
            .map(|pos| pos + 1)
            .unwrap_or(0);
        self.coef.truncate(len);
    }
    fn extend(&mut self, len: usize)
    where
        T: Zero,
    {
        if self.len() < len {
            self.coef.resize_with(len, T::zero);
        }
    }
    /** construct polynomial

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 2, 1], 5);
    assert_eq!(p.to_string(), "x^2+2*x+3");
    ```
    */
    pub fn new(coef: Vec<T>, prime: T) -> Self
    where
        T: Zero,
        for<'x> &'x T: Rem<Output = T>,
    {
        let mut poly = Self {
            coef: coef.into_iter().map(|x| &x % &prime).collect(),
            prime,
        };
        poly.trim_zero();
        poly
    }
    /** construct polynomial from monomial $`cx^d`$ ($`c`$=coefficent, $`d`$=degree)

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::from_monomial(3, 2, 5);
    let q = PolynomialOverP::<i32>::new(vec![0, 0, 3], 5);
    assert_eq!(p, q);
    ```
    */
    pub fn from_monomial(coefficent: T, degree: usize, prime: T) -> Self
    where
        T: Zero,
    {
        let coef = if coefficent.is_zero() {
            Vec::new()
        } else {
            let mut coef = Vec::with_capacity(degree + 1);
            for _ in 0..degree {
                coef.push(T::zero());
            }
            coef.push(coefficent);
            coef
        };
        Self { coef, prime }
    }
    fn add_assign_ref(&mut self, other: &Self)
    where
        T: Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
        for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Rem<Output = T>,
    {
        if self.is_zero() {
            *self = other.clone();
            return;
        } else if other.is_zero() {
            return;
        }
        let len = self.len();
        self.extend(other.len());
        let prime = self.prime.clone();
        self.coef
            .iter_mut()
            .zip(other.coef.iter())
            .for_each(|(l, r)| {
                *l = add_mod::<T>(&*l, r, &prime);
            });
        if len == other.len() {
            self.trim_zero()
        }
    }
    fn neg_impl(self) -> Self
    where
        T: Clone + Neg<Output = T>,
    {
        Self {
            coef: self.coef.into_iter().map(|v| -v).collect(),
            prime: self.prime,
        }
    }
    fn neg_ref(&self) -> Self
    where
        T: Clone,
        for<'x> &'x T: Neg<Output = T>,
    {
        Self {
            coef: self.coef.iter().map(|v| -v).collect(),
            prime: self.prime.clone(),
        }
    }
    fn sub_assign_ref(&mut self, other: &Self)
    where
        T: Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
        for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Rem<Output = T>,
    {
        if self.is_zero() {
            *self = -other;
            return;
        } else if other.is_zero() {
            return;
        }
        let len = self.len();
        self.extend(other.len());
        let prime = self.prime.clone();
        self.coef
            .iter_mut()
            .zip(other.coef.iter())
            .for_each(|(l, r)| *l = sub_mod::<T>(&*l, r, &prime));
        if len == other.len() {
            self.trim_zero()
        }
    }
    fn mul_impl(&self, other: &Self) -> Self
    where
        T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
        for<'x> &'x T:
            Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
    {
        if self.is_zero() || other.is_zero() {
            let prime = if self.prime.is_zero() {
                other.prime.clone()
            } else {
                self.prime.clone()
            };
            return Self::new(Vec::new(), prime);
        }
        let mut coef = vec![T::zero(); self.len() + other.len() - 1];
        let prime = self.prime.clone();
        self.coef
            .iter()
            .enumerate()
            .for_each(|(i, c)| mul_aux::<T>(&mut coef[i..], c, &other.coef, &prime));
        Self {
            coef,
            prime: self.prime.clone(),
        }
    }
    /** Make one polynomial

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::one(5);
    assert_eq!(p, PolynomialOverP::<i32>::new(vec![1], 5));
    ```
    */
    pub fn one(prime: T) -> Self
    where
        T: One,
    {
        Self {
            coef: vec![T::one()],
            prime,
        }
    }
    /** evaluate polynomial by Horner's method

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 2, 1], 7); // 3+2x+x^2
    assert_eq!(p.eval(&1), 6);
    assert_eq!(p.eval(&2), 4);
    assert_eq!(p.eval(&3), 4);
    ```
    */
    pub fn eval<'a>(&self, x: &'a T) -> T
    where
        T: Sized + Clone + Zero,
        for<'x> &'x T: Add<Output = T> + Mul<Output = T> + Rem<Output = T>,
    {
        if self.coef.is_empty() {
            return T::zero();
        }
        let mut sum = self.lc().unwrap().clone();
        for i in (0..self.len() - 1).rev() {
            sum = &(&(&sum * x) + &self.coef[i]) % &self.prime;
        }
        sum
    }
    /** derivative

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<usize>::new(vec![1, 2, 3, 2, 1], 7); // 1+2x+3x^2+2x^3+x^4
    assert_eq!(p.derivative(), PolynomialOverP::<usize>::new(vec![2, 6, 6, 4], 7));
    let q = PolynomialOverP::<usize>::new(vec![1, 2, 3, 4, 3, 2, 1], 5); // 1+2x+3x^2+4x^3+3x^4+2x^5+x^6
    assert_eq!(q.derivative(), PolynomialOverP::<usize>::new(vec![2, 1, 2, 2, 0, 1], 5));
    ```
    */
    #[must_use]
    pub fn derivative(self) -> Self
    where
        T: Clone + Zero + for<'x> AddAssign<&'x T> + TryFrom<usize>,
        for<'x> &'x T: Mul<Output = T> + Rem<Output = T>,
        <T as TryFrom<usize>>::Error: Debug,
    {
        let Self { coef, prime } = self;
        let coef = coef
            .into_iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| mul_mod::<T>(&T::try_from(i).unwrap(), &c, &prime))
            .collect();
        PolynomialOverP::<T>::new(coef, prime)
    }
    /** make polynomial monic

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let mut p = PolynomialOverP::<i32>::new(vec![1, 2, 3], 5);
    p.monic();
    let q = PolynomialOverP::<i32>::new(vec![2, 4, 1], 5);
    assert_eq!(p, q);
    ```
    */
    pub fn monic(&mut self)
    where
        T: Clone + Eq + Zero + One + RingNormalize,
        for<'x> &'x T:
            Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Rem<Output = T>,
    {
        if let Some(lc) = self.lc() {
            let prime = self.prime.clone();
            let lc_inv = modulo_inverse::<T>(lc.clone(), prime.clone()).unwrap();
            self.coef
                .iter_mut()
                .for_each(|v| *v = mul_mod::<T>(&*v, &lc_inv, &prime));
        }
    }
    /** polynomial division

    ```
    use num_traits::Zero;
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4, 1, 5], 17);
    let q = PolynomialOverP::<i32>::new(vec![2, 7, 1], 17);
    let mut r = p.clone();
    let d = r.division(&q);
    assert!((d * q + r - p).is_zero());
    ```
    */
    #[allow(unknown_lints, clippy::return_self_not_must_use)]
    pub fn division(&mut self, other: &Self) -> Self
    where
        T: Clone
            + Ord
            + Eq
            + Zero
            + One
            + for<'x> AddAssign<&'x T>
            + for<'x> SubAssign<&'x T>
            + RingNormalize,
        for<'x> &'x T:
            Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Rem<Output = T>,
    {
        let g_deg = other.deg().expect("Division by zero");
        if self.deg() < other.deg() {
            let prime = if self.prime.is_zero() {
                other.prime.clone()
            } else {
                self.prime.clone()
            };
            return Self::new(Vec::new(), prime);
        }
        let prime = self.prime.clone();
        let lc_inv = modulo_inverse::<T>(other.lc().unwrap().clone(), prime.clone()).unwrap();
        let mut coef = vec![T::zero(); self.len() - other.len() + 1];
        while self.deg() >= other.deg() {
            let d = self.deg().unwrap() - g_deg;
            let c = mul_mod::<T>(self.lc().unwrap(), &lc_inv, &prime);
            for i in 0..other.len() - 1 {
                self.coef[i + d] = &(&self.coef[i + d] - &(&c * &other.coef[i])) % &prime;
            }
            self.coef.pop(); // new deg < prev deg
            self.trim_zero();
            coef[d] = c;
        }
        Self {
            coef,
            prime: self.prime.clone(),
        }
    }
    /** calclate $`\displaystyle \sqrt[p]{f}`$

    Return `Some` if input polynomial is $`\displaystyle \sum_{i=0}^{n} c_{ip}x^{ip}`$.
    Otherwize returns `None`.
    ```
    use num_traits::Zero;
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let p = PolynomialOverP::<usize>::new(vec![3, 0, 1, 0, 4], 2);
    let q = p.pth_root().unwrap();
    assert_eq!(q, PolynomialOverP::<usize>::new(vec![3, 1, 4], 2));
    let p = PolynomialOverP::<usize>::new(vec![2, 0, 0, 7, 0, 0, 1], 3);
    let q = p.pth_root().unwrap();
    assert_eq!(q, PolynomialOverP::<usize>::new(vec![2, 7, 1], 3));
    let p = PolynomialOverP::<usize>::new(vec![1, 0, 0, 0, 1], 2);
    let q = p.pth_root().unwrap();
    assert_eq!(q, PolynomialOverP::<usize>::new(vec![1, 0, 1], 2));
    let p = PolynomialOverP::<usize>::new(vec![3, 0, 1, 0, 2, 9], 2);
    let q = p.pth_root();
    assert_eq!(q, None);
    ```
    */
    pub fn pth_root(mut self) -> Option<Self>
    where
        T: Clone + Eq + Zero + TryFrom<usize> + TryInto<usize>,
        for<'x> &'x T: Rem<Output = T>,
        usize: std::convert::TryFrom<T>,
        <T as TryFrom<usize>>::Error: Debug,
    {
        //     i % p != 0 => c == 0
        // <=> !(i % p !=0) || c == 0
        // <=> i % p == 0 || c == 0
        let b = self
            .coef
            .iter()
            .enumerate()
            .all(|(i, c)| (&T::try_from(i).unwrap() % &self.prime).is_zero() || c.is_zero());
        if !b {
            return None;
        }
        let n = self.coef.len();
        let p = match usize::try_from(self.prime.clone()) {
            Ok(x) => x,
            Err(_) => return None,
        };
        let mut i = 1;
        let mut j = p;
        while j < n {
            let (c1, c2) = self.coef.split_at_mut(j);
            std::mem::swap(&mut c1[i], &mut c2[0]);
            i += 1;
            j += p;
        }
        self.coef.truncate(i);
        Some(Self::new(self.coef, self.prime))
    }
    /** to positive

    ```
    use polynomial_over_finite_prime_field::PolynomialOverP;
    let mut p = PolynomialOverP::<i32>::new(vec![-3, 1, -4, 1, -5], 7);
    p.to_positive();
    let q = PolynomialOverP::<i32>::new(vec![4, 1, 3, 1, 2], 7);
    assert_eq!(p, q);
    ```
    */
    pub fn to_positive(&mut self)
    where
        T: Zero + Ord + for<'x> AddAssign<&'x T>,
    {
        let z = T::zero();
        let p = &self.prime;
        self.coef.iter_mut().for_each(|x| {
            if *x < z {
                *x += p;
            }
        });
    }
}

fn mul_aux<T>(sum: &mut [T], coef: &T, vec: &[T], prime: &T)
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    sum.iter_mut()
        .zip(vec.iter())
        .for_each(|(l, r)| *l = &(&(coef * r) + &*l) % prime);
}

impl<T> std::fmt::Display for PolynomialOverP<T>
where
    T: std::cmp::Eq + std::fmt::Display + Zero + One,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let vec = &self.coef;
        if vec.is_empty() {
            return write!(f, "0");
        }
        let mut is_first = true;
        for (i, c) in vec.iter().enumerate().rev() {
            if c.is_zero() {
                continue;
            }
            if is_first {
                is_first = false;
            } else {
                write!(f, "+")?
            }
            if c.is_one() {
                match i {
                    0 => write!(f, "1")?,
                    1 => write!(f, "x")?,
                    _ => write!(f, "x^{i}")?,
                }
            } else {
                match i {
                    0 => write!(f, "{c}")?,
                    1 => write!(f, "{c}*x")?,
                    _ => write!(f, "{c}*x^{i}")?,
                }
            }
        }
        Ok(())
    }
}

impl<T> RingNormalize for PolynomialOverP<T>
where
    T: Sized + Clone + Eq + Zero + One + RingNormalize,
    for<'x> &'x T:
        Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Rem<Output = T>,
{
    fn leading_unit(&self) -> Self {
        if let Some(x) = self.lc() {
            Self::from_monomial(x.clone(), 0, self.prime.clone())
        } else {
            Self::one(self.prime.clone())
        }
    }
    fn normalize_mut(&mut self) {
        self.monic();
    }
}

impl<T> Zero for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Rem<Output = T>,
{
    fn zero() -> Self {
        Self {
            coef: Vec::new(),
            prime: T::zero(),
        }
    }
    fn is_zero(&self) -> bool {
        self.deg().is_none()
    }
}
