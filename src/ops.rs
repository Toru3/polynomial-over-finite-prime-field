use crate::{PolynomialOverP, Sized};
use num_traits::{One, Zero};
use ring_algorithm::RingNormalize;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};
// AddAssign
impl<'a, T> AddAssign<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    fn add_assign(&mut self, other: &Self) {
        self.add_assign_ref(other);
    }
}
impl<T> AddAssign for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    fn add_assign(&mut self, other: Self) {
        *self += &other
    }
}

// Add
impl<'a, T> Add for &'a PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn add(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f += other;
        f
    }
}
impl<'a, T> Add<PolynomialOverP<T>> for &'a PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn add(self, other: PolynomialOverP<T>) -> Self::Output {
        let mut f = self.clone();
        f += &other;
        f
    }
}
impl<'a, T> Add<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    type Output = Self;
    fn add(mut self, other: &Self) -> Self::Output {
        self += other;
        self
    }
}
impl<T> Add for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T>,
{
    type Output = Self;
    fn add(mut self, other: PolynomialOverP<T>) -> Self::Output {
        self += &other;
        self
    }
}

// Neg
impl<T> Neg for PolynomialOverP<T>
where
    T: Sized + Clone + Neg<Output = T>,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.neg_impl()
    }
}
impl<'a, T> Neg for &'a PolynomialOverP<T>
where
    T: Sized + Clone,
    for<'x> &'x T: Neg<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn neg(self) -> Self::Output {
        self.neg_ref()
    }
}

// SubAssign
impl<'a, T> SubAssign<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T>,
{
    fn sub_assign(&mut self, other: &Self) {
        self.sub_assign_ref(other)
    }
}
impl<T> SubAssign for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T>,
{
    fn sub_assign(&mut self, other: Self) {
        *self -= &other
    }
}

// Sub
impl<'a, T> Sub for &'a PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn sub(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f -= other;
        f
    }
}
impl<'a, T> Sub<PolynomialOverP<T>> for &'a PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn sub(self, other: PolynomialOverP<T>) -> Self::Output {
        let mut f = self.clone();
        f -= &other;
        f
    }
}
impl<'a, T> Sub<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T>,
{
    type Output = Self;
    fn sub(mut self, other: &Self) -> Self::Output {
        self -= other;
        self
    }
}
impl<T> Sub for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T>,
{
    type Output = Self;
    fn sub(mut self, other: PolynomialOverP<T>) -> Self::Output {
        self -= &other;
        self
    }
}

// Mul
impl<'a, T> Mul for &'a PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn mul(self, other: Self) -> Self::Output {
        self.mul_impl(other)
    }
}
impl<'a, T> Mul<PolynomialOverP<T>> for &'a PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn mul(self, other: PolynomialOverP<T>) -> Self::Output {
        self * &other
    }
}
impl<T> Mul for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}
impl<'a, T> Mul<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    type Output = Self;
    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

// MulAssign
impl<'a, T> MulAssign<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    fn mul_assign(&mut self, other: &Self) {
        *self = &*self * other;
    }
}
impl<T> MulAssign for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    fn mul_assign(&mut self, other: Self) {
        *self = &*self * &other;
    }
}

// Div
impl<'a, T> Div for &'a PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = PolynomialOverP<T>;
    fn div(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f.division(other)
    }
}
impl<'a, T> Div<PolynomialOverP<T>> for &'a PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = PolynomialOverP<T>;
    fn div(self, other: PolynomialOverP<T>) -> Self::Output {
        let mut f = self.clone();
        f.division(&other)
    }
}
impl<T> Div for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = Self;
    fn div(mut self, other: Self) -> Self::Output {
        self.division(&other)
    }
}
impl<'a, T> Div<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = Self;
    fn div(mut self, other: &Self) -> Self::Output {
        self.division(other)
    }
}

// DivAssign
impl<'a, T> DivAssign<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    fn div_assign(&mut self, other: &Self) {
        *self = &*self / other;
    }
}
impl<T> DivAssign for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    fn div_assign(&mut self, other: Self) {
        *self = &*self / &other;
    }
}

// RemAssign
impl<'a, T> RemAssign<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    fn rem_assign(&mut self, other: &Self) {
        self.division(other);
    }
}
impl<'a, T> RemAssign for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    fn rem_assign(&mut self, other: Self) {
        self.division(&other);
    }
}

// Rem
impl<'a, T> Rem for &'a PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = PolynomialOverP<T>;
    fn rem(self, other: Self) -> Self::Output {
        let mut t = self.clone();
        t %= other;
        t
    }
}
impl<'a, T> Rem<PolynomialOverP<T>> for &'a PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = PolynomialOverP<T>;
    fn rem(self, other: PolynomialOverP<T>) -> Self::Output {
        let mut t = self.clone();
        t %= other;
        t
    }
}
impl<'a, T> Rem<&'a PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = Self;
    fn rem(mut self, other: &Self) -> Self::Output {
        self %= other;
        self
    }
}
impl<T> Rem for PolynomialOverP<T>
where
    T: Sized
        + Clone
        + Ord
        + Eq
        + Zero
        + One
        + for<'x> AddAssign<&'x T>
        + for<'x> SubAssign<&'x T>
        + RingNormalize,
    for<'x> &'x T: Add<Output=T> + Sub<Output=T> + Mul<Output=T> + Div<Output=T> + Rem<Output=T>,
{
    type Output = Self;
    fn rem(mut self, other: Self) -> Self::Output {
        self %= &other;
        self
    }
}
