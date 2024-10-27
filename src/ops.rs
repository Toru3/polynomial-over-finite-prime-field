use crate::{PolynomialOverP, Sized};
use num_traits::{One, Zero};
use ring_algorithm::RingNormalize;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};
#[auto_impl_ops::auto_ops]
impl<T> AddAssign<&PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Rem<Output = T>,
{
    fn add_assign(&mut self, other: &Self) {
        self.add_assign_ref(other);
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
impl<T> Neg for &PolynomialOverP<T>
where
    T: Sized + Clone,
    for<'x> &'x T: Neg<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn neg(self) -> Self::Output {
        self.neg_ref()
    }
}

#[auto_impl_ops::auto_ops]
impl<T> SubAssign<&PolynomialOverP<T>> for PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T: Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Rem<Output = T>,
{
    fn sub_assign(&mut self, other: &Self) {
        self.sub_assign_ref(other)
    }
}

#[auto_impl_ops::auto_ops]
impl<T> Mul for &PolynomialOverP<T>
where
    T: Sized + Clone + Ord + Zero + for<'x> AddAssign<&'x T> + for<'x> SubAssign<&'x T>,
    for<'x> &'x T:
        Add<Output = T> + Sub<Output = T> + Neg<Output = T> + Mul<Output = T> + Rem<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn mul(self, other: Self) -> Self::Output {
        self.mul_impl(other)
    }
}

#[auto_impl_ops::auto_ops]
impl<T> Div for &PolynomialOverP<T>
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
    for<'x> &'x T:
        Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Rem<Output = T>,
{
    type Output = PolynomialOverP<T>;
    fn div(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f.division(other)
    }
}

#[auto_impl_ops::auto_ops]
impl<T> RemAssign<&PolynomialOverP<T>> for PolynomialOverP<T>
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
    for<'x> &'x T:
        Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Rem<Output = T>,
{
    fn rem_assign(&mut self, other: &Self) {
        self.division(other);
    }
}
