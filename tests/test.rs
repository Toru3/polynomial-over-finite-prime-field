use num_traits::Zero;
use polynomial_over_finite_prime_field::*;

#[test]
fn add() {
    let p = PolynomialOverP::<i32>::new(vec![1, 4, 1, 4, 1, 3, 5, 6], 7);
    let q = PolynomialOverP::<i32>::new(vec![2, 2, 3, 6, 0, 6, 7, 9], 7);
    assert_eq!(
        p + q,
        PolynomialOverP::<i32>::new(vec![3, 6, 4, 3, 1, 2, 5, 1], 7)
    );
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4], 11);
    let q = PolynomialOverP::<i32>::new(vec![-3, -1, -4], 11);
    assert!((p + q).is_zero());
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4], 17);
    let q = PolynomialOverP::<i32>::new(vec![1, 2, -4], 17);
    assert_eq!(p + q, PolynomialOverP::<i32>::new(vec![4, 3], 17));
    let p = PolynomialOverP::<i32>::new(vec![-3, 1, 4], 7);
    let q = PolynomialOverP::<i32>::new(vec![3, -1, 1], 7);
    assert_eq!(p + q, PolynomialOverP::<i32>::new(vec![0, 0, 5], 7));
    let p = PolynomialOverP::<i32>::new(vec![-3, 1, 4], 5);
    let q = PolynomialOverP::<i32>::new(vec![3, -1, 1], 5);
    assert_eq!(p + q, PolynomialOverP::<i32>::new(vec![], 5));
}

#[test]
fn neg() {
    let p = PolynomialOverP::<i32>::new(vec![1, 4, 1, 4, 1, 3, 5, 6], 7);
    assert_eq!(
        -p,
        PolynomialOverP::<i32>::new(vec![-1, -4, -1, -4, -1, -3, -5, -6], 7)
    );
}

#[test]
fn sub() {
    let p = PolynomialOverP::<i32>::new(vec![1, 4, 1, 4, 1, 3, 5, 6], 7);
    let q = PolynomialOverP::<i32>::new(vec![2, 2, 3, 6, 0, 6, 7, 9], 7);
    let r = PolynomialOverP::<i32>::new(vec![-1, 2, -2, -2, 1, -3, 5, 4], 7);
    assert_eq!(p - q, r);
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4], 11);
    let q = p.clone();
    assert!((p - q).is_zero());
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4], 17);
    let q = PolynomialOverP::<i32>::new(vec![1, 2, 4], 17);
    assert_eq!(p - q, PolynomialOverP::<i32>::new(vec![2, -1], 17));
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4], 7);
    let q = PolynomialOverP::<i32>::new(vec![3, 1, -1], 7);
    assert_eq!(p - q, PolynomialOverP::<i32>::new(vec![0, 0, 5], 7));
    let p = PolynomialOverP::<i32>::new(vec![3, 1, 4], 5);
    let q = PolynomialOverP::<i32>::new(vec![3, 1, -1], 5);
    assert_eq!(p - q, PolynomialOverP::<i32>::new(vec![], 5));
}

#[test]
fn mul() {
    let p = PolynomialOverP::<i32>::new(vec![1, 2, 3], 5);
    let q = PolynomialOverP::<i32>::new(vec![4, 5, 6], 5);
    let r = PolynomialOverP::<i32>::new(vec![4, 3, 3, 2, 3], 5);
    assert_eq!(p * q, r);
}

#[test]
fn eq() {
    let p = PolynomialOverP::<i32>::new(vec![1, 2, 3], 5);
    let q = PolynomialOverP::<i32>::new(vec![-4, -3, -2], 5);
    assert_eq!(p, q);
    let p = PolynomialOverP::<i32>::new(vec![1, 2, 3], 5);
    let q = PolynomialOverP::<i32>::new(vec![1, 2, 3], 7);
    assert_ne!(p, q);
}
