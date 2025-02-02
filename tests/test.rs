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

#[test]
fn square_free_decomposition_1() {
    let p = PolynomialOverP::<i32>::new(vec![1, 3, 4, 3, 1], 5);
    let v = p.square_free_decomposition();
    let q0 = PolynomialOverP::<i32>::new(vec![1, 1], 5);
    let q1 = PolynomialOverP::<i32>::new(vec![1, 1, 1], 5);
    let w = vec![(q1, 1), (q0, 2)];
    assert_eq!(v, w);
}

#[test]
fn square_free_decomposition_2() {
    let p = PolynomialOverP::<i32>::new(vec![2, 1, 0, 3, 1, 0, 0, 4, 1, 1], 5);
    let v = p.square_free_decomposition();
    let q0 = PolynomialOverP::<i32>::new(vec![1, 1], 5);
    let q1 = PolynomialOverP::<i32>::new(vec![1, 1, 1], 5);
    let q2 = PolynomialOverP::<i32>::new(vec![2, 1, 1], 5);
    let w = vec![(q2, 1), (q1, 2), (q0, 3)];
    assert_eq!(v, w);
}

#[test]
fn square_free_decomposition_3() {
    let p = PolynomialOverP::<i32>::new(vec![2, 0, 4, 4, 2, 0, 1, 4, 4, 2, 3, 1], 5);
    let v = p.square_free_decomposition();
    let q0 = PolynomialOverP::<i32>::new(vec![1, 1], 5);
    let q1 = PolynomialOverP::<i32>::new(vec![1, 1, 1], 5);
    let q2 = PolynomialOverP::<i32>::new(vec![2, 1, 1], 5);
    let w = vec![(q2, 1), (q1, 2), (q0, 5)];
    assert_eq!(v, w);
}

#[test]
fn square_free_decomposition_4() {
    let p = PolynomialOverP::<i32>::new(vec![0, 1, 0, 0, 0, 1], 2);
    let v = p.square_free_decomposition();
    let q0 = PolynomialOverP::<i32>::new(vec![0, 1], 2);
    let q1 = PolynomialOverP::<i32>::new(vec![1, 1], 2);
    let w = vec![(q0, 1), (q1, 4)];
    assert_eq!(v, w);
}

#[test]
fn distinct_degree_factorization_1() {
    let p = PolynomialOverP::<i32>::new(vec![1, 2, 1, 1, 1], 3);
    let v = p.distinct_degree_factorization();
    let q1 = PolynomialOverP::<i32>::new(vec![2, 0, 1], 3);
    let q2 = PolynomialOverP::<i32>::new(vec![2, 1, 1], 3);
    let w = vec![Some(q1), Some(q2)];
    assert_eq!(v, w);
}
