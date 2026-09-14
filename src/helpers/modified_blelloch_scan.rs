use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::prelude::*;

pub fn prefix_sum<T, F, G>(values: &mut [T], op: F, same_subarray: G)
where
    T: Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync,
    G: Fn(usize, usize) -> bool + Send + Sync,
{
    let n = values.len();
    assert!(n >= 1);
    let access = values.into_par_access();
    let mut stride = 1;
    while stride < n {
        let start = stride - 1;
        (start..n - stride)
            .into_par_iter()
            .step_by(2 * stride)
            .for_each(|j| {
                let k = j + stride;
                if same_subarray(j, k) {
                    unsafe {
                        let j_ptr = access.get_unsync(j);
                        let k_ptr = access.get_unsync(k);
                        *k_ptr = op(&*k_ptr, &*j_ptr);
                    }
                }
            });
        stride *= 2;
    }
    stride /= 2;
    while stride >= 1 {
        let start = 2 * stride - 1;
        (start..n - stride)
            .into_par_iter()
            .step_by(2 * stride)
            .for_each(|j| {
                let k = j + stride;
                if same_subarray(j, k) {
                    unsafe {
                        let j_ptr = access.get_unsync(j);
                        let k_ptr = access.get_unsync(k);
                        *k_ptr = op(&*k_ptr, &*j_ptr);
                    }
                }
            });
        stride /= 2;
    }
}

pub fn suffix_sum<T, F, G>(values: &mut [T], op: F, same_subarray: G)
where
    T: Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync,
    G: Fn(usize, usize) -> bool + Send + Sync,
{
    let n = values.len();
    assert!(n >= 1);
    let access = values.into_par_access();
    let mut stride = 1;
    while stride < n {
        (0..n - stride)
            .into_par_iter()
            .step_by(2 * stride)
            .for_each(|j| {
                let k = j + stride;
                if same_subarray(j, k) {
                    unsafe {
                        let j_ptr = access.get_unsync(j);
                        let k_ptr = access.get_unsync(k);
                        *j_ptr = op(&*j_ptr, &*k_ptr);
                    }
                }
            });
        stride *= 2;
    }
    stride /= 2;
    while stride >= 1 {
        (stride..n - stride)
            .into_par_iter()
            .step_by(2 * stride)
            .for_each(|j| {
                let k = j + stride;
                if same_subarray(j, k) {
                    unsafe {
                        let j_ptr = access.get_unsync(j);
                        let k_ptr = access.get_unsync(k);
                        *j_ptr = op(&*j_ptr, &*k_ptr);
                    }
                }
            });
        stride /= 2;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..500 {
            let mut a: Vec<i32> = (0..n).map(|_| rand::random_range(0..1_000_000)).collect();
            let mut b = a.clone();
            let mut c = a.clone();
            let mut d = a.clone();

            let mut l = rand::random_range(0..n);
            let mut r = rand::random_range(0..n);
            if l > r {
                (l, r) = (r, l);
            }

            prefix_sum(&mut b, |&x, &y| x + y, |i, j| l <= i && j <= r);

            for i in l..r {
                a[i + 1] += a[i];
            }

            assert_eq!(a, b);

            suffix_sum(&mut c, |&x, &y| x + y, |i, j| l <= i && j <= r);
            for i in (l..r).rev() {
                d[i] += d[i + 1];
            }
            assert_eq!(c, d);
        }
    }
}
