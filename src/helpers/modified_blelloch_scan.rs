//use rayon::prelude::*;

pub fn prefix_sum<T, F, G>(values: &mut [T], op: F, same_subarray: G)
where
    T: Clone + Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync,
    G: Fn(usize, usize) -> bool + Send + Sync,
{
    let n = values.len();
    if n == 0 {
        return;
    }

    let mut stride = 1;
    while stride < n {
        let mut j = stride - 1;
        while j + stride < n {
            let k = j + stride;
            if same_subarray(j, k) {
                values[k] = op(&values[k], &values[j]);
            }
            j += 2 * stride;
        }
        stride *= 2;
    }

    stride /= 2;
    while stride >= 1 {
        let mut j = 2 * stride - 1;
        while j + stride < n {
            let k = j + stride;
            if same_subarray(j, k) {
                values[k] = op(&values[k], &values[j]);
            }
            j += 2 * stride;
        }
        stride /= 2;
    }
}

pub fn suffix_sum<T, F, G>(values: &mut [T], op: F, same_subarray: G)
where
    T: Clone + Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync,
    G: Fn(usize, usize) -> bool + Send + Sync,
{
    let n = values.len();
    if n == 0 {
        return;
    }

    let mut stride2 = 1;
    while stride2 < n {
        let mut j = 0;
        while j + stride2 < n {
            let k = j + stride2;
            if same_subarray(j, k) {
                values[j] = op(&values[j], &values[k]);
            }
            j += 2 * stride2;
        }
        stride2 *= 2;
    }

    stride2 /= 2;
    while stride2 >= 1 {
        let mut j = stride2;
        while j + stride2 < n {
            let k = j + stride2;
            if same_subarray(j, k) {
                values[j] = op(&values[j], &values[k]);
            }
            j += 2 * stride2;
        }
        stride2 /= 2;
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
