//use crate::helpers::modified_blelloch_scan::*;
use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::prelude::*;

#[allow(clippy::too_many_arguments)]
pub fn ancestors_commutative_associative_idempotent<T, F>(
    values: &[T],
    identity: T,
    op: F,
    parent: &[usize],
    euler_tour: &[usize],
    et_time_in: &[usize],
    et_time_out: &[usize],
    p: usize,
) -> Vec<T>
where
    T: Clone + Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync,
{
    let n = values.len();
    assert!(n >= 1 && p >= 1);

    let b = (2 * n - 1).div_ceil(p);

    let mut ancestor_agg = values.to_vec();
    {
        let access = ancestor_agg.into_par_access();
        (0..p).into_par_iter().for_each(|i| {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                return;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);

            for &node in &euler_tour[start_idx..end_idx] {
                let par = parent[node];
                if par != node && et_time_in[par] / b == et_time_out[par] / b {
                    unsafe {
                        let node_ptr = access.get_unsync(node);
                        let par_ptr = access.get_unsync(par);
                        *node_ptr = op(&*node_ptr, &*par_ptr);
                    }
                }
            }
        });
    }

    let mut euler_tour_max: Vec<T> = euler_tour
        .par_iter()
        .map(|&node| values[node].clone())
        .collect();

    let mut inverse_sparse_table = vec![vec![identity.clone(); p]; p.ilog2() as usize + 2];

    {
        let access = inverse_sparse_table.into_par_access();
        (0..p).into_par_iter().for_each(|i| {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                return;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);
            for &node in &euler_tour[start_idx..end_idx] {
                let num_middle_right = ((et_time_out[node] + 1) / b) as isize - (i + 1) as isize;
                if num_middle_right >= 1 {
                    let lg = (num_middle_right as usize).ilog2() as usize;
                    unsafe {
                        let row_ptr = access.get_unsync(lg);
                        (*row_ptr)[i + 1] = op(&(*row_ptr)[i + 1], &values[node]);
                    }
                }
            }
        });
        (0..p).into_par_iter().for_each(|i| {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                return;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);
            for &node in &euler_tour[start_idx..end_idx] {
                let num_middle_left = i as isize - (et_time_in[node] / b + 1) as isize;
                if num_middle_left >= 1 {
                    let lg = (num_middle_left as usize).ilog2() as usize;
                    let col = i - (1 << lg);
                    unsafe {
                        let row_ptr = access.get_unsync(lg);
                        (*row_ptr)[col] = op(&(*row_ptr)[col], &values[node]);
                    }
                }
            }
        });
    }

    for l in (1..inverse_sparse_table.len()).rev() {
        let (lower, upper) = inverse_sparse_table.split_at_mut(l);
        let current_row = &mut lower[l - 1];
        let prev_row = &upper[0];

        current_row.par_iter_mut().enumerate().for_each(|(j, val)| {
            if j < prev_row.len() {
                *val = op(val, &prev_row[j]);
            }
            if j >= (1 << (l - 1)) {
                let prev_idx = j - (1 << (l - 1));
                if prev_idx < prev_row.len() {
                    *val = op(val, &prev_row[prev_idx]);
                }
            }
        });
    }

    let mut prefix_agg = vec![identity.clone(); 2 * n - 1];
    let mut suffix_agg = vec![identity.clone(); 2 * n - 1];
    {
        let pref_access = prefix_agg.into_par_access();
        let suf_access = suffix_agg.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            let l = et_time_in[i];
            let r = et_time_out[i];
            assert!(l <= r);
            if l / b == r / b {
                return;
            }
            unsafe {
                *suf_access.get_unsync(l) = values[i].clone();
                *pref_access.get_unsync(r) = values[i].clone();
            }
        });
    }

    euler_tour_max
        .par_chunks_mut(b)
        .enumerate()
        .for_each(|(i, chunk)| {
            let start_idx = i * b;
            let end_idx = start_idx + chunk.len();
            let block_val = &inverse_sparse_table[0][i];

            for val in chunk.iter_mut() {
                *val = op(val, block_val);
            }

            let mut running_agg = identity.clone();
            for (offset, j) in (start_idx..end_idx).enumerate() {
                running_agg = op(&running_agg, &suffix_agg[j]);
                chunk[offset] = op(&chunk[offset], &running_agg);
            }

            let mut running_agg = identity.clone();
            for (offset, j) in (start_idx..end_idx).enumerate().rev() {
                running_agg = op(&running_agg, &prefix_agg[j]);
                chunk[offset] = op(&chunk[offset], &running_agg);
            }
        });

    {
        let access = ancestor_agg.into_par_access();
        (0..p).into_par_iter().for_each(|i| {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                return;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);
            for j in start_idx..end_idx {
                let node = euler_tour[j];
                if et_time_in[node] == j {
                    unsafe {
                        let node_ptr = access.get_unsync(node);
                        *node_ptr = op(&*node_ptr, &euler_tour_max[j]);
                    }
                }
            }
        });
    }

    ancestor_agg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=80 {
            for p in 1..=n + 5 {
                println!("n,P: {} {}", n, p);

                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i);
                    adjacency_list[parent[i]].push(i);
                }

                let value: Vec<i32> = (0..n)
                    .map(|_| rand::random_range(-1_000_000..=1_000_000))
                    .collect();

                let mut et_time_in = vec![0; n];
                let mut et_time_out = vec![0; n];
                let mut euler_tour = vec![0; 2 * n - 1];
                let mut ancestor_agg_naive = value.clone();

                {
                    let mut timer = 0;
                    fn dfs(
                        node: usize,
                        timer: &mut usize,
                        adjacency_list: &[Vec<usize>],
                        et_time_in: &mut [usize],
                        et_time_out: &mut [usize],
                        euler_tour: &mut [usize],
                        ancestor_agg_naive: &mut [i32],
                    ) {
                        et_time_in[node] = *timer;
                        euler_tour[*timer] = node;
                        *timer += 1;
                        for &child in &adjacency_list[node] {
                            ancestor_agg_naive[child] =
                                ancestor_agg_naive[child].max(ancestor_agg_naive[node]);
                            dfs(
                                child,
                                timer,
                                adjacency_list,
                                et_time_in,
                                et_time_out,
                                euler_tour,
                                ancestor_agg_naive,
                            );
                            euler_tour[*timer] = node;
                            *timer += 1;
                        }
                        et_time_out[node] = *timer - 1;
                    }

                    dfs(
                        0,
                        &mut timer,
                        &adjacency_list,
                        &mut et_time_in,
                        &mut et_time_out,
                        &mut euler_tour,
                        &mut ancestor_agg_naive,
                    );
                    assert_eq!(timer, 2 * n - 1);
                }

                for i in 0..n {
                    assert_eq!(euler_tour[et_time_in[i]], i);
                    assert_eq!(euler_tour[et_time_out[i]], i);
                    assert!(et_time_in[i] <= et_time_out[i]);
                }

                let ancestor_agg = ancestors_commutative_associative_idempotent(
                    &value,
                    i32::MIN,
                    |a, b| *a.max(b),
                    &parent,
                    &euler_tour,
                    &et_time_in,
                    &et_time_out,
                    p,
                );

                assert_eq!(ancestor_agg, ancestor_agg_naive);
            }
        }
    }
}
