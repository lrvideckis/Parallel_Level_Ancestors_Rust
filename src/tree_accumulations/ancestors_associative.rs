use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::prelude::*;

#[allow(clippy::too_many_arguments)]
pub fn ancestors_associative<T, F>(
    values: &[T],
    identity: T,
    op: F,
    parent: &[usize],
    time_in: &[usize],
    time_out: &[usize],
    euler_tour: &[usize],
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
    let access = ancestor_agg.into_par_access();
    (0..p).into_par_iter().for_each(|i| {
        let start_idx = i * b;
        if start_idx >= 2 * n - 1 {
            return;
        }
        let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);
        for j in start_idx..end_idx {
            let node = euler_tour[j];
            let par = parent[node];
            if time_in[node] == j && par != node {
                if time_in[par] / b == time_out[par] / b {
                    unsafe {
                        let node_ptr = access.get_unsync(node);
                        let par_ptr = access.get_unsync(par);
                        *node_ptr = op(&*node_ptr, &*par_ptr);
                    }
                }
            }
        }
    });

    {
        let mut prefix_agg = vec![identity.clone(); 2 * n - 1];
        let mut suffix_agg = vec![identity.clone(); 2 * n - 1];
        let pref_access = prefix_agg.into_par_access();
        let suf_access = suffix_agg.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            let l = time_in[i];
            let r = time_out[i];
            assert!(l <= r);
            if l / b == r / b {
                return;
            }
            unsafe {
                *suf_access.get_unsync(l) = values[i].clone();
                *pref_access.get_unsync(r) = values[i].clone();
            }
        });

        (0..p).into_par_iter().for_each(|i| {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                return;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);

            {
                let mut running_agg = identity.clone();
                for j in start_idx..end_idx {
                    let sub_node = euler_tour[j];
                    if time_in[sub_node] == j {
                        unsafe {
                            let node_ptr = access.get_unsync(sub_node);
                            *node_ptr = op(&*node_ptr, &running_agg);
                        }
                    }
                    running_agg = op(&suffix_agg[j], &running_agg);
                }
            }
            {
                let mut running_agg = identity.clone();
                for j in (start_idx..end_idx).rev() {
                    running_agg = op(&prefix_agg[j], &running_agg);
                    let sub_node = euler_tour[j];
                    if time_in[sub_node] == j {
                        unsafe {
                            let node_ptr = access.get_unsync(sub_node);
                            *node_ptr = op(&*node_ptr, &running_agg);
                        }
                    }
                }
            }
        });
    }

    {
        let mut disjoint_rmq = vec![vec![identity.clone(); p]; (p.ilog2() as usize) + 1];

        for i in 0..p {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                continue;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);

            for l_idx in (start_idx..end_idx).rev() {
                let node = euler_tour[l_idx];
                let r = time_out[node];
                if time_in[node] == l_idx {
                    let block_l = i + 1;
                    let block_r = r / b;
                    if block_l < block_r {
                        let mut lg = 0;
                        if block_l + 1 < block_r {
                            lg = (block_l ^ (block_r - 1)).ilog2() as usize;
                        }
                        disjoint_rmq[lg][block_l] = op(&disjoint_rmq[lg][block_l], &values[node]);
                    }
                }
            }
        }

        for i in 0..p {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                continue;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);

            for j in start_idx..end_idx {
                let node = euler_tour[j];
                let l = time_in[node];
                let r = time_out[node];
                if r == j {
                    let block_l = l / b + 1;
                    let block_r = i;
                    if block_l + 1 < block_r {
                        let lg = (block_l ^ (block_r - 1)).ilog2() as usize;
                        disjoint_rmq[lg][block_r - 1] =
                            op(&disjoint_rmq[lg][block_r - 1], &values[node]);
                    }
                }
            }
        }

        for i in 0..disjoint_rmq.len() {
            let mut stride = 1;
            while stride < p {
                let mut j = stride - 1;
                while j + stride < p {
                    let k = j + stride;
                    if (j >> i) == (k >> i) && ((j >> i) % 2 == 0) {
                        disjoint_rmq[i][k] = op(&disjoint_rmq[i][k], &disjoint_rmq[i][j]);
                    }
                    j += 2 * stride;
                }
                stride *= 2;
            }
            stride /= 2;
            while stride >= 1 {
                let mut j = 2 * stride - 1;
                while j + stride < p {
                    let k = j + stride;
                    if (j >> i) == (k >> i) && ((j >> i) % 2 == 0) {
                        disjoint_rmq[i][k] = op(&disjoint_rmq[i][k], &disjoint_rmq[i][j]);
                    }
                    j += 2 * stride;
                }
                stride /= 2;
            }

            let mut stride2 = 1;
            while stride2 < p {
                let mut j = 0;
                while j + stride2 < p {
                    let k = j + stride2;
                    if (j >> i) == (k >> i) && ((j >> i) % 2 == 1) {
                        disjoint_rmq[i][j] = op(&disjoint_rmq[i][j], &disjoint_rmq[i][k]);
                    }
                    j += 2 * stride2;
                }
                stride2 *= 2;
            }
            stride2 /= 2;
            while stride2 >= 1 {
                let mut j = stride2;
                while j + stride2 < p {
                    let k = j + stride2;
                    if (j >> i) == (k >> i) && ((j >> i) % 2 == 1) {
                        disjoint_rmq[i][j] = op(&disjoint_rmq[i][j], &disjoint_rmq[i][k]);
                    }
                    j += 2 * stride2;
                }
                stride2 /= 2;
            }
        }

        let rmq_agg: Vec<T> = (0..p)
            .into_par_iter()
            .map(|j| {
                let mut acc = identity.clone();
                for i in 0..disjoint_rmq.len() {
                    acc = op(&acc, &disjoint_rmq[i][j]);
                }
                acc
            })
            .collect();

        let access = ancestor_agg.into_par_access();
        (0..p).into_par_iter().for_each(|i| {
            let start_idx = i * b;
            if start_idx >= 2 * n - 1 {
                return;
            }
            let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);
            for j in start_idx..end_idx {
                let node = euler_tour[j];
                if time_in[node] == j {
                    unsafe {
                        let node_ptr = access.get_unsync(node);
                        *node_ptr = op(&*node_ptr, &rmq_agg[i]);
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
        type Matrix = [[u64; 2]; 2];
        const IDENTITY: Matrix = [[1, 0], [0, 1]];
        let mult = |a: &Matrix, b: &Matrix, mod_val: u64| -> Matrix {
            let mut res: Matrix = [[0, 0], [0, 0]];
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        res[i][j] = (res[i][j] + a[i][k] * b[k][j]) % mod_val;
                    }
                }
            }
            res
        };

        for n in 1..=80 {
            for p in 1..=2 * n + 5 {
                println!("n,P: {} {}", n, p);

                let mod_val = rand::random_range(10..1_000_000_000) as u64;

                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i) as usize;
                    adjacency_list[parent[i]].push(i);
                }

                let mut value = vec![[[0_u64; 2]; 2]; n];
                for i in 0..n {
                    value[i][0][0] = rand::random_range(0..mod_val) as u64;
                    value[i][0][1] = rand::random_range(0..mod_val) as u64;
                    value[i][1][0] = rand::random_range(0..mod_val) as u64;
                    value[i][1][1] = rand::random_range(0..mod_val) as u64;
                }

                let mut time_in = vec![0; n];
                let mut time_out = vec![0; n];
                let mut euler_tour = vec![0; 2 * n - 1];
                let mut depth = vec![0; n];
                let mut ancestor_agg_naive = vec![IDENTITY; n];

                {
                    let mut timer = 0;
                    fn dfs(
                        node: usize,
                        timer: &mut usize,
                        adjacency_list: &Vec<Vec<usize>>,
                        time_in: &mut Vec<usize>,
                        time_out: &mut Vec<usize>,
                        euler_tour: &mut Vec<usize>,
                        depth: &mut Vec<usize>,
                        value: &Vec<Matrix>,
                        ancestor_agg_naive: &mut Vec<Matrix>,
                        mod_val: u64,
                    ) {
                        time_in[node] = *timer;
                        euler_tour[*timer] = node;
                        *timer += 1;
                        for &child in &adjacency_list[node] {
                            depth[child] = 1 + depth[node];

                            let a = value[child];
                            let b = ancestor_agg_naive[node];
                            let mut res: Matrix = [[0, 0], [0, 0]];
                            for i in 0..2 {
                                for j in 0..2 {
                                    for k in 0..2 {
                                        res[i][j] = (res[i][j] + a[i][k] * b[k][j]) % mod_val;
                                    }
                                }
                            }
                            ancestor_agg_naive[child] = res;

                            dfs(
                                child,
                                timer,
                                adjacency_list,
                                time_in,
                                time_out,
                                euler_tour,
                                depth,
                                value,
                                ancestor_agg_naive,
                                mod_val,
                            );
                            euler_tour[*timer] = node;
                            *timer += 1;
                        }
                        time_out[node] = *timer - 1;
                    }

                    ancestor_agg_naive[0] = value[0];
                    dfs(
                        0,
                        &mut timer,
                        &adjacency_list,
                        &mut time_in,
                        &mut time_out,
                        &mut euler_tour,
                        &mut depth,
                        &value,
                        &mut ancestor_agg_naive,
                        mod_val,
                    );
                    assert_eq!(timer, 2 * n - 1);
                }

                for i in 0..n {
                    assert_eq!(euler_tour[time_in[i]], i);
                    assert_eq!(euler_tour[time_out[i]], i);
                    assert!(time_in[i] <= time_out[i]);
                }

                let ancestor_agg = ancestors_associative(
                    &value,
                    IDENTITY,
                    |a, b| mult(a, b, mod_val),
                    &parent,
                    &time_in,
                    &time_out,
                    &euler_tour,
                    p,
                );

                assert_eq!(ancestor_agg, ancestor_agg_naive);
            }
        }
    }
}
