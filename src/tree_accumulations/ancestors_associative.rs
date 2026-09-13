//use crate::helpers::sparse_table::SparseTable;
//use paradis_core::{BoundedParAccess, IntoParAccess};
//use rayon::prelude::*;

/*
pub fn ancestors_associative<T, F>(
    values: &[T],
    op: F,
    parent: &[usize],
    time_in: &[usize],
    time_out: &[usize],
    pre_order: &[usize],
    p: usize,
) -> Vec<T>
where
    T: Clone + Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync,
{
    let n = values.len();
    assert!(n >= 1 && p >= 1);

}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        type Matrix = [[i64; 2]; 2];
        const IDENTITY: Matrix = [[1, 0], [0, 1]];
        let mult = |a: Matrix, b: Matrix, mod_val: i64| -> Matrix {
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

        for n in 1..=100 {
            let max_p = 2 * n + 5;
            for p in 1..=max_p {
                println!("n,P: {} {}", n, p);

                let mod_val = rand::random_range(10..1_000_000_000) as i64;

                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i) as usize;
                    adjacency_list[parent[i]].push(i);
                }

                let mut value = vec![[[0_i64; 2]; 2]; n];
                for i in 0..n {
                    value[i][0][0] = rand::random_range(0..mod_val) as i64;
                    value[i][0][1] = rand::random_range(0..mod_val) as i64;
                    value[i][1][0] = rand::random_range(0..mod_val) as i64;
                    value[i][1][1] = rand::random_range(0..mod_val) as i64;
                }

                let mut time_in = vec![0; n];
                let mut time_out = vec![0; n];
                let mut euler_tour = vec![0; 2 * n - 1];
                let mut depth = vec![0; n];

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
                    ) {
                        time_in[node] = *timer;
                        euler_tour[*timer] = node;
                        *timer += 1;
                        for &child in &adjacency_list[node] {
                            depth[child] = 1 + depth[node];
                            dfs(
                                child,
                                timer,
                                adjacency_list,
                                time_in,
                                time_out,
                                euler_tour,
                                depth,
                            );
                            euler_tour[*timer] = node;
                            *timer += 1;
                        }
                        time_out[node] = *timer - 1;
                    }

                    dfs(
                        0,
                        &mut timer,
                        &adjacency_list,
                        &mut time_in,
                        &mut time_out,
                        &mut euler_tour,
                        &mut depth,
                    );
                    assert_eq!(timer, 2 * n - 1);
                }

                for i in 0..n {
                    assert_eq!(euler_tour[time_in[i]], i);
                    assert_eq!(euler_tour[time_out[i]], i);
                    assert!(time_in[i] <= time_out[i]);
                }

                let b = (2 * n - 1 + p - 1) / p;
                let mut ancestor_agg = value.clone();

                for i in 0..p {
                    let start_idx = i * b;
                    if start_idx >= 2 * n - 1 {
                        continue;
                    }
                    let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);

                    for j in start_idx..end_idx {
                        let node = euler_tour[j];
                        let par = parent[node];
                        if time_in[node] == j && par != node {
                            if time_in[par] / b == time_out[par] / b {
                                ancestor_agg[node] =
                                    mult(ancestor_agg[node], ancestor_agg[par], mod_val);
                            }
                        }
                    }
                }

                {
                    let mut prefix_agg = vec![IDENTITY; 2 * n - 1];
                    let mut suffix_agg = vec![IDENTITY; 2 * n - 1];

                    for i in 0..n {
                        let l = time_in[i];
                        let r = time_out[i];
                        assert!(l <= r);
                        if l / b == r / b {
                            continue;
                        }
                        assert!(suffix_agg[l] == IDENTITY);
                        suffix_agg[l] = value[i];
                        assert!(prefix_agg[r] == IDENTITY);
                        prefix_agg[r] = value[i];
                    }

                    for i in 0..p {
                        let start_idx = i * b;
                        if start_idx >= 2 * n - 1 {
                            continue;
                        }
                        let end_idx = std::cmp::min((i + 1) * b, 2 * n - 1);

                        {
                            let mut running_agg = IDENTITY;
                            for j in start_idx..end_idx {
                                let sub_node = euler_tour[j];
                                if time_in[sub_node] == j {
                                    ancestor_agg[sub_node] =
                                        mult(ancestor_agg[sub_node], running_agg, mod_val);
                                }
                                running_agg = mult(suffix_agg[j], running_agg, mod_val);
                            }
                        }

                        {
                            let mut running_agg = IDENTITY;
                            for j in (start_idx..end_idx).rev() {
                                running_agg = mult(prefix_agg[j], running_agg, mod_val);
                                let sub_node = euler_tour[j];
                                if time_in[sub_node] == j {
                                    ancestor_agg[sub_node] =
                                        mult(ancestor_agg[sub_node], running_agg, mod_val);
                                }
                            }
                        }
                    }
                }

                {
                    let mut disjoint_rmq = vec![vec![IDENTITY; p]; (p.ilog2() as usize) + 1];

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
                                    disjoint_rmq[lg][block_l] =
                                        mult(disjoint_rmq[lg][block_l], value[node], mod_val);
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
                                        mult(disjoint_rmq[lg][block_r - 1], value[node], mod_val);
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
                                    disjoint_rmq[i][k] =
                                        mult(disjoint_rmq[i][k], disjoint_rmq[i][j], mod_val);
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
                                    disjoint_rmq[i][k] =
                                        mult(disjoint_rmq[i][k], disjoint_rmq[i][j], mod_val);
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
                                    disjoint_rmq[i][j] =
                                        mult(disjoint_rmq[i][j], disjoint_rmq[i][k], mod_val);
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
                                    disjoint_rmq[i][j] =
                                        mult(disjoint_rmq[i][j], disjoint_rmq[i][k], mod_val);
                                }
                                j += 2 * stride2;
                            }
                            stride2 /= 2;
                        }
                    }

                    let mut rmq_agg = vec![IDENTITY; p];
                    for i in 0..disjoint_rmq.len() {
                        for j in 0..p {
                            rmq_agg[j] = mult(rmq_agg[j], disjoint_rmq[i][j], mod_val);
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
                            if time_in[node] == j {
                                ancestor_agg[node] = mult(ancestor_agg[node], rmq_agg[i], mod_val);
                            }
                        }
                    }
                }

                fn verify_dfs(
                    node: usize,
                    agg_up: Matrix,
                    adjacency_list: &Vec<Vec<usize>>,
                    ancestor_agg: &Vec<Matrix>,
                    value: &Vec<Matrix>,
                    mod_val: i64,
                ) {
                    for i in 0..2 {
                        for j in 0..2 {
                            assert_eq!(ancestor_agg[node][i][j], agg_up[i][j]);
                        }
                    }

                    let mult_fn = |a: Matrix, b: Matrix| -> Matrix {
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

                    for &child in &adjacency_list[node] {
                        verify_dfs(
                            child,
                            mult_fn(value[child], agg_up),
                            adjacency_list,
                            ancestor_agg,
                            value,
                            mod_val,
                        );
                    }
                }

                verify_dfs(0, value[0], &adjacency_list, &ancestor_agg, &value, mod_val);
            }
        }
    }
}
