//use crate::range_minimum_query::RMQ;
//use rayon::prelude::*;
//use std::cmp::{max, min};
//use rand::{thread_rng, Rng};

pub fn calculate_subtree_max<T>(
    values: &[T],
    parent: &[usize],
    time_in: &[usize],
    time_out: &[usize],
    pre_order: &[usize],
    p: usize,
) -> Vec<T>
where
    T: Ord + Clone + Send + Sync,
{
    let n = values.len();
    if n == 0 {
        return vec![];
    }

    let mut subtree_max = values.to_vec();

    let mut prefix_of_block = vec![values[0].clone(); n];
    let mut suffix_of_block = vec![values[0].clone(); n];
    for i in 0..n {
        prefix_of_block[i] = values[pre_order[i]].clone();
        suffix_of_block[i] = values[pre_order[i]].clone();
    }

    let b = (n + p - 1) / p;
    if b == 0 {
        return subtree_max;
    }

    for i in 0..p {
        let start_idx = i * b;
        if start_idx >= n {
            continue;
        }
        let end_idx = std::cmp::min((i + 1) * b, n);
        for j in (start_idx + 1)..end_idx {
            prefix_of_block[j] =
                std::cmp::max(prefix_of_block[j].clone(), prefix_of_block[j - 1].clone());
        }
        for j in (start_idx..(end_idx - 1)).rev() {
            suffix_of_block[j] =
                std::cmp::max(suffix_of_block[j].clone(), suffix_of_block[j + 1].clone());
        }
    }

    for i in 0..p {
        let start_idx = i * b;
        if start_idx >= n {
            continue;
        }
        let end_idx = std::cmp::min((i + 1) * b, n);
        for j in (start_idx..end_idx).rev() {
            let node = pre_order[j];
            let par = parent[node];
            if time_in[par] / b == (time_out[par] - 1) / b {
                subtree_max[par] =
                    std::cmp::max(subtree_max[par].clone(), subtree_max[node].clone());
            }
        }
    }

    let mut sparse_table = vec![vec![values[0].clone(); p]];
    for i in 0..p {
        if i * b < n {
            sparse_table[0][i] = suffix_of_block[i * b].clone();
        }
    }

    let mut i_row = 0;
    while (2 << i_row) <= p {
        let curr_size = p - (2 << i_row) + 1;
        let mut next_row = vec![values[0].clone(); curr_size];
        for j in 0..curr_size {
            next_row[j] = std::cmp::max(
                sparse_table[i_row][j].clone(),
                sparse_table[i_row][j + (1 << i_row)].clone(),
            );
        }
        sparse_table.push(next_row);
        i_row += 1;
    }

    let sparse_table_query = |l: usize, r: usize, sparse_table: &[Vec<T>]| -> T {
        let lg = usize::BITS as usize - 1 - (r - l).leading_zeros() as usize;
        std::cmp::max(
            sparse_table[lg][l].clone(),
            sparse_table[lg][r - (1 << lg)].clone(),
        )
    };

    for i in 0..n {
        let l = time_in[i];
        let r = time_out[i];
        assert!(l < r);

        if l / b == (r - 1) / b {
            continue;
        }

        subtree_max[i] = std::cmp::max(suffix_of_block[l].clone(), prefix_of_block[r - 1].clone());
        if l / b + 1 < r / b {
            subtree_max[i] = std::cmp::max(
                subtree_max[i].clone(),
                sparse_table_query(l / b + 1, r / b, &sparse_table),
            );
        }
    }

    subtree_max
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subtree_max_stress_test() {
        for n in 1..=100 {
            for p in 1..=(n + 5) {
                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i);
                    adjacency_list[parent[i]].push(i);
                }

                let values: Vec<i32> = (0..n).map(|_| rand::random_range(0..1_000_000)).collect();

                let mut time_in = vec![0; n];
                let mut time_out = vec![0; n];
                let mut pre_order = vec![0; n];

                {
                    let mut timer = 0;
                    fn dfs(
                        node: usize,
                        timer: &mut usize,
                        time_in: &mut [usize],
                        time_out: &mut [usize],
                        pre_order: &mut [usize],
                        adjacency_list: &[Vec<usize>],
                    ) {
                        time_in[node] = *timer;
                        pre_order[*timer] = node;
                        *timer += 1;
                        for &child in &adjacency_list[node] {
                            dfs(child, timer, time_in, time_out, pre_order, adjacency_list);
                        }
                        time_out[node] = *timer;
                    }
                    dfs(
                        0,
                        &mut timer,
                        &mut time_in,
                        &mut time_out,
                        &mut pre_order,
                        &adjacency_list,
                    );
                }

                let subtree_max_block =
                    calculate_subtree_max(&values, &parent, &time_in, &time_out, &pre_order, p);

                let mut expected = vec![0; n];
                fn compute_naive(
                    node: usize,
                    adjacency_list: &[Vec<usize>],
                    values: &[i32],
                    expected: &mut [i32],
                ) -> i32 {
                    let mut curr = values[node];
                    for &child in &adjacency_list[node] {
                        curr = std::cmp::max(
                            curr,
                            compute_naive(child, adjacency_list, values, expected),
                        );
                    }
                    expected[node] = curr;
                    curr
                }
                compute_naive(0, &adjacency_list, &values, &mut expected);

                assert_eq!(subtree_max_block, expected);
            }
        }
    }
}
