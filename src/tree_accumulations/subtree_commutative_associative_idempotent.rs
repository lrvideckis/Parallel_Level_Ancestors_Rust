use crate::sparse_tables::sparse_table::SparseTable;
use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::prelude::*;

pub fn subtree_commutative_associative_idempotent<T, F>(
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

    let mut prefix_of_block: Vec<T> = (0..n)
        .into_par_iter()
        .map(|i| values[pre_order[i]].clone())
        .collect();

    let mut suffix_of_block = prefix_of_block.clone();

    let b = n.div_ceil(p);
    assert!(b >= 1);

    prefix_of_block.par_chunks_mut(b).for_each(|chunk| {
        let len = chunk.len();
        for j in 1..len {
            chunk[j] = op(&chunk[j], &chunk[j - 1]);
        }
    });
    suffix_of_block.par_chunks_mut(b).for_each(|chunk| {
        let len = chunk.len();
        for j in (1..len).rev() {
            chunk[j - 1] = op(&chunk[j - 1], &chunk[j]);
        }
    });

    let mut subtree_max = values.to_vec();
    let access = subtree_max.into_par_access();
    (0..p).into_par_iter().for_each(|i| {
        let start_idx = i * b;
        if start_idx >= n {
            return;
        }
        let end_idx = std::cmp::min((i + 1) * b, n);
        for j in (start_idx..end_idx).rev() {
            let node = pre_order[j];
            let par = parent[node];
            if time_in[par] / b == (time_out[par] - 1) / b {
                unsafe {
                    let par_ref = access.get_unsync(par);
                    let node_ref = access.get_unsync(node);
                    *par_ref = op(par_ref, node_ref);
                }
            }
        }
    });

    let block_values: Vec<T> = (0..p)
        .into_par_iter()
        .map(|i| {
            if i * b < n {
                suffix_of_block[i * b].clone()
            } else {
                values[0].clone()
            }
        })
        .collect();

    let sparse_table = SparseTable::new(block_values, |x, y| op(x, y));

    subtree_max.par_iter_mut().enumerate().for_each(|(i, val)| {
        let l = time_in[i];
        let r = time_out[i];
        assert!(l < r);

        if l / b == (r - 1) / b {
            return;
        }

        *val = op(&suffix_of_block[l], &prefix_of_block[r - 1]);
        if l / b + 1 < r / b {
            *val = op(val, &sparse_table.query(l / b + 1..r / b));
        }
    });

    subtree_max
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
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
                let mut subtree_max_naive = values.clone();

                {
                    let mut timer = 0;
                    fn dfs(
                        node: usize,
                        timer: &mut usize,
                        time_in: &mut [usize],
                        time_out: &mut [usize],
                        pre_order: &mut [usize],
                        subtree_max_naive: &mut [i32],
                        adjacency_list: &[Vec<usize>],
                    ) {
                        time_in[node] = *timer;
                        pre_order[*timer] = node;
                        *timer += 1;
                        for &child in &adjacency_list[node] {
                            dfs(
                                child,
                                timer,
                                time_in,
                                time_out,
                                pre_order,
                                subtree_max_naive,
                                adjacency_list,
                            );
                            subtree_max_naive[node] =
                                std::cmp::max(subtree_max_naive[node], subtree_max_naive[child]);
                        }
                        time_out[node] = *timer;
                    }
                    dfs(
                        0,
                        &mut timer,
                        &mut time_in,
                        &mut time_out,
                        &mut pre_order,
                        &mut subtree_max_naive,
                        &adjacency_list,
                    );
                }

                let subtree_max = subtree_commutative_associative_idempotent(
                    &values,
                    |&x, &y| std::cmp::max(x, y),
                    &parent,
                    &time_in,
                    &time_out,
                    &pre_order,
                    p,
                );

                assert_eq!(subtree_max, subtree_max_naive);
            }
        }
    }
}
