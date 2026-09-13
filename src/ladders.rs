use crate::logn_query_methods::binary_lifting::BinaryLifting;
use crate::tree_accumulations::subtree_max::calculate_subtree_max;
use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::iter::once;
use rayon::prelude::*;
use rayon_scan::ScanParallelIterator;

pub struct Ladders {
    level: Vec<usize>,
    ladder: Vec<usize>,
    deepest_leaf: Vec<usize>,
    leaf_to_size: Vec<usize>,
    leaf_to_start: Vec<usize>,
}

impl Ladders {
    pub fn new(
        parent: &[usize],
        level: &[usize],
        time_in: &[usize],
        time_out: &[usize],
        pre_order: &[usize],
        p: usize,
    ) -> Self {
        let n = parent.len();
        assert!(n >= 1 && p >= 1);

        let values: Vec<(usize, usize)> = (0..n).into_par_iter().map(|i| (level[i], i)).collect();

        let deepest_leaf = calculate_subtree_max(&values, parent, time_in, time_out, pre_order, p);
        let deepest_leaf: Vec<usize> = deepest_leaf.into_par_iter().map(|t| t.1).collect();

        let mut leaf_to_size = vec![0; n];
        let access = leaf_to_size.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            if parent[i] == i || deepest_leaf[parent[i]] != deepest_leaf[i] {
                unsafe {
                    let target_ref = access.get_unsync(deepest_leaf[i]);
                    *target_ref = 2 * (level[deepest_leaf[i]] - level[i] + 1);
                }
            }
        });
        let leaf_to_size = leaf_to_size;

        assert_eq!(leaf_to_size.iter().sum::<usize>(), 2 * n);

        let leaf_to_start: Vec<usize> = once(0)
            .chain(leaf_to_size.par_iter().cloned())
            .scan(|a, b| *a + *b, 0)
            .collect();

        let mut ladder_to_start = vec![0; 2 * n + 1];
        let access = ladder_to_start.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            if deepest_leaf[i] == i {
                let idx = leaf_to_start[i] + leaf_to_size[i];
                unsafe {
                    let target_ref = access.get_unsync(idx);
                    *target_ref = leaf_to_size[i];
                }
            }
        });

        let ladder_to_start: Vec<usize> = ladder_to_start
            .into_par_iter()
            .scan(|a, b| *a + *b, 0)
            .collect();
        assert_eq!(ladder_to_start[2 * n], 2 * n);

        let mut ladder = vec![usize::MAX; 2 * n];
        let access = ladder.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            if deepest_leaf[i] == i {
                let idx = leaf_to_start[i];
                unsafe {
                    let target_ref = access.get_unsync(idx);
                    *target_ref = i;
                }
            }
        });

        let jump_pointers = BinaryLifting::new(parent, level, p);

        let block_size = (2 * n).div_ceil(p);

        let ladder_temp = ladder.clone();
        ladder
            .par_chunks_mut(block_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                let start_index = i * block_size;
                let leaf = ladder_temp[ladder_to_start[start_index]];
                let k = std::cmp::min(start_index - ladder_to_start[start_index], level[leaf]);
                chunk[0] = jump_pointers.kth_parent(leaf, k);
                for j in 1..chunk.len() {
                    if chunk[j] == usize::MAX {
                        chunk[j] = parent[chunk[j - 1]];
                    }
                }
            });

        Self {
            level: level.to_vec(),
            deepest_leaf,
            ladder,
            leaf_to_size,
            leaf_to_start,
        }
    }

    pub fn kth_parent(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        let leaf = self.deepest_leaf[v];
        let difference = self.level[leaf] - self.level[v];
        // fails when ladder is not long enough
        assert!(difference + k < self.leaf_to_size[leaf]);
        self.ladder[self.leaf_to_start[leaf] + difference + k]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=80 {
            for p in 1..=(n + 5) {
                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                let mut level = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i);
                    level[i] = 1 + level[parent[i]];
                    adjacency_list[parent[i]].push(i);
                }

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

                let ladder = Ladders::new(&parent, &level, &time_in, &time_out, &pre_order, p);
                for i in 0..n {
                    assert!(ladder.kth_parent(i, 0) == i);
                    if i > 0 {
                        assert!(ladder.kth_parent(i, 1) == parent[i]);
                    }
                    let mut u = i;
                    let mut ancestors = vec![];
                    for _ in 0..=level[i] {
                        ancestors.push(u);
                        u = parent[u];
                    }
                    assert!(ancestors.len() == level[i] + 1);
                    for k in 0..=level[i] / 2 {
                        assert!(ladder.kth_parent(ancestors[k], k) == ancestors[2 * k]);
                    }
                }
            }
        }
    }
}
