use crate::n_build_logn_query_methods::jump_pointers::JumpPointersModedLevel;
use crate::tree_accumulations::subtree_max::calculate_subtree_max;
//use paradis_core::{BoundedParAccess, IntoParAccess};
//use rayon::prelude::*;

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

        let mut ladder = vec![usize::MAX; 2 * n];
        let mut leaf_to_size = vec![0; n];
        let mut ladder_to_start = vec![0; 2 * n + 1];

        let values: Vec<(usize, usize)> = (0..n).map(|i| (level[i], i)).collect();

        let deepest_leaf = calculate_subtree_max(&values, parent, time_in, time_out, pre_order, p);
        let deepest_leaf: Vec<usize> = deepest_leaf.into_iter().map(|t| t.1).collect();

        leaf_to_size[deepest_leaf[0]] = 2 * (level[deepest_leaf[0]] + 1);
        for i in 0..n {
            if parent[i] != i {
                let par = parent[i];
                let node = i;
                if deepest_leaf[par] != deepest_leaf[node] {
                    leaf_to_size[deepest_leaf[node]] =
                        2 * (level[deepest_leaf[node]] - level[node] + 1);
                }
            }
        }

        let total_sum: usize = leaf_to_size.iter().sum();
        assert_eq!(total_sum, 2 * n);

        let mut leaf_to_start = vec![0; n];
        let mut acc = 0;
        for i in 0..n {
            leaf_to_start[i] = acc;
            acc += leaf_to_size[i];
        }

        for i in 0..n {
            if deepest_leaf[i] == i {
                ladder[leaf_to_start[i]] = i;
            }
        }

        for i in 0..n {
            if deepest_leaf[i] == i {
                ladder_to_start[leaf_to_start[i] + leaf_to_size[i]] = leaf_to_size[i];
            }
        }

        let mut acc = 0;
        for val in ladder_to_start.iter_mut() {
            acc += *val;
            *val = acc;
        }
        assert_eq!(ladder_to_start[2 * n], 2 * n);

        let jump_pointers = JumpPointersModedLevel::new(parent, level, p);

        let block_size = (2 * n + p - 1) / p;

        for i in 0..p {
            let start_index = i * block_size;
            if start_index >= 2 * n {
                continue;
            }
            let end_index = std::cmp::min(2 * n, start_index + block_size);

            let ladder_idx = ladder_to_start[start_index];
            let max_jump = level[ladder[ladder_idx]];
            let kth_offset = std::cmp::min(start_index - ladder_to_start[start_index], max_jump);
            ladder[start_index] = jump_pointers.query(ladder[ladder_idx], kth_offset);

            for j in (start_index + 1)..end_index {
                if ladder[j] == usize::MAX {
                    let prev = ladder[j - 1];
                    ladder[j] = parent[prev];
                }
            }
        }

        for i in 0..n {
            let leaf = deepest_leaf[i];
            let diff = level[leaf] - level[i];
            assert_eq!(ladder[leaf_to_start[leaf] + diff], i);
        }

        Self {
            level: level.to_vec(),
            deepest_leaf,
            ladder,
            leaf_to_size,
            leaf_to_start,
        }
    }

    pub fn query(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        let anc_d = self.level[v] - k;
        let leaf = self.deepest_leaf[v];
        assert!(self.level[leaf] - anc_d < self.leaf_to_size[leaf]);
        return self.ladder[self.leaf_to_start[leaf] + self.level[leaf] - anc_d];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ladders_stress_test() {
        for n in 1..=100 {
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
                    let mut u = i;
                    let mut ancestors = vec![];
                    for _ in 0..=level[i] {
                        ancestors.push(u);
                        u = parent[u];
                    }
                    assert!(ancestors.len() == level[i] + 1);
                    for k in 0..=level[i] / 2 {
                        assert!(ladder.query(ancestors[k], k) == ancestors[2 * k]);
                    }
                }
            }
        }
    }
}
