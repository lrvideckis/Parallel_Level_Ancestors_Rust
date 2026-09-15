use crate::helpers::sparse_table::SparseTable;
use rayon::prelude::*;

pub struct Method2 {
    depth: Vec<usize>,
    euler_tour: Vec<usize>,
    time_in: Vec<usize>,
    b: usize,
    block_min_level: Vec<usize>,
    block_table: Vec<Vec<Vec<usize>>>,
    prefix_of_block: Vec<usize>,
    suffix_of_block: Vec<usize>,
    sparse_table: SparseTable<usize, fn(&usize, &usize) -> usize>,
}

impl Method2 {
    pub fn new(depth: &[usize], euler_tour: &[usize], time_in: &[usize], p: usize) -> Self {
        let m = euler_tour.len();
        assert!(m >= 1);
        assert!(p >= 1);

        let b = m.div_ceil(p);

        let mut block_min_level = vec![0; p];
        let mut block_table = vec![vec![]; p];

        block_min_level
            .par_iter_mut()
            .zip(block_table.par_iter_mut())
            .enumerate()
            .for_each(|(i, (min_lvl, table))| {
                let start_idx = i * b;
                if start_idx >= m {
                    return;
                }
                let end_idx = (start_idx + b).min(m);

                let mut local_min = usize::MAX;
                let mut local_max = 0;
                for j in start_idx..end_idx {
                    let d = depth[euler_tour[j]];
                    local_min = local_min.min(d);
                    local_max = local_max.max(d);
                }
                *min_lvl = local_min;

                if local_min != usize::MAX {
                    let range_d = local_max - local_min + 1;
                    let mut local_table = vec![vec![]; range_d];
                    for &node in &euler_tour[start_idx..end_idx] {
                        let lev = depth[node] - local_min;
                        local_table[lev].push(node);
                    }
                    *table = local_table;
                }
            });

        let mut prefix_of_block: Vec<usize> = euler_tour.par_iter().map(|&v| depth[v]).collect();
        let mut suffix_of_block = prefix_of_block.clone();

        prefix_of_block.par_chunks_mut(b).for_each(|chunk| {
            for j in 1..chunk.len() {
                chunk[j] = chunk[j].min(chunk[j - 1]);
            }
        });

        suffix_of_block.par_chunks_mut(b).for_each(|chunk| {
            for j in (1..chunk.len()).rev() {
                chunk[j - 1] = chunk[j - 1].min(chunk[j]);
            }
        });

        let num_blocks = m.div_ceil(b);
        let block_suffs: Vec<usize> = (0..num_blocks).map(|i| suffix_of_block[i * b]).collect();

        let min_fn: fn(&usize, &usize) -> usize = |&x, &y| x.min(y);
        let sparse_table = SparseTable::new(block_suffs, min_fn);

        Self {
            depth: depth.to_vec(),
            euler_tour: euler_tour.to_vec(),
            time_in: time_in.to_vec(),
            b,
            block_min_level,
            block_table,
            prefix_of_block,
            suffix_of_block,
            sparse_table,
        }
    }

    fn query(&self, l: usize, r: usize) -> usize {
        assert!(l / self.b < (r - 1) / self.b);
        let mut res = self.suffix_of_block[l].min(self.prefix_of_block[r - 1]);
        if l / self.b + 1 < r / self.b {
            res = res.min(self.sparse_table.query((l / self.b + 1)..(r / self.b)));
        }
        res
    }

    pub fn kth_parent(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.depth[v]);
        let anc_d = self.depth[v] - k;
        let tv = self.time_in[v];

        if self.prefix_of_block[tv] <= anc_d {
            let block_id = tv / self.b;
            let lev = anc_d - self.block_min_level[block_id];
            let nodes_on_level = &self.block_table[block_id][lev];

            let mut start = 0;
            let mut end = nodes_on_level.len();
            while start + 1 < end {
                let mid = (start + end) / 2;
                if self.time_in[nodes_on_level[mid]] <= tv {
                    start = mid;
                } else {
                    end = mid;
                }
            }
            return nodes_on_level[start];
        }

        let mut start = 0;
        let mut end = (tv / self.b) * self.b;
        while start + 1 < end {
            let mid = (start + end) / 2;
            if self.query(mid, tv + 1) <= anc_d {
                start = mid;
            } else {
                end = mid;
            }
        }
        self.euler_tour[start]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=80 {
            for p in 1..=(2 * n + 5) {
                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i);
                    adjacency_list[parent[i]].push(i);
                }

                let mut depth = vec![0; n];
                let mut euler_tour = vec![0; 2 * n - 1];
                let mut time_in = vec![0; n];
                let mut time_out = vec![0; n];
                let mut timer = 0;

                fn dfs(
                    v: usize,
                    adj: &[Vec<usize>],
                    depth: &mut [usize],
                    euler_tour: &mut [usize],
                    time_in: &mut [usize],
                    time_out: &mut [usize],
                    timer: &mut usize,
                ) {
                    time_in[v] = *timer;
                    euler_tour[*timer] = v;
                    *timer += 1;

                    for &child in &adj[v] {
                        depth[child] = 1 + depth[v];
                        dfs(child, adj, depth, euler_tour, time_in, time_out, timer);
                        euler_tour[*timer] = v;
                        *timer += 1;
                    }
                    time_out[v] = *timer;
                }

                dfs(
                    0,
                    &adjacency_list,
                    &mut depth,
                    &mut euler_tour,
                    &mut time_in,
                    &mut time_out,
                    &mut timer,
                );
                assert_eq!(timer, 2 * n - 1);

                let ancestor = Method2::new(&depth, &euler_tour, &time_in, p);

                for i in 0..n {
                    let mut kth_parent_naive = i;
                    for k in 0..=depth[i] {
                        assert_eq!(kth_parent_naive, ancestor.kth_parent(i, k));
                        kth_parent_naive = parent[kth_parent_naive];
                    }
                }
            }
        }
    }
}
