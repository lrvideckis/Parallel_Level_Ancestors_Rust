use rayon::prelude::*;
use rayon::iter::once;
use rayon_scan::ScanParallelIterator;
use paradis_core::{BoundedParAccess, IntoParAccess};

pub struct JumpPointersModedLevel {
    parent: Vec<usize>,
    level: Vec<usize>,
    mod_val: usize,
    mod_value_with_least_nodes: usize,
    special_nodes: Vec<usize>,
    index_in_special_nodes: Vec<usize>,
    jump: Vec<Vec<usize>>,
}

// approach described here:
// https://codeforces.com/blog/entry/100826#comment-1013981
impl JumpPointersModedLevel {
    pub fn new(parent: &[usize], level: &[usize], p: usize) -> Self {
        let n = parent.len();
        assert!(n >= 1 && p >= 1);

        let mod_val = n.div_ceil(p);

        let mut frequency = vec![0; mod_val * p];
        let access = frequency.into_par_access();
        (0..p).into_par_iter().for_each(|i| {
            for j in 0..mod_val {
                let node = i * mod_val + j;
                if node < n {
                    let idx = (level[node] % mod_val) * p + i;
                    unsafe {
                        *access.get_unsync(idx) += 1;
                    }
                }
            }
        });

        let frequency: Vec<usize> = once(0)
            .chain(frequency.into_par_iter())
            .scan(|a, b| *a + *b, 0)
            .collect();

        assert!(frequency[mod_val * p] == n);

        let frequency_count = |mod_value: usize| -> usize {
            frequency[(mod_value + 1) * p] - frequency[mod_value * p]
        };

        let mut mod_value_with_least_nodes = 0;
        for j in 1..mod_val {
            if frequency_count(mod_value_with_least_nodes) > frequency_count(j) {
                mod_value_with_least_nodes = j;
            }
        }

        let index_in_special_nodes: Vec<usize> = level
            .par_iter()
            .map(|&lvl| (lvl % mod_val == mod_value_with_least_nodes) as usize)
            .collect();

        let index_in_special_nodes: Vec<usize> = once(0)
            .chain(index_in_special_nodes.into_par_iter())
            .scan(|a, b| *a + *b, 0)
            .collect();

        let special_nodes_len = index_in_special_nodes[n];

        assert!(special_nodes_len <= p);

        let mut special_nodes = vec![0; special_nodes_len];
        let access = special_nodes.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            if level[i] % mod_val == mod_value_with_least_nodes {
                let idx = index_in_special_nodes[i];
                unsafe {
                    *access.get_unsync(idx) = i;
                }
            }
        });

        let bit_width_p: usize = (p.ilog2() + 1).try_into().unwrap();
        let num_special = special_nodes.len();
        let mut jump = vec![vec![0; num_special]; bit_width_p];

        jump[0].par_iter_mut().enumerate().for_each(|(i, slot)| {
            let mut node = special_nodes[i];
            if level[node] < mod_val {
                *slot = i;
                return;
            }
            for _ in 0..mod_val {
                node = parent[node];
            }
            *slot = index_in_special_nodes[node];
        });

        for i in 1..bit_width_p {
            let (prev_rows, current_rows) = jump.split_at_mut(i);
            let prev_row = &prev_rows[i - 1];
            let curr_row = &mut current_rows[0];
            curr_row.par_iter_mut().enumerate().for_each(|(j, val)| {
                *val = prev_row[prev_row[j]];
            });
        }

        Self {
            parent: parent.to_vec(),
            level: level.to_vec(),
            mod_val,
            mod_value_with_least_nodes,
            special_nodes,
            index_in_special_nodes,
            jump,
        }
    }

    pub fn query(&self, mut v: usize, mut k: usize) -> usize {
        assert!(k <= self.level[v]);
        while k > 0 && self.level[v] % self.mod_val != self.mod_value_with_least_nodes {
            v = self.parent[v];
            k -= 1;
        }
        if k >= self.mod_val {
            v = self.index_in_special_nodes[v];
            let divided_value = k / self.mod_val;
            for bit in 0.. {
                if (1 << bit) > divided_value {
                    break;
                }
                if (divided_value & (1 << bit)) != 0 {
                    v = self.jump[bit][v];
                    k -= (1 << bit) * self.mod_val;
                }
            }
            v = self.special_nodes[v];
        }
        assert!(k < self.mod_val);
        while k > 0 {
            v = self.parent[v];
            k -= 1;
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test_kth_ancestor() {
        for n in 1..=100 {
            for p in 1..=(n + 5) {
                let mut parent = vec![0; n];
                let mut level = vec![0; n];

                for i in 1..n {
                    parent[i] = rand::random_range(0..i);
                    level[i] = 1 + level[parent[i]];
                }

                let structure = JumpPointersModedLevel::new(&parent, &level, p);

                for i in 0..n {
                    let mut kth_parent_naive = i;
                    for k in 0..=level[i] {
                        let ans = structure.query(i, k);
                        assert_eq!(ans, kth_parent_naive);
                        kth_parent_naive = parent[kth_parent_naive];
                    }
                }
            }
        }
    }
}
