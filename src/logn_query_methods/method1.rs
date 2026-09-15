use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::iter::once;
use rayon::prelude::*;
use rayon_scan::ScanParallelIterator;

// On Finding Lowest Common Ancestors: Simplification and Parallelization
// Authors: Baruch Schieber and Uzi Vishkin

pub struct Method1 {
    parent: Vec<usize>,
    level: Vec<usize>,
    pub inlabel: Vec<usize>,
    head: Vec<usize>,
    pub ascendant: Vec<usize>,
    path: Vec<usize>,
    head_to_start: Vec<usize>,
}

impl Method1 {
    pub fn new(
        parent: &[usize],
        level: &[usize],
        time_in: &[usize],
        time_out: &[usize],
        pre_order: &[usize],
        et_time_in: &[usize],
        et_time_out: &[usize],
    ) -> Self {
        let n = level.len();
        assert!(n >= 1);
        let inlabel: Vec<usize> = (0..n)
            .into_par_iter()
            .map(|i| time_out[i] & (1usize << (time_in[i] ^ time_out[i]).ilog2()).wrapping_neg())
            .collect();
        let mut head = vec![0; n + 1];
        let access = head.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            if i == parent[i] || inlabel[parent[i]] != inlabel[i] {
                unsafe {
                    *access.get_unsync(inlabel[i]) = i;
                }
            }
        });
        let mut head_to_start = vec![0; n + 1];
        let access = head_to_start.into_par_access();
        (0..n).into_par_iter().for_each(|v| {
            if pre_order[inlabel[v] - 1] == v {
                let head_v = head[inlabel[v]];
                let val = level[v] - level[head_v] + 1;
                unsafe {
                    *access.get_unsync(head_v) = val;
                }
            }
        });
        let head_to_start: Vec<usize> = once(0)
            .chain(head_to_start.par_iter().cloned())
            .scan(|a, b| *a + *b, 0)
            .collect();
        assert_eq!(head_to_start[n], n);
        let mut path = vec![0; n];
        let access = path.into_par_access();
        (0..n).into_par_iter().for_each(|v| {
            let head_v = head[inlabel[v]];
            let idx = head_to_start[head_v] + level[v] - level[head_v];
            unsafe {
                *access.get_unsync(idx) = v;
            }
        });
        let mut diff = vec![0isize; 2 * n];
        let access_diff = diff.into_par_access();
        (0..n).into_par_iter().for_each(|i| {
            if i == parent[i] || inlabel[parent[i]] != inlabel[i] {
                let lsb_val = (inlabel[i] & inlabel[i].wrapping_neg()) as isize;
                unsafe {
                    *access_diff.get_unsync(et_time_in[i]) += lsb_val;
                    *access_diff.get_unsync(et_time_out[i]) -= lsb_val;
                }
            }
        });
        let diff: Vec<isize> = diff.into_par_iter().scan(|a, b| *a + *b, 0).collect();
        let ascendant: Vec<usize> = (0..n)
            .into_par_iter()
            .map(|v| diff[et_time_in[v]] as usize)
            .collect();
        Self {
            parent: parent.to_vec(),
            level: level.to_vec(),
            inlabel,
            head,
            path,
            head_to_start,
            ascendant,
        }
    }

    pub fn kth_parent(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        let mut u = v;
        loop {
            let head_u = self.head[self.inlabel[u]];
            if self.level[v] - k >= self.level[head_u] {
                return self.path
                    [self.head_to_start[head_u] + self.level[v] - k - self.level[head_u]];
            }
            u = self.parent[head_u];
        }
    }

    // O(log(log(n))) but with bad constant factor, mostly just to show it's possible
    pub fn kth_parent_binary_search_paths(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        let anc_d = self.level[v] - k;
        let n = self.level.len();
        let mut start: i32 = -1;
        let mut end: i32 = n.ilog2() as i32;
        while start + 1 < end {
            let mid = (start + end) / 2;
            let b = (self.ascendant[v] & (1usize << mid).wrapping_neg()).isolate_lowest_one();
            let curr_inlabel = (self.inlabel[v] & b.wrapping_neg()) | b;
            let curr_head = self.head[curr_inlabel];
            let dist_to_go = anc_d as isize - self.level[curr_head] as isize;
            if dist_to_go >= 0 {
                end = mid;
            } else {
                start = mid;
            }
        }
        let b = (self.ascendant[v] & (1usize << end).wrapping_neg()).isolate_lowest_one();
        let curr_inlabel = (self.inlabel[v] & b.wrapping_neg()) | b;
        let curr_head = self.head[curr_inlabel];
        let dist_to_go = anc_d - self.level[curr_head];
        self.path[self.head_to_start[curr_head] + dist_to_go]
    }

    pub fn lowest_common_ancestor(&self, mut u: usize, mut v: usize) -> usize {
        let j = self.inlabel[u] ^ self.inlabel[v];
        if j != 0 {
            let j = self.ascendant[u] & self.ascendant[v] & (1usize << j.ilog2()).wrapping_neg();
            let k = self.ascendant[u] ^ j;
            if k != 0 {
                let k = 1usize << k.ilog2();
                u = self.parent[self.head[(self.inlabel[u] & k.wrapping_neg()) | k]];
            }
            let k = self.ascendant[v] ^ j;
            if k != 0 {
                let k = 1usize << k.ilog2();
                v = self.parent[self.head[(self.inlabel[v] & k.wrapping_neg()) | k]];
            }
        }
        if self.level[u] < self.level[v] { u } else { v }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=300 {
            let mut adjacency_list = vec![vec![]; n];
            let mut parent = vec![0; n];
            let mut level = vec![0; n];
            for i in 1..n {
                parent[i] = rand::random_range(0..i);
                level[i] = 1 + level[parent[i]];
                adjacency_list[parent[i]].push(i);
            }
            let mut pre_order = vec![0; n];
            let mut time_in = vec![0; n];
            let mut time_out = vec![0; n];
            let mut et_time_in = vec![0; n];
            let mut et_time_out = vec![0; n];
            let mut timer = 0;
            let mut et_timer = 0;
            fn dfs(
                v: usize,
                adj: &[Vec<usize>],
                time_in: &mut [usize],
                time_out: &mut [usize],
                pre_order: &mut [usize],
                timer: &mut usize,
                et_time_in: &mut [usize],
                et_time_out: &mut [usize],
                et_timer: &mut usize,
            ) {
                time_in[v] = *timer;
                pre_order[*timer] = v;
                *timer += 1;

                et_time_in[v] = *et_timer;
                *et_timer += 1;

                for &child in &adj[v] {
                    dfs(
                        child,
                        adj,
                        time_in,
                        time_out,
                        pre_order,
                        timer,
                        et_time_in,
                        et_time_out,
                        et_timer,
                    );
                }
                time_out[v] = *timer;

                et_time_out[v] = *et_timer;
                *et_timer += 1;
            }

            dfs(
                0,
                &adjacency_list,
                &mut time_in,
                &mut time_out,
                &mut pre_order,
                &mut timer,
                &mut et_time_in,
                &mut et_time_out,
                &mut et_timer,
            );
            assert_eq!(timer, n);
            assert_eq!(et_timer, 2 * n);

            let ancestor = Method1::new(
                &parent,
                &level,
                &time_in,
                &time_out,
                &pre_order,
                &et_time_in,
                &et_time_out,
            );

            for i in 0..n {
                let mut kth_parent_naive = i;
                for k in 0..=level[i] {
                    assert_eq!(kth_parent_naive, ancestor.kth_parent(i, k));
                    assert_eq!(
                        kth_parent_naive,
                        ancestor.kth_parent_binary_search_paths(i, k)
                    );
                    kth_parent_naive = parent[kth_parent_naive];
                }
            }

            let naive_lca = |mut u: usize, mut v: usize| -> usize {
                while u != v {
                    if level[u] > level[v] {
                        u = parent[u];
                    } else {
                        v = parent[v];
                    }
                }
                u
            };

            for u in 0..n {
                for v in 0..n {
                    assert_eq!(naive_lca(u, v), ancestor.lowest_common_ancestor(u, v));
                }
            }
        }
    }
}
