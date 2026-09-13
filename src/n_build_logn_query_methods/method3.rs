use rayon::prelude::*;

pub struct Method3 {
    level: Vec<usize>,
    pre_order: Vec<usize>,
    time_in: Vec<usize>,
    seg_tree: Vec<usize>,
}

impl Method3 {
    pub fn new(level: &[usize], pre_order: &[usize], time_in: &[usize]) -> Self {
        let n = level.len();

        let mut seg_tree = vec![0; 2 * n];
        seg_tree[n..2 * n]
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, val)| {
                *val = level[pre_order[i]];
            });

        let mut level_size = if n > 0 { 1 << (n as u32).ilog2() } else { 0 };
        let mut prev = n;

        while level_size >= 1 {
            let updates: Vec<(usize, usize)> = (level_size..prev)
                .into_par_iter()
                .map(|i| (i, std::cmp::min(seg_tree[2 * i], seg_tree[2 * i + 1])))
                .collect();
            for (i, val) in updates {
                seg_tree[i] = val;
            }
            prev = level_size;
            level_size /= 2;
        }

        Self {
            level: level.to_vec(),
            pre_order: pre_order.to_vec(),
            time_in: time_in.to_vec(),
            seg_tree,
        }
    }

    // seg tree walk in this style: https://codeforces.com/blog/entry/118682
    pub fn query(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        let mut l = 1;
        let mut r = self.time_in[v] + 1;
        while l < r {
            let u = r + self.level.len();
            let b = std::cmp::min(u.isolate_lowest_one(), r - l).ilog2() as usize;
            let m = r - (1 << b);
            if self.seg_tree[(u - 1) >> b] > self.level[v] - k {
                r = m;
            } else {
                l = m + 1;
            }
        }
        self.pre_order[l - 1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=2000 {
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
            let mut timer = 0;

            fn dfs(
                v: usize,
                adj: &[Vec<usize>],
                time_in: &mut [usize],
                pre_order: &mut [usize],
                timer: &mut usize,
            ) {
                time_in[v] = *timer;
                pre_order[*timer] = v;
                *timer += 1;
                for &child in &adj[v] {
                    dfs(child, adj, time_in, pre_order, timer);
                }
            }

            dfs(0, &adjacency_list, &mut time_in, &mut pre_order, &mut timer);
            assert_eq!(timer, n);

            let ancestor = Method3::new(&level, &pre_order, &time_in);

            for i in 0..n {
                let mut kth_parent_naive = i;
                for k in 0..=level[i] {
                    assert_eq!(kth_parent_naive, ancestor.query(i, k));
                    kth_parent_naive = parent[kth_parent_naive];
                }
            }
        }
    }
}
