use rayon::prelude::*;

pub struct Method1 {
    parent: Vec<usize>,
    level: Vec<usize>,
    inlabel: Vec<usize>,
    head: Vec<usize>,
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
    ) -> Self {
        let n = level.len();
        assert!(n >= 1);
        let inlabel: Vec<usize> = (0..n)
            .into_par_iter()
            .map(|i| time_out[i] & (1usize << (time_in[i] ^ time_out[i]).ilog2()).wrapping_neg())
            .collect();
        let mut head = vec![0; n + 1];
        for i in 0..n {
            if i == parent[i] || inlabel[parent[i]] != inlabel[i] {
                head[inlabel[i]] = i;
            }
        }
        let mut head_to_start = vec![0; n + 1];
        for v in 0..n {
            if pre_order[inlabel[v]] == v {
                let head_v = head[inlabel[v]];
                head_to_start[head_v] = level[v] - level[head_v] + 1;
            }
        }
        let mut sum = 0;
        for i in 0..=n {
            let val = head_to_start[i];
            head_to_start[i] = sum;
            sum += val;
        }
        assert_eq!(sum, n);
        let mut path = vec![0; n];
        for v in 0..n {
            let head_v = head[inlabel[v]];
            let idx = head_to_start[head_v] + level[v] - level[head_v];
            path[idx] = v;
        }
        Self {
            parent: parent.to_vec(),
            level: level.to_vec(),
            inlabel,
            head,
            path,
            head_to_start,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=1000 {
            let mut adjacency_list = vec![vec![]; n];
            let mut parent = vec![0; n];
            let mut level = vec![0; n];
            for i in 1..n {
                parent[i] = rand::random_range(0..i);
                level[i] = 1 + level[parent[i]];
                adjacency_list[parent[i]].push(i);
            }
            let mut pre_order = vec![0; n + 1];
            let mut time_in = vec![0; n];
            let mut time_out = vec![0; n];
            let mut timer = 0;
            fn dfs(
                v: usize,
                adj: &[Vec<usize>],
                time_in: &mut [usize],
                time_out: &mut [usize],
                pre_order: &mut [usize],
                timer: &mut usize,
            ) {
                time_in[v] = *timer;
                *timer += 1;
                pre_order[*timer] = v;
                for &child in &adj[v] {
                    dfs(child, adj, time_in, time_out, pre_order, timer);
                }
                time_out[v] = *timer;
            }
            dfs(
                0,
                &adjacency_list,
                &mut time_in,
                &mut time_out,
                &mut pre_order,
                &mut timer,
            );
            assert_eq!(timer, n);
            let ancestor = Method1::new(&parent, &level, &time_in, &time_out, &pre_order);
            for i in 0..n {
                let mut kth_parent_naive = i;
                for k in 0..=level[i] {
                    assert_eq!(kth_parent_naive, ancestor.kth_parent(i, k));
                    kth_parent_naive = parent[kth_parent_naive];
                }
            }
        }
    }
}
