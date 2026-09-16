use crate::helpers::ladders::Ladders;
use crate::logn_query_methods::method1::Method1; // Adjust module path as needed
use rayon::prelude::*;

pub struct Variation3 {
    level: Vec<usize>,
    method1: Method1,
    jump: Vec<Vec<usize>>,
    ladders: Ladders,
}

impl Variation3 {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        parent: &[usize],
        level: &[usize],
        time_in: &[usize],
        time_out: &[usize],
        pre_order: &[usize],
        et_time_in: &[usize],
        et_time_out: &[usize],
        p: usize,
    ) -> Self {
        let n = parent.len();
        assert!(n >= 1 && p >= 1);

        let method1 = Method1::new(
            parent,
            level,
            time_in,
            time_out,
            pre_order,
            et_time_in,
            et_time_out,
        );

        let ladders = Ladders::new(parent, level, time_in, time_out, pre_order, p, |len| {
            2 * len
        });

        let block_size = n.div_ceil(p);
        let mut jump: Vec<Vec<usize>> = vec![vec![]; n];

        jump.par_chunks_mut(block_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                for (j, item) in chunk.iter_mut().enumerate() {
                    let original_idx = i * block_size + j;
                    let mut u = parent[pre_order[original_idx]];
                    item.push(u);
                    let mut k = 1;
                    while k <= (original_idx + 1).isolate_lowest_one() {
                        u = ladders.kth_parent(u, std::cmp::min(k, level[u]));
                        //push even when u goes above root so that we can verify total number of jump
                        //pointers
                        item.push(u);
                        k *= 2;
                    }
                }
            });

        assert_eq!(
            jump.iter().map(|v| v.len()).sum::<usize>(),
            3 * n - (n.count_ones() as usize)
        );

        Self {
            level: level.to_vec(),
            method1,
            jump,
            ladders,
        }
    }

    pub fn kth_parent(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        if k == 0 {
            v
        } else {
            let i = k.ilog2();
            let j = (self.method1.ascendant[v] & (1usize << i).wrapping_neg()).trailing_zeros();
            let la_bt = (self.method1.inlabel[v] & (1usize << j).wrapping_neg()) | (1 << j);
            let dist_to_go = self.level[self.jump[la_bt - 1][0]] - (self.level[v] - k) + 1;
            let jump_node = self.jump[la_bt - 1][dist_to_go.ilog2() as usize];
            self.ladders
                .kth_parent(jump_node, self.level[jump_node] - (self.level[v] - k))
        }
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
                let mut euler_tour = vec![0; 2 * n];
                let mut et_time_in = vec![0; n];
                let mut et_time_out = vec![0; n];

                {
                    let mut timer = 0;
                    let mut et_timer = 0;
                    let mut timer_euler_tour = 1;
                    fn dfs(
                        node: usize,
                        timer: &mut usize,
                        et_timer: &mut usize,
                        timer_euler_tour: &mut usize,
                        time_in: &mut [usize],
                        time_out: &mut [usize],
                        pre_order: &mut [usize],
                        euler_tour: &mut [usize],
                        et_time_in: &mut [usize],
                        et_time_out: &mut [usize],
                        adjacency_list: &[Vec<usize>],
                    ) {
                        euler_tour[*timer_euler_tour] = node;
                        *timer_euler_tour += 1;

                        time_in[node] = *timer;
                        pre_order[*timer] = node;
                        *timer += 1;

                        et_time_in[node] = *et_timer;
                        *et_timer += 1;

                        for &child in &adjacency_list[node] {
                            dfs(
                                child,
                                timer,
                                et_timer,
                                timer_euler_tour,
                                time_in,
                                time_out,
                                pre_order,
                                euler_tour,
                                et_time_in,
                                et_time_out,
                                adjacency_list,
                            );
                            euler_tour[*timer_euler_tour] = node;
                            *timer_euler_tour += 1;
                        }
                        time_out[node] = *timer;

                        et_time_out[node] = *et_timer;
                        *et_timer += 1;
                    }
                    dfs(
                        0,
                        &mut timer,
                        &mut et_timer,
                        &mut timer_euler_tour,
                        &mut time_in,
                        &mut time_out,
                        &mut pre_order,
                        &mut euler_tour,
                        &mut et_time_in,
                        &mut et_time_out,
                        &adjacency_list,
                    );
                }

                let variation3 = Variation3::new(
                    &parent,
                    &level,
                    &time_in,
                    &time_out,
                    &pre_order,
                    &et_time_in,
                    &et_time_out,
                    p,
                );

                for i in 0..n {
                    let mut kth_parent_naive = i;
                    for k in 0..=level[i] {
                        let ans = variation3.kth_parent(i, k);
                        assert_eq!(ans, kth_parent_naive);
                        kth_parent_naive = parent[kth_parent_naive];
                    }
                }
            }
        }
    }
}
