use crate::helpers::ladders::Ladders;
use rayon::prelude::*;

pub struct Variation2 {
    level: Vec<usize>,
    et_time_in: Vec<usize>,
    jump: Vec<Vec<usize>>,
    ladders: Ladders,
}

impl Variation2 {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        parent: &[usize],
        level: &[usize],
        time_in: &[usize],
        time_out: &[usize],
        pre_order: &[usize],
        euler_tour: &[usize],
        et_time_in: &[usize],
        p: usize,
    ) -> Self {
        let n = parent.len();
        assert!(n >= 1 && p >= 1);

        let ladders = Ladders::new(parent, level, time_in, time_out, pre_order, p, |len| {
            2 * len
        });

        let block_size = (2 * n).div_ceil(p);
        let mut jump: Vec<Vec<usize>> = vec![vec![]; 2 * n];
        jump.par_chunks_mut(block_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                for (j, item) in chunk.iter_mut().enumerate() {
                    let original_idx = i * block_size + j;
                    if original_idx == 0 {
                        continue;
                    }
                    let mut u = parent[euler_tour[original_idx]];
                    item.push(u);
                    let mut k = 1;
                    while k < original_idx.isolate_lowest_one() {
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
            2 * (2 * n - 1) - ((2 * n - 1).count_ones() as usize)
        );

        Self {
            level: level.to_vec(),
            et_time_in: et_time_in.to_vec(),
            jump,
            ladders,
        }
    }

    pub fn kth_parent(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        if k == 0 {
            v
        } else {
            let i = self.et_time_in[v];
            let j = (i + k - 1) & k.next_power_of_two().wrapping_neg();
            assert!(i.abs_diff(j) < k);
            let dist_to_go = self.level[self.jump[j][0]] - (self.level[v] - k) + 1;
            let jump_node = self.jump[j][dist_to_go.ilog2() as usize];
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
                let mut et_time_in = vec![0; 2 * n];

                {
                    let mut timer = 0;
                    let mut timer_euler_tour = 1;
                    fn dfs(
                        node: usize,
                        timer: &mut usize,
                        timer_euler_tour: &mut usize,
                        time_in: &mut [usize],
                        time_out: &mut [usize],
                        pre_order: &mut [usize],
                        euler_tour: &mut [usize],
                        et_time_in: &mut [usize],
                        adjacency_list: &[Vec<usize>],
                    ) {
                        euler_tour[*timer_euler_tour] = node;
                        et_time_in[node] = *timer_euler_tour;
                        *timer_euler_tour += 1;
                        time_in[node] = *timer;
                        pre_order[*timer] = node;
                        *timer += 1;
                        for &child in &adjacency_list[node] {
                            dfs(
                                child,
                                timer,
                                timer_euler_tour,
                                time_in,
                                time_out,
                                pre_order,
                                euler_tour,
                                et_time_in,
                                adjacency_list,
                            );
                            euler_tour[*timer_euler_tour] = node;
                            *timer_euler_tour += 1;
                        }
                        time_out[node] = *timer;
                    }
                    dfs(
                        0,
                        &mut timer,
                        &mut timer_euler_tour,
                        &mut time_in,
                        &mut time_out,
                        &mut pre_order,
                        &mut euler_tour,
                        &mut et_time_in,
                        &adjacency_list,
                    );
                }

                let variation2 = Variation2::new(
                    &parent,
                    &level,
                    &time_in,
                    &time_out,
                    &pre_order,
                    &euler_tour,
                    &et_time_in,
                    p,
                );
                for i in 0..n {
                    let mut kth_parent_naive = i;
                    for k in 0..=level[i] {
                        let ans = variation2.kth_parent(i, k);
                        assert_eq!(ans, kth_parent_naive);
                        kth_parent_naive = parent[kth_parent_naive];
                    }
                }
            }
        }
    }
}
