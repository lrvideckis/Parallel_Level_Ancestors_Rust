use crate::ladders::Ladders;
//use rayon::prelude::*;

pub struct Variation1 {
    level: Vec<usize>,
    euler_tour_index: Vec<usize>,
    jump: Vec<[usize; 2]>,
    ladders: Ladders,
}

impl Variation1 {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        parent: &[usize],
        level: &[usize],
        time_in: &[usize],
        time_out: &[usize],
        pre_order: &[usize],
        euler_tour: &[usize],
        euler_tour_index: &[usize],
        p: usize,
    ) -> Self {
        let n = parent.len();
        assert!(n >= 1 && p >= 1);

        let ladders = Ladders::new(parent, level, time_in, time_out, pre_order, p);

        let mut jump: Vec<[usize; 2]> = vec![[0, 0]; 2 * n];

        let mut i = 1;
        while i < 2 * n {
            let u = if level[euler_tour[i - 1]] > level[euler_tour[i]] {
                euler_tour[i - 1]
            } else {
                euler_tour[i]
            };
            jump[i] = [parent[u], parent[parent[u]]];
            i += 2;
        }

        for rnz in 1.. {
            let mut i = 1 << rnz;
            if i >= 2 * n {
                break;
            }
            while i < 2 * n {
                let l_ch = i - (1 << (rnz - 1));
                let r_ch = i + (1 << (rnz - 1));
                for j in 0..=1 {
                    let mut node = jump[l_ch][j];
                    if r_ch < 2 * n && level[node] < level[jump[r_ch][j]] {
                        node = jump[r_ch][j];
                    }
                    jump[i][j] =
                        ladders.query(node, std::cmp::min(1 << (rnz - 1 + j), level[node]));
                }
                i += 2 << rnz;
            }
        }

        Self {
            level: level.to_vec(),
            euler_tour_index: euler_tour_index.to_vec(),
            jump,
            ladders,
        }
    }

    pub fn query(&self, v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        if k == 0 {
            v
        } else {
            let i = self.euler_tour_index[v];
            let j = k.ilog2();
            let la_bt = (i & (1_usize << j).wrapping_neg()) | (1 << j);
            assert!(i.abs_diff(la_bt) <= k);
            let jump_node = if self.level[self.jump[la_bt][1]] >= self.level[v] - k {
                self.jump[la_bt][1]
            } else {
                self.jump[la_bt][0]
            };
            self.ladders
                .query(jump_node, self.level[jump_node] - (self.level[v] - k))
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
                        adjacency_list: &[Vec<usize>],
                    ) {
                        euler_tour[*timer_euler_tour] = node;
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
                        &adjacency_list,
                    );
                }

                let mut euler_tour_index = vec![0; n];
                for i in 1..2 * n {
                    euler_tour_index[euler_tour[i]] = i;
                }

                let variation1 = Variation1::new(
                    &parent,
                    &level,
                    &time_in,
                    &time_out,
                    &pre_order,
                    &euler_tour,
                    &euler_tour_index,
                    p,
                );
                for i in 0..n {
                    let mut kth_parent_naive = i;
                    for k in 0..=level[i] {
                        let ans = variation1.query(i, k);
                        assert_eq!(ans, kth_parent_naive);
                        kth_parent_naive = parent[kth_parent_naive];
                    }
                }
            }
        }
    }
}
