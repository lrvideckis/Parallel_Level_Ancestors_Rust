use crate::helpers::ladders::Ladders;
use paradis_core::{BoundedParAccess, IntoParAccess};
use rayon::prelude::*;

pub struct Hagerup {
    level: Vec<usize>,
    euler_tour_index: Vec<usize>,
    jump: Vec<usize>,
    ladders: Ladders,
    kappa: usize,
}

impl Hagerup {
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
        kappa: usize,
    ) -> Self {
        let n = parent.len();
        assert!(n >= 1 && kappa >= 1 && p >= 1);

        let ladders = Ladders::new(parent, level, time_in, time_out, pre_order, p, |len| {
            ((len * (2 * kappa + 3)) / kappa).max(len + kappa)
        });

        let mut jump = vec![0; 2 * n];
        let access = jump.into_par_access();

        (1..2 * n).into_par_iter().step_by(2).for_each(|i| {
            let node = euler_tour[i];
            let jump_node = ladders.kth_parent(node, kappa.min(level[node]));
            unsafe {
                *access.get_unsync(i) = jump_node;
            }
        });

        let mut i = 2;
        while i < 2 * n {
            (i..2 * n).into_par_iter().step_by(2 * i).for_each(|j| {
                let prev_jump = unsafe { *access.get_unsync(j - i / 2) };
                let node = euler_tour[j];
                let step = kappa * i;
                let jump_node = if level[node] <= step {
                    0
                } else {
                    let target_level = level[node] - step;
                    let dist = level[prev_jump] - target_level;
                    ladders.kth_parent(prev_jump, dist)
                };
                unsafe {
                    *access.get_unsync(j) = jump_node;
                }
            });
            i *= 2;
        }

        Self {
            level: level.to_vec(),
            euler_tour_index: euler_tour_index.to_vec(),
            jump,
            ladders,
            kappa,
        }
    }

    pub fn kth_parent(&self, mut v: usize, k: usize) -> usize {
        assert!(k <= self.level[v]);
        let anc_d = self.level[v] - k;
        if k >= (self.kappa + 1) {
            let l = 1_usize << (k / (self.kappa + 1)).ilog2();
            v = self.jump[(self.euler_tour_index[v] & l.wrapping_neg()) | l];
        }
        self.ladders.kth_parent(v, self.level[v] - anc_d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_test() {
        for n in 1..=50 {
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

            for kappa in 1..=7 {
                for p in 1..=(n + 5) {
                    let hagerup = Hagerup::new(
                        &parent,
                        &level,
                        &time_in,
                        &time_out,
                        &pre_order,
                        &euler_tour,
                        &euler_tour_index,
                        p,
                        kappa,
                    );
                    for i in 0..n {
                        let mut kth_parent_naive = i;
                        for k in 0..=level[i] {
                            assert_eq!(hagerup.kth_parent(i, k), kth_parent_naive);
                            kth_parent_naive = parent[kth_parent_naive];
                        }
                    }
                }
            }
        }
    }
}
