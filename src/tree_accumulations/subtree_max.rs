use rayon::prelude::*;
use std::cmp::max;

pub struct RMQ<T, F> {
    t: Vec<Vec<T>>,
    op: F,
}

impl<T: Clone + Send + Sync, F: Fn(&T, &T) -> T + Send + Sync> RMQ<T, F> {
    pub fn new(a: &[T], op: F) -> Self {
        let mut t = vec![a.to_vec(); 1];
        let mut i = 0;
        while (2 << i) <= a.len() {
            let prev = &t[i];
            let shift = 1 << i;
            let layer_len = prev.len() - shift;

            let next_layer: Vec<T> = (0..layer_len)
                .into_par_iter()
                .map(|j| op(&prev[j], &prev[j + shift]))
                .collect();

            t.push(next_layer);
            i += 1;
        }
        Self { t, op }
    }

    pub fn query(&self, range: std::ops::Range<usize>) -> T {
        let lg = range.len().ilog2() as usize;
        (self.op)(&self.t[lg][range.start], &self.t[lg][range.end - (1 << lg)])
    }
}

pub fn calculate_subtree_max<T, F>(
    values: &[T],
    time_in: &[usize],
    time_out: &[usize],
    pre_order: &[usize],
    op: F,
) -> Vec<T>
where
    T: Clone + Send + Sync,
    F: Fn(&T, &T) -> T + Send + Sync + Clone,
{
    let n = values.len();
    let a: Vec<T> = pre_order.iter().map(|&node| values[node].clone()).collect();
    let rmq = RMQ::new(&a, op);

    (0..n)
        .into_par_iter()
        .map(|i| rmq.query(time_in[i]..time_out[i]))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn unit_and_stress_test() {
        for n in 1..=100 {
            for _p in 1..=(n + 5) {
                let mut adjacency_list = vec![vec![]; n];
                let mut parent = vec![0; n];
                for i in 1..n {
                    parent[i] = rand::random_range(0..i);
                    adjacency_list[parent[i]].push(i);
                }

                let values: Vec<i32> = (0..n).map(|_| rand::random_range(0..1_000_000)).collect();

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

                let subtree_max_rmq =
                    calculate_subtree_max(&values, &time_in, &time_out, &pre_order, |&x, &y| {
                        max(x, y)
                    });

                let mut expected = vec![0; n];
                fn compute_naive(
                    node: usize,
                    adjacency_list: &[Vec<usize>],
                    values: &[i32],
                    expected: &mut [i32],
                ) -> i32 {
                    let mut curr = values[node];
                    for &child in &adjacency_list[node] {
                        curr = max(curr, compute_naive(child, adjacency_list, values, expected));
                    }
                    expected[node] = curr;
                    curr
                }
                compute_naive(0, &adjacency_list, &values, &mut expected);

                assert_eq!(subtree_max_rmq, expected);
            }
        }
    }
}
