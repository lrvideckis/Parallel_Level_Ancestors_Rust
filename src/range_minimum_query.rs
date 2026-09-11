use rayon::prelude::*;

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn unit_test() {
        let a = [1, 3, 2, 4, 5];
        let rmq = RMQ::new(&a, |&x, &y| std::cmp::min(x, y));
        let n = a.len();
        for i in 0..n {
            let mut naive = a[i];
            for j in i..n {
                naive = std::cmp::min(naive, a[j]);
                assert_eq!(naive, rmq.query(i..j + 1));
            }
        }
    }

    #[test]
    fn stress_test() {
        for n in 0..500 {
            let a: Vec<i32> = (0..n).map(|_| rand::random_range(0..1_000_000)).collect();

            let rmq = RMQ::new(&a, |&x, &y| std::cmp::min(x, y));

            let n = a.len();
            for i in 0..n {
                let mut naive = a[i];
                for j in i..n {
                    naive = std::cmp::min(naive, a[j]);
                    assert_eq!(naive, rmq.query(i..j + 1));
                }
            }
        }
    }

    #[test]
    fn benchmark_rmq_build() {
        let n = 1_000_000;
        let a: Vec<i32> = (0..n).map(|_| rand::random_range(0..1_000_000)).collect();

        let start = Instant::now();
        let _rmq_multi = RMQ::new(&a, |&x, &y| std::cmp::min(x, y));
        let multi_duration = start.elapsed();

        let single_thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();

        let start = Instant::now();
        let _rmq_single = single_thread_pool.install(|| {
            RMQ::new(&a, |&x, &y| std::cmp::min(x, y))
        });
        let single_duration = start.elapsed();

        println!("\n--- RMQ Build Benchmark (N = {}) ---", n);
        println!("Multi-threaded time:  {:?}", multi_duration);
        println!("Single-threaded time: {:?}", single_duration);

        let speedup = single_duration.as_secs_f64() / multi_duration.as_secs_f64();
        println!("Speedup: {:.2}x\n", speedup);
    }
}
