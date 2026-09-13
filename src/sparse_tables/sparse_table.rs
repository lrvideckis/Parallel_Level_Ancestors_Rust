use rayon::prelude::*;

pub struct SparseTable<T, F> {
    t: Vec<Vec<T>>,
    op: F,
}

impl<T: Clone + Send + Sync, F: Fn(&T, &T) -> T + Send + Sync> SparseTable<T, F> {
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

    #[test]
    fn stress_test() {
        for n in 0..500 {
            let a: Vec<i32> = (0..n).map(|_| rand::random_range(0..1_000_000)).collect();

            let sparse_table = SparseTable::new(&a, |&x, &y| std::cmp::min(x, y));

            let n = a.len();
            for i in 0..n {
                let mut naive = a[i];
                for j in i..n {
                    naive = std::cmp::min(naive, a[j]);
                    assert_eq!(naive, sparse_table.query(i..j + 1));
                }
            }
        }
    }
}
