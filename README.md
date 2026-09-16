An implementation of https://link.springer.com/chapter/10.1007/978-3-032-29003-8_29 in rust.

The purpose is to help in understanding and to show correctness of the paper, not to be the most optimized implementation. So this implementation follows the details in the paper as closely as possible.

---

Note about parallel complexity: the paper assumes PRAM model where something like `(0..n).into_par_iter()` is O(n) work, O(1) span. But rayon crate implements this in O(n) work, O(logn) span due to recursive tree splitting. I'm pretty sure this implementation is actually O(log^2 n) span due to this, so it is a bit fudged. Let's just sweep this under the rug.

---

Note about parallel list ranking: The paper assumes the tree is given as an edge list where for every node u, the set of edges (u <-> some child of u) form a subarray in the edge list. Then parallel list ranking can calculate:

- a pre order traversal
- arrays time_in, time_out such that subarray [time_in[u], time_out[u]) of pre_order corresponds to u's subtree
- level a.k.a. depth
- euler tour
- arrays time_in, time_out such that subarray [et_time_in[u], et_time_out[u]) of euler_tour corresponds to u's subtree

I chose not to implement parallel list ranking as it seems rather complicated. Instead these arrays are calculated naively in the tests as needed.
 
---

Note about use of unsafe: Consider the following problem: a permutation of 0..=n-1 is stored in p. Calculate the inverse permutation in p_inv. Sequentially it is trivial:

```
for i in 0..n {
    p_inv[p[i]] = i;
}
```

Now how to do this in parallel? We know there are no race conditions when writing to p_inv, but the borrow checker doesn't know this. So we need unsafe. I decided on the paradis crate:

```
let mut p_inv = vec![0, n];
let access = p_inv.into_par_access();
(0..n).into_par_iter().for_each(|i| {
    unsafe {
        *access.get_unsync(p[i]) = i;
    }
});
```
There are some places in the paper where we can prove there are no race conditions, but I couldn't think of how to implement it in safe rust.
