use crate::ruranges_structs::Interval;

pub fn nearest_next_nonoverlapping_k(
    is1: &mut [Interval], // sorted by .end
    is2: &mut [Interval], // sorted by .start
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // Sort is1 by end, is2 by start for the "next" intervals
    radsort::sort_by_key(is1, |i| i.end);
    radsort::sort_by_key(is2, |i| i.start);

    let len_l = is1.len();
    let len_r = is2.len();

    // We'll store up to k results per interval
    let mut results_idx1 = vec![-1; k * len_l];
    let mut results_idx2 = vec![-1; k * len_l];
    let mut results_dist = vec![i64::MAX; k * len_l];

    let mut j = 0; // pointer into is2

    for i in 0..len_l {
        let x = is1[i].end;

        // original_i = the interval's original index in the unsorted array
        let original_i = is1[i].idx as usize;

        // Advance j until is2[j].start >= x
        while j < len_r && is2[j].start < x {
            j += 1;
        }

        // We only want k "next" intervals for this interval
        let max_j = (j + k).min(len_r);

        // Fill results for those up to k intervals
        // (count is 0..(max_j - j))
        for (count, jj) in (j..max_j).enumerate() {
            // Place the result at offset = original_i * k + count
            let offset = original_i * k + count;
            results_idx1[offset] = is1[i].idx;         // who we are
            results_idx2[offset] = is2[jj].idx;        // our neighbor
            results_dist[offset] = is2[jj].start - x + 1; // distance
        }
    }

    (results_idx1, results_idx2, results_dist)
}

pub fn nearest_previous_nonoverlapping_k(
    is1: &mut [Interval], // sorted by .start
    is2: &mut [Interval], // sorted by .end
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // Sort is1 by start, is2 by end for the "previous" intervals
    radsort::sort_by_key(is1, |i| i.start);
    radsort::sort_by_key(is2, |i| i.end);

    let len_l = is1.len();
    let len_r = is2.len();

    let mut results_idx1 = vec![-1; k * len_l];
    let mut results_idx2 = vec![-1; k * len_l];
    let mut results_dist = vec![i64::MAX; k * len_l];

    let mut j = 0; // pointer into is2

    for i in 0..len_l {
        let x = is1[i].start;

        // original_i = the interval's original index in the unsorted array
        let original_i = is1[i].idx as usize;

        // Move j until is2[j].end >= x
        // so all intervals in [0..j) have .end < x
        while j < len_r && is2[j].end <= x {
            j += 1;
        }

        // j is now the first index in is2 where .end >= x
        // so the "valid" previous intervals are is2[0..j]
        let num_valid = j;
        let start_idx = if num_valid > k { num_valid - k } else { 0 };

        // We want the last k among those valid intervals
        for (count, jj) in (start_idx..num_valid).enumerate() {
            // Place the result at offset = original_i * k + count
            let offset = original_i * k + count;
            results_idx1[offset] = is1[i].idx;  // who we are
            results_idx2[offset] = is2[jj].idx; // our neighbor
            // distance is x - is2[jj].end (+1 if you like)
            results_dist[offset] = (x - is2[jj].end + 1).abs();
        }
    }

    (results_idx1, results_idx2, results_dist)
}


pub fn nearest_nonoverlapping_k(
    is1: &mut [Interval],
    is2: &mut [Interval],
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let (nn_idx1, nn_idx2, nn_dist) = nearest_next_nonoverlapping_k(is1, is2, k);
    let (mut np_idx1, mut np_idx2, mut np_dist) = nearest_previous_nonoverlapping_k(is1, is2, k);

    // println!("next {:?}", nn_idx1);
    // println!("next {:?}", nn_idx2);
    // println!("next {:?}", nn_dist);

    // println!("prev {:?}", np_idx1);
    // println!("prev {:?}", np_idx2);
    // println!("prev {:?}", np_dist);

    let outlen = k * is1.len();

    let mut results_idx1 = Vec::with_capacity(outlen);
    let mut results_idx2 = Vec::with_capacity(outlen);
    let mut results_dist = Vec::with_capacity(outlen);
    
    for i in 0..outlen {
        if (nn_dist[i] == i64::MAX) & (np_dist[i] == i64::MAX) {
        } else if nn_dist[i] < np_dist[i] {
            results_idx1.push(nn_idx1[i]);
            results_idx2.push(nn_idx2[i]);
            results_dist.push(nn_dist[i]);
        } else {
            results_idx1.push(np_idx1[i]);
            results_idx2.push(np_idx2[i]);
            results_dist.push(np_dist[i]);
        }
    }
    (results_idx1, results_idx2, results_dist)
}