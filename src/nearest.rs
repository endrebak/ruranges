
use std::time::Instant;

use radsort::sort_by_key;

use crate::{overlaps::{self, sweep_line_overlaps}, ruranges_structs::Interval, sorts::{align_interval_collections_on_chromosome, build_intervals}};

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
        let original_i = is1[i]._idx as usize;

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
        let original_i = is1[i]._idx as usize;

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
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i32],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {

    let start = Instant::now();
    let mut is1 = build_intervals(chrs, starts, ends, idxs);
    let mut is2 = build_intervals(chrs2, starts2, ends2, idxs2);
    let duration = start.elapsed();
    println!("Time elapsed building intervals: {:?}", duration);

    let outlen = k * is1.len();
    let mut results_idx1 = Vec::with_capacity(outlen);
    let mut results_idx2 = Vec::with_capacity(outlen);
    let mut results_dist = Vec::with_capacity(outlen);
    let start = Instant::now();
    let d = align_interval_collections_on_chromosome(&mut is1, &mut is2);
    let duration = start.elapsed();
    println!("Time elapsed aligning chroms: {:?}", duration);

    for (_chr, (mut is1, mut is2)) in d {
        let (nn_idx1, nn_idx2, nn_dist) = nearest_next_nonoverlapping_k(&mut is1, &mut is2, k);
        let (np_idx1, np_idx2, np_dist) = nearest_previous_nonoverlapping_k(&mut is1, &mut is2, k);

        for i in 0..(np_idx1.len()) {
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
    }
    let duration = start.elapsed();
    println!("Time elapsed finding nearest nonoverlapping intervals: {:?}", duration);
    (results_idx1, results_idx2, results_dist)
}


pub fn nearest_k(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i32],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let start = Instant::now();
    let (oidx1, oidx2) = sweep_line_overlaps(
        chrs,
        starts,
        ends,
        idxs,
        chrs2,
        starts2,
        ends2,
        idxs2,
    );
    let duration = start.elapsed();
    println!("Time elapsed: {:?}", duration);

    #[derive(Debug, Clone)]
    struct Triplet {
        pub idx: i64,
        pub idx2: i64,
        pub n: i64,
    }
    let start = Instant::now();

    let mut overlaps: Vec<Triplet> = Vec::with_capacity(oidx1.len());
    for i in 0..oidx1.len() {
        overlaps.push(Triplet {
            idx: oidx1[i], idx2: oidx2[i], n: 0,
        })
    }
    sort_by_key(&mut overlaps, |i| i.idx);
    let duration = start.elapsed();
    println!("Time elapsed: {:?}", duration);

    let (results_idx1, results_idx2, dist) = nearest_nonoverlapping_k(
        chrs,
        starts,
        ends,
        idxs,
        chrs2,
        starts2,
        ends2,
        idxs2,
        k,
    );

    let start = Instant::now();
    let mut nearest: Vec<Triplet> = Vec::with_capacity(results_idx1.len());
    for i in 0..results_idx1.len() {
        nearest.push(Triplet {
            idx: results_idx1[i], idx2: results_idx2[i], n: dist[i],
        })
    }
    let duration = start.elapsed();
    println!("Time elapsed: {:?}", duration);


    sort_by_key(&mut nearest, |i| i.idx);
    println!("nearests {:?}", &nearest[..nearest.len().min(10)]);
    println!("overlaps {:?}", &overlaps[..overlaps.len().min(10)]);

    let mut idx1: Vec<i64> = Vec::with_capacity(k * starts.len());
    let mut idx2: Vec<i64> = Vec::with_capacity(k * starts.len());
    let mut dists: Vec<i64> = Vec::with_capacity(k * starts.len());

    let start = Instant::now();

    let mut o_i = 0;
    let mut n_i = 0;

    for query_i in 0..starts.len() {
        let q_idx = query_i as i64;
        let mut count = 0;  // count of elements added for this query

        // Find ranges for both overlaps and nearest that match current query
        let mut o_end = o_i;
        while o_end < overlaps.len() && overlaps[o_end].idx == q_idx {
            o_end += 1;
        }

        let mut n_end = n_i;
        while n_end < nearest.len() && nearest[n_end].idx == q_idx {
            n_end += 1;
        }

        // Two-pointer merge until we have k elements or exhaust both lists
        while count < k && (o_i < o_end || n_i < n_end) {
            let take_overlap = if o_i < o_end {
                if n_i < n_end {
                    // Both available - take overlap since it's closer
                    true
                } else {
                    // Only overlap available
                    true
                }
            } else {
                // Only nearest available
                false
            };

            if take_overlap {
                idx1.push(overlaps[o_i].idx);
                idx2.push(overlaps[o_i].idx2);
                dists.push(overlaps[o_i].n);
                o_i += 1;
            } else {
                idx1.push(nearest[n_i].idx);
                idx2.push(nearest[n_i].idx2);
                dists.push(nearest[n_i].n);
                n_i += 1;
            }
            count += 1;
        }

        // Advance pointers to start of next query
        o_i = o_end;
        n_i = n_end;
    }
    let duration = start.elapsed();
    println!("Time elapsed: {:?}", duration);

    (idx1, idx2, dists)
}
