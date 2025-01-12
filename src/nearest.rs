
use radsort::sort_by_key;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;
use std::time::Instant;

// Assume these are your own structures/modules
use crate::{
    ruranges_structs::{Event, Nearest},
    sorts::build_sorted_events_single_position,
};

/// This is your existing function signature:
pub fn sweep_line_k_nearest(
    chrs: &[i64],
    pos: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    pos2: &[i64],
    idxs2: &[i64],
    invert_pos: bool,
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // 1) Build event arrays:
    //    - from Set1: we only care about "end" positions
    //    - from Set2: we only care about "start" positions
    let ends_events = build_sorted_events_single_position(chrs, pos, idxs, false, true, invert_pos, 0);
    let starts_events = build_sorted_events_single_position(chrs2, pos2, idxs2, true, false, invert_pos, 0);

    // We'll store the results in these vectors:
    let mut out_idxs1 = Vec::new();
    let mut out_idxs2 = Vec::new();
    let mut out_dists = Vec::new();
    let mut out_nearest: Vec<Nearest> = Vec::new();

    // 2) We'll do a two-pointer approach:
    let mut j = 0; // pointer into starts_events

    // For each "end" event in set1 (the intervals we want to find neighbors for):
    for e1 in ends_events.iter() {
        // We'll skip forward in starts_events until we find the first position
        // that is on the same chr and pos >= e1.pos. We also skip if we are behind in chr.
        while j < starts_events.len() {
            let e2 = &starts_events[j];

            // If we've passed the chromosome or if e2 is already at or beyond e1:
            //   break the while loop if e2 is on the same chr and e2.pos >= e1.pos
            //   or if e2.chr > e1.chr (we're definitely at a future chr now).
            if e2.chr > e1.chr {
                // We've moved beyond the chromosome entirely, so we don't advance j further.
                break;
            }
            if e2.chr == e1.chr && e2.pos >= e1.pos {
                // found the first position in set2 that's >= e1.pos on the same chromosome
                break;
            }

            j += 1;
        }

        // Now, gather up to k distinct positions from j forward,
        // but collect all events that share that same position.
        let mut distinct_count = 0;
        let mut local_ptr = j;
        let mut last_pos = None;

        while local_ptr < starts_events.len() {
            let e2 = &starts_events[local_ptr];
            // If we've changed chromosomes, we're done for this e1.
            if e2.chr != e1.chr {
                break;
            }

            // Check if we have exceeded the distinct positions limit
            if distinct_count >= k {
                break;
            }

            // If this is the *first* or a *new* position:
            if last_pos.map_or(true, |p| p != e2.pos) {
                distinct_count += 1;
                last_pos = Some(e2.pos);
            }

            // We include *all* events with this same e2.pos
            // because they count collectively as just one distinct position.
            let current_pos = e2.pos;
            while local_ptr < starts_events.len() {
                let e2_same = &starts_events[local_ptr];
                // stop if chromosome changes or position changes
                if e2_same.chr != e1.chr || e2_same.pos != current_pos {
                    break;
                }
                // Distance formula. Adjust to match your "bookended=1" scheme:
                // For example:
                //   dist = (start2 - end1) + 1
                // so if start2 == end1+1 => dist=1 (adjacent)
                // if start2 == end1+2 => dist=2, etc.
                // That means:
                let dist = (e2_same.pos - e1.pos) + 1;
                
                out_nearest.push(
                    Nearest {
                        idx: e1.idx,
                        idx2: e2_same.idx,
                        distance: dist,
                    }
                );

                local_ptr += 1;
            }
        }
    }

    sort_by_key(&mut out_nearest, |e| e.distance);
    sort_by_key(&mut out_nearest, |e| e.idx);

    for e in out_nearest {
        out_idxs1.push(e.idx);
        out_idxs2.push(e.idx2);
        out_dists.push(e.distance);
    }

    // Return the 3 vectors
    (out_idxs1, out_idxs2, out_dists)
}

/// Combine left & right search results, then pick up to `k` distinct distances for each idx1.
pub fn pick_k_distances_combined(
    l_idxs1: &[i64],
    l_idxs2: &[i64],
    l_dists: &[i64],
    r_idxs1: &[i64],
    r_idxs2: &[i64],
    r_dists: &[i64],
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // 1) Combine them into a single Vec of tuples.
    //    We don't strictly need to store whether they're left or right, 
    //    unless you want that info later.
    
    let combined = merge_sorted_by_idx_dist(l_idxs1, l_idxs2, l_dists, r_idxs1, r_idxs2, r_dists);

    // 3) We'll now pick up to `k` distinct distances per `idx1`.
    let mut out_idxs1 = Vec::new();
    let mut out_idxs2 = Vec::new();
    let mut out_dists = Vec::new();

    let mut i = 0;
    while i < combined.len() {
        let current_idx1 = combined[i].0;

        // Gather all rows for this `idx1` in the slice [start..i_end)
        let start = i;
        while i < combined.len() && combined[i].0 == current_idx1 {
            i += 1;
        }
        // Now the sub-range [start..i) includes all tuples for that `idx1`

        let mut distinct_count = 0;
        let mut last_dist = None;

        // 4) Within this sub-range, iterate in ascending dist order
        for j in start..i {
            let dist = combined[j].2;
            // If we've hit a *new* distance, increment our distinct distance count
            if Some(dist) != last_dist {
                distinct_count += 1;
                last_dist = Some(dist);
            }
            // If we've exceeded `k` distinct distances, stop collecting for this idx1
            if distinct_count > k {
                break;
            }
            // Otherwise, add it to our output
            out_idxs1.push(combined[j].0);
            out_idxs2.push(combined[j].1);
            out_dists.push(dist);
        }
    }

    (out_idxs1, out_idxs2, out_dists)
}


/// Merge two sorted lists (by (idx1, dist)) into a single sorted Vec.
/// left_* and right_* each represent (idx1, idx2, dist) from two searches,
/// sorted by (idx1, dist).
fn merge_sorted_by_idx_dist(
    l_idxs1: &[i64],
    l_idxs2: &[i64],
    l_dists: &[i64],
    r_idxs1: &[i64],
    r_idxs2: &[i64],
    r_dists: &[i64],
) -> Vec<(i64, i64, i64)> {
    let mut combined = Vec::with_capacity(l_idxs1.len() + r_idxs1.len());
    let mut i = 0;
    let mut j = 0;

    while i < l_idxs1.len() && j < r_idxs1.len() {
        // Left tuple
        let left = (l_idxs1[i], l_idxs2[i], l_dists[i]);
        // Right tuple
        let right = (r_idxs1[j], r_idxs2[j], r_dists[j]);

        // Compare (idx1, dist)
        if left.0 < right.0 {
            combined.push(left);
            i += 1;
        } else if left.0 > right.0 {
            combined.push(right);
            j += 1;
        } else {
            // idx1 are equal => compare dist
            if left.2 <= right.2 {
                combined.push(left);
                i += 1;
            } else {
                combined.push(right);
                j += 1;
            }
        }
    }

    // If any remaining items on the left side, push them
    while i < l_idxs1.len() {
        combined.push((l_idxs1[i], l_idxs2[i], l_dists[i]));
        i += 1;
    }
    // If any remaining items on the right side, push them
    while j < r_idxs1.len() {
        combined.push((r_idxs1[j], r_idxs2[j], r_dists[j]));
        j += 1;
    }

    combined
}
