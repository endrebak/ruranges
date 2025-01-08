use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;
use std::time::Instant;

// Assume these are your own structures/modules
use crate::{
    ruranges_structs::{Event, Nearest},
    sorts,
};

pub fn sweep_line_k_nearest(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    k: usize,
    overlaps: bool,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let start = Instant::now();

    // -- Original return vectors --
    let mut overlaps_left = Vec::new();
    let mut overlaps_right = Vec::new();
    let mut distances = Vec::new();

    // -- Temporary storage for nearest candidates --
    let mut found_nearest: FxHashMap<i64, Vec<Nearest>> = FxHashMap::default();

    if chrs.is_empty() || chrs2.is_empty() {
        // If either set is empty, just return empty vectors.
        return (overlaps_left, overlaps_right, distances);
    }

    // Build combined / sorted event list
    let events = sorts::build_sorted_events(chrs, starts, ends, idxs, chrs2, starts2, ends2, idxs2, 0);

    let duration = start.elapsed();
    println!("Time elapsed building events: {:?}", duration);

    // Active sets for normal overlap logic
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    // For "k nearest previous"
    let mut right_buffer: VecDeque<Event> = VecDeque::with_capacity(k);

    // For "k nearest next"
    // left_idx => (end_position_of_left, how_many_right_starts_remaining)
    let mut next_counts: FxHashMap<i64, (i64, usize)> = FxHashMap::default();

    let mut current_chr: i64 = events.first().unwrap().chr;

    // -- Process events in ascending order of position --
    for e in events {
        if e.chr != current_chr {
            // If we switched chromosome, clear everything
            active1.clear();
            active2.clear();
            next_counts.clear();
            right_buffer.clear();
            current_chr = e.chr;
        }

        if e.is_start {
            // Interval is STARTING
            if e.first_set {
                // -- LEFT interval is starting --
                // 1) "k nearest previous" logic
                for right_event in right_buffer.iter() {
                    // Instead of pushing directly into overlaps_left/etc.,
                    // accumulate in found_nearest:
                    found_nearest.entry(e.idx).or_default().push(Nearest {
                        distance: e.pos - right_event.pos + 1,
                        idx: right_event.idx,
                    });
                }

                if overlaps {
                    // 2) Overlapping logic: all currently active2 intervals
                    for &idx2 in &active2 {
                        found_nearest.entry(e.idx).or_default().push(Nearest {
                            distance: 0, // distance=0 for overlaps
                            idx: idx2,
                        });
                    }
                }

                // Mark this left interval as active
                active1.insert(e.idx);
            } else {
                // -- RIGHT interval is starting --
                // 1) "k nearest next" logic: for each left interval that
                //    requested the next k starts from the right set
                if !next_counts.is_empty() {
                    let mut to_remove = Vec::new();
                    for (left_idx, (pos, count_left)) in next_counts.iter_mut() {
                        found_nearest.entry(*left_idx).or_default().push(Nearest {
                            distance: e.pos - *pos + 1,
                            idx: e.idx,
                        });

                        // Decrement how many more right starts we pair with
                        *count_left -= 1;
                        if *count_left == 0 {
                            to_remove.push(*left_idx);
                        }
                    }
                    // remove left indices that have used up their quota
                    for left_idx in to_remove {
                        next_counts.remove(&left_idx);
                    }
                }

                if overlaps {
                    // 2) Overlapping logic: all currently active1 intervals
                    for &idx1 in &active1 {
                        found_nearest.entry(idx1).or_default().push(Nearest {
                            distance: 0,
                            idx: e.idx,
                        });
                    }
                }

                // Mark this right interval as active
                active2.insert(e.idx);

                // We do NOT add this right interval to right_buffer here,
                // because we only track ENDING intervals in right_buffer
                // for the "k nearest previous" logic.
            }
        } else {
            // Interval is ENDING
            if e.first_set {
                // -- LEFT interval is ending --
                // Remove from active1
                active1.remove(&e.idx);

                // Mark that we want "k nearest next" from the right set
                next_counts.insert(e.idx, (e.pos, k));
            } else {
                // -- RIGHT interval is ending --
                // Remove from active2
                active2.remove(&e.idx);

                // "k nearest previous": we add this ended right interval
                // to right_buffer (rolling buffer of size k)
                if right_buffer.len() == k {
                    right_buffer.pop_front();
                }
                right_buffer.push_back(e);
            }
        }
    }

    // ------------------------------------------------
    // Post-processing: sort & keep only the top k nearest
    // ------------------------------------------------
    for (left_idx, mut nearest_vec) in found_nearest {
        // 1) Sort each vector by distance ascending
        nearest_vec.sort_unstable_by_key(|n| n.distance);

        // 2) Keep only the top k (smallest distance)
        //    (If nearest_vec.len() < k, .take(k) just yields all)
        for candidate in nearest_vec.into_iter().take(k) {
            // 3) Push into the original vectors
            overlaps_left.push(left_idx);
            overlaps_right.push(candidate.idx);
            distances.push(candidate.distance);
        }
    }

    let duration = start.elapsed();
    println!("Time elapsed finding pairs: {:?}", duration);

    // Return the original vectors, now populated with
    // the (up to) k-nearest intervals (and/or overlaps).
    (overlaps_left, overlaps_right, distances)
}
