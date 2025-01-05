use std::time::Instant;

use rustc_hash::FxHashSet;

use crate::sorts;

/// Returns all overlapping pairs (idx1, idx2) between intervals in set1 and set2.
/// This uses a line-sweep / active-set approach.
///
/// Algorithm steps:
///   1. Build a list of events (start & end) for each interval in both sets.
///   2. Sort events by coordinate. Where coordinates tie, put start before end.
///   3. Maintain active sets (for set1 and set2). For a start event in set1,
///      record overlap with all active in set2, then insert into active1. Etc.
///   4. Return the list of all cross-set overlaps.
pub fn sweep_line_overlaps(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i32],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
) -> (Vec<i64>, Vec<i64>) {
    let start = Instant::now();

    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() | chrs2.is_empty() {
        return (overlaps, overlaps2);
    };

    let events = sorts::build_sorted_events(chrs, starts, ends, idxs, chrs2, starts2, ends2, idxs2);
    let duration = start.elapsed();
    println!("Time elapsed building events: {:?}", duration);

    // Active sets
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    let mut current_chr: i32 = events.first().unwrap().chr;

    // Process events in ascending order of position
    for e in events {
        if e.chr != current_chr {
            active1.clear();
            active2.clear();
            current_chr = e.chr;
        }

        if e.is_start {
            // Interval is starting
            if e.first_set {
                // Overlaps with all currently active intervals in set2
                for &idx2 in active2.iter() {
                    overlaps.push(e.idx);
                    overlaps2.push(idx2);
                }
                // Now add it to active1
                active1.insert(e.idx);
            } else {
                // Overlaps with all currently active intervals in set1
                for &idx1 in active1.iter() {
                    overlaps.push(idx1);
                    overlaps2.push(e.idx);
                }
                // Now add it to active2
                active2.insert(e.idx);
            }
        } else {
            // Interval is ending
            if e.first_set {
                active1.remove(&e.idx);
            } else {
                active2.remove(&e.idx);
            }
        }
    }

    let duration = start.elapsed();
    println!("Time elapsed finding overlaps: {:?}", duration);

    (overlaps, overlaps2)
}
