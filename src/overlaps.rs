use rustc_hash::FxHashSet;

use crate::ruranges_structs::MinEvent;
use crate::sorts::build_sorted_events_single_collection_separate_outputs;
use crate::sorts;


/// Perform a four-way merge sweep to find cross overlaps.
pub fn sweep_line_overlaps(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    slack: i64,
) -> (Vec<i64>, Vec<i64>) {
    // Build the four sorted lists of positions:
    //
    //   - sorted_starts:   starts from set 1
    //   - sorted_ends:     ends   from set 1
    //   - sorted_starts2:  starts from set 2
    //   - sorted_ends2:    ends   from set 2
    //
    // Internally, each MinEvent has:
    //   { chr: i64, pos: i64, idx: i64 }
    //
    // For example, `sorts::build_sorted_positions` might produce these four lists
    // ordered by (chr, pos) ascending, with the appropriate slack adjustments
    // if that's part of the logic. Make sure each MinEvent has the correct chr.
    let (sorted_starts, sorted_ends) = build_sorted_events_single_collection_separate_outputs(chrs, starts, ends, idxs, slack);
    let (sorted_starts2, sorted_ends2) = build_sorted_events_single_collection_separate_outputs(chrs2, starts2, ends2, idxs2, 0);
    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    // If either set is empty, there can be no cross overlaps.
    if sorted_starts.is_empty() && sorted_starts2.is_empty() {
        return (overlaps, overlaps2);
    }

    let k = 2;

    // Active sets
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    // Indices for each of the four sorted lists
    let mut i1 = 0usize; // pointer into sorted_starts  (set 1, is_start = true)
    let mut i2 = 0usize; // pointer into sorted_ends    (set 1, is_start = false)
    let mut i3 = 0usize; // pointer into sorted_starts2 (set 2, is_start = true)
    let mut i4 = 0usize; // pointer into sorted_ends2   (set 2, is_start = false)

    // A small helper to get the "next" event from one of the four lists if available
    // We also embed the additional info: (is_start, first_set)
    // - (is_start = true,  first_set = true)  => from sorted_starts
    // - (is_start = false, first_set = true)  => from sorted_ends
    // - (is_start = true,  first_set = false) => from sorted_starts2
    // - (is_start = false, first_set = false) => from sorted_ends2
    //
    // We return `None` if the index is out of range for that list.
    fn get_event(
        list: &[MinEvent],
        idx: usize,
        is_start: bool,
        first_set: bool,
    ) -> Option<(i64, i64, i64, bool, bool)> {
        if idx < list.len() {
            let p = &list[idx];
            Some((p.chr, p.pos, p.idx, is_start, first_set))
        } else {
            None
        }
    }

    // Grab the chromosome of the first event we’ll process, if any
    // We do this by picking the lexicographically smallest among all heads:
    let first = [
        get_event(&sorted_starts, i1, true, true),
        get_event(&sorted_ends, i2, false, true),
        get_event(&sorted_starts2, i3, true, false),
        get_event(&sorted_ends2, i4, false, false),
    ]
    .into_iter()
    .flatten() // drop the Nones
    .min_by_key(|(chr, pos, _, _, _)| (*chr, *pos));

    // If there’s no first event at all, we’re done
    if first.is_none() {
        return (overlaps, overlaps2);
    }
    let mut current_chr = first.unwrap().0; // unwrap and take the chr of the first event

    // While there are still elements in any of the 4 lists...
    while i1 < sorted_starts.len()
        || i2 < sorted_ends.len()
        || i3 < sorted_starts2.len()
        || i4 < sorted_ends2.len()
    {
        // Find the smallest (chr, pos) among the 4 candidates
        // [i1-th of sorted_starts, i2-th of sorted_ends, i3-th of sorted_starts2, i4-th of sorted_ends2].
        let mut candidate: Option<(i64, i64, i64, bool, bool)> = None; // (chr, pos, idx, is_start, first_set)
        let mut which_list: i8 = 0; // 1 => i1, 2 => i2, 3 => i3, 4 => i4

        // Check each list’s head if available, track the min
        for (lst_id, ev) in [
            (1, get_event(&sorted_starts, i1, true, true)),
            (2, get_event(&sorted_ends, i2, false, true)),
            (3, get_event(&sorted_starts2, i3, true, false)),
            (4, get_event(&sorted_ends2, i4, false, false)),
        ] {
            if let Some(e) = ev {
                // e has the shape (chr, pos, idx, is_start, first_set)
                if let Some(current) = candidate {
                    // compare e to current
                    if (e.0, e.1) < (current.0, current.1) {
                        candidate = Some(e);
                        which_list = lst_id;
                    }
                } else {
                    candidate = Some(e);
                    which_list = lst_id;
                }
            }
        }

        // If no candidate was found, we’re done
        let (chr, pos, idx, is_start, first_set) = match candidate {
            None => break,
            Some(x) => x,
        };

        // If we moved to a new chromosome, clear active sets
        if chr != current_chr {
            active1.clear();
            active2.clear();
            current_chr = chr;
        }

        // Apply the sweep-line logic
        if is_start {
            // Interval is starting
            if first_set {
                add_nearest_intervals_to_the_left(
                    chr,                // current chromosome
                    pos,
                    &sorted_ends2,      // all "end" events from set2
                    i4,                 // pointer to the next "end" event in set2
                    k,                  // how many left intervals you want
                );
                // Overlaps with all currently active intervals in set2
                for &idx2 in active2.iter() {
                    overlaps.push(idx);
                    overlaps2.push(idx2);
                }
                // Now add it to active1
                active1.insert(idx);
            } else {
                // Overlaps with all currently active intervals in set1
                for &idx1 in active1.iter() {
                    overlaps.push(idx1);
                    overlaps2.push(idx);
                }
                // Now add it to active2
                active2.insert(idx);
            }
        } else {
            // Interval is ending
            if first_set {
                active1.remove(&idx);
            } else {
                active2.remove(&idx);
            }
        }

        // Advance the pointer for whichever list we took an element from
        match which_list {
            1 => i1 += 1,
            2 => i2 += 1,
            3 => i3 += 1,
            4 => i4 += 1,
            _ => {}
        }
    }

    (overlaps, overlaps2)
}



#[derive(Debug, Clone, Hash)]
pub struct Nearest {
    pub distance: i64,
    pub idx2: i64,
}


fn add_nearest_intervals_to_the_left(
    chr: i64,
    current_interval_pos: i64,
    sorted_ends2: &[MinEvent],
    i4: usize,
    k: usize,
) {
    let mut distances: Vec<Nearest> = Vec::with_capacity(k);
    // We'll look backward through the intervals that ended in set2,
    // which are at positions [0..i4] in sorted_ends2.
    // We want up to k intervals on the same chromosome (chr).
    let mut count = 0;

    // i4 is the "next to process" in ends2, so the last *finished* event is i4 - 1.
    // Make sure we don't panic if i4 == 0.
    let mut idx = i4;
    while idx > 0 && count < k {
        idx -= 1;
        let ev = &sorted_ends2[idx];
        // If we only want intervals on the same chromosome, break if we see a mismatch.
        if ev.chr != chr {
            break;
        }
        // ev.idx is the index of the ended interval in set2
        idxs2.push(ev.idx);
        distances.push(ev.pos - (current_interval_pos + 1));

        count += 1;
    }
}



/// Returns all overlapping pairs (idx1, idx2) between intervals in set1 and set2.
/// This uses a line-sweep / active-set approach.
///
/// Algorithm steps:
///   1. Build a list of events (start & end) for each interval in both sets.
///   2. Sort events by coordinate. Where coordinates tie, put start before end.
///   3. Maintain active sets (for set1 and set2). For a start event in set1,
///      record overlap with all active in set2, then insert into active1. Etc.
///   4. Return the list of all cross-set overlaps.
pub fn sweep_line_overlaps2(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    slack: i64,
) -> (Vec<i64>, Vec<i64>) {

    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() | chrs2.is_empty() {
        return (overlaps, overlaps2);
    };

    let events = sorts::build_sorted_events(chrs, starts, ends, idxs, chrs2, starts2, ends2, idxs2, slack);

    // Active sets
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    let mut current_chr: i64 = events.first().unwrap().chr;

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

    (overlaps, overlaps2)
}