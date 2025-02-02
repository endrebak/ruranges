use rustc_hash::FxHashSet;

use crate::ruranges_structs::MinEvent;
use crate::sorts::build_sorted_events_single_collection_separate_outputs;
use crate::sorts;


#[derive(Debug)]
pub struct Position {
    pub chr: i64,  // you’ll need the chromosome in here to sort by (chr, pos)
    pub pos: i64,
    pub idx: i64,
}

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
    // Internally, each Position has:
    //   { chr: i64, pos: i64, idx: i64 }
    //
    // For example, `sorts::build_sorted_positions` might produce these four lists
    // ordered by (chr, pos) ascending, with the appropriate slack adjustments
    // if that's part of the logic. Make sure each Position has the correct chr.
    let (sorted_starts, sorted_ends) = build_sorted_events_single_collection_separate_outputs(chrs, starts, ends, idxs, slack);
    let (sorted_starts2, sorted_ends2) = build_sorted_events_single_collection_separate_outputs(chrs2, starts2, ends2, idxs2, 0);
    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    // If either set is empty, there can be no cross overlaps.
    if sorted_starts.is_empty() && sorted_starts2.is_empty() {
        return (overlaps, overlaps2);
    }

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
        get_event(&sorted_starts2, i2, true, false),
        get_event(&sorted_ends, i3, false, true),
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
        let mut which_list = 0; // 1 => i1, 2 => i2, 3 => i3, 4 => i4

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

use rustc_hash::{FxHashMap};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

/// A simple struct to store pairs of intervals and their distance.
/// We implement `Ord`/`PartialOrd` so that a BinaryHeap can treat
/// it like a *max*-heap by distance. This way, popping removes the
/// *largest* distance first, ensuring only the smallest k are kept.
#[derive(Debug, Eq, PartialEq)]
pub struct NearestInterval {
    pub idx1: i64,
    pub idx2: i64,
    pub distance: i64,
}

impl Ord for NearestInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        // By default, Rust's BinaryHeap is a max-heap. We invert the usual comparison
        // so that *larger* distances are popped first, leaving the smallest distances
        // in the heap if we exceed size k.
        other.distance.cmp(&self.distance)
    }
}

impl PartialOrd for NearestInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


/// A helper function to push `(idx1, idx2, distance)` into the correct min-heap
/// for each interval.  We keep only the smallest k distances by popping from
/// the top if size > k.
fn push_heap(
    heaps: &mut FxHashMap<i64, BinaryHeap<NearestInterval>>,
    interval_idx: i64,
    new_item: NearestInterval,
    k: usize,
) {

    println!("pushing: {:?}", new_item);
    let h = heaps.entry(interval_idx).or_insert_with(BinaryHeap::new);
    h.push(new_item);
    if h.len() > k {
        h.pop(); // pop the largest distance, keeping only k smallest
    }
}

/// Returns three vectors of equal length: (idx1, idx2, distance).
/// For each interval from set1, we keep a size-k min-heap of nearest intervals from set2,
/// and vice versa. We store "overlapping intervals," plus the 'last k ended' from set2
/// to the left, plus the 'next k starts' from set2 to the right, etc.
/// Once an interval ends, we pop its min-heap into the output vectors.
pub fn sweep_line_overlaps_merged_with_heaps(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    slack: i64,
    k: usize,        // how many intervals we keep in each min-heap
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // These will be our final outputs.
    let mut overlaps_1 = Vec::new();
    let mut overlaps_2 = Vec::new();
    let mut distances = Vec::new();

    if chrs.is_empty() || chrs2.is_empty() {
        return (overlaps_1, overlaps_2, distances);
    }

    // Build your sorted events. For simplicity, we assume your existing
    // `sorts::build_sorted_events(...)` returns an ordered list of Events.
    let events = sorts::build_sorted_events(
        chrs, starts, ends, idxs,
        chrs2, starts2, ends2, idxs2,
        slack,
    );

    // *** CHANGED: Instead of a single "active1" or "active2", we still keep
    //              sets of active intervals, but also keep a separate
    //              heap mapping for each *currently active* interval.
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    // *** CHANGED: For each *interval* (whether in set1 or set2) that is active,
    //              we store a size-k min-heap of nearest intervals from the *other* set.
    //              We can keep two separate maps for clarity.
    let mut heaps1: FxHashMap<i64, BinaryHeap<NearestInterval>> = FxHashMap::default();
    let mut heaps2: FxHashMap<i64, BinaryHeap<NearestInterval>> = FxHashMap::default();

    // Keep track of last k ended intervals from set2
    let mut last_ended_set2: VecDeque<i64> = VecDeque::with_capacity(k);

    let n = events.len();
    let mut current_chr = events[0].chr;

    // We'll iterate by index so we can do a forward look-up for "next k" intervals.
    for i in 0..n {
        let e = &events[i];

        // If we switched chromosomes, reset everything
        if e.chr != current_chr {
            // clear active sets
            active1.clear();
            active2.clear();

            // clear any leftover heaps
            heaps1.clear();
            heaps2.clear();

            // if you do NOT want last-ended-set2 intervals to carry over between chromosomes
            last_ended_set2.clear();

            current_chr = e.chr;
        }

        // For convenience, define a small closure that helps us compute distance
        // between a set1 event and a set2 event. Here, we simply do absolute difference
        // in positions. You can adapt to midpoint differences, etc.
        let distance = |pos1: i64, pos2: i64| (pos1 - pos2).abs();

        if e.is_start {
            // A start event
            if e.first_set {
                // e belongs to set1

                // 1) Overlaps: for all active2 intervals, push into heaps
                //    We push into the set1's heap, and also the set2's heap for each active2.
                for &idx2 in active2.iter() {
                    // you might want to retrieve the 'pos' of that idx2 from your original arrays
                    // or store it in a small structure. Here we guess the 'starts2[idx2]' as position:
                    let dist_val = distance(e.pos, starts2[idx2 as usize]);
                    let item = NearestInterval {
                        idx1: e.idx,
                        idx2,
                        distance: dist_val,
                    };

                    // push into set1's heap
                    push_heap(&mut heaps1, e.idx, item, k);

                    // push a symmetrical item into set2's heap
                    // (idx2's perspective: the nearest set1 intervals)
                    let item2 = NearestInterval {
                        idx1: e.idx,
                        idx2,
                        distance: dist_val,
                    };
                    push_heap(&mut heaps2, idx2, item2, k);
                }

                // 2) "Last k ended set2" logic: push each ended_idx2 into set1's heap
                //    (and symmetrically, push the pair to set2's heap so that from
                //     that set2 interval's perspective, we track this set1 as well).
                for &ended_idx2 in &last_ended_set2 {
                    let dist_val = distance(e.pos, starts2[ended_idx2 as usize]);
                    let item = NearestInterval {
                        idx1: e.idx,
                        idx2: ended_idx2,
                        distance: dist_val,
                    };
                    push_heap(&mut heaps1, e.idx, item, k);

                    let item2 = NearestInterval {
                        idx1: e.idx,
                        idx2: ended_idx2,
                        distance: dist_val,
                    };
                    push_heap(&mut heaps2, ended_idx2, item2, k);
                }

                // Mark this set1 interval as active, initialize its heap
                active1.insert(e.idx);
                heaps1.entry(e.idx).or_insert_with(BinaryHeap::new);

            } else {
                // e belongs to set2

                // Overlaps: for all active1 intervals, push into heaps
                for &idx1 in active1.iter() {
                    let dist_val = distance(starts[idx1 as usize], e.pos);
                    let item = NearestInterval {
                        idx1,
                        idx2: e.idx,
                        distance: dist_val,
                    };
                    push_heap(&mut heaps2, e.idx, item, k);

                    // symmetrical entry for the set1 perspective
                    let item2 = NearestInterval {
                        idx1,
                        idx2: e.idx,
                        distance: dist_val,
                    };
                    push_heap(&mut heaps1, idx1, item2, k);
                }

                // Mark set2 interval as active
                active2.insert(e.idx);
                heaps2.entry(e.idx).or_insert_with(BinaryHeap::new);
            }

        } else {
            // An interval is *ending*
            if e.first_set {
                // A set1 interval is ending
                active1.remove(&e.idx);

                // 3) "Next k intervals in set2" logic:
                //    We'll do a simple forward search from i+1 to find up to k intervals
                //    that start in set2 (on the same chromosome, with pos >= e.pos).
                let mut found_count = 0;
                let mut j = i + 1;
                while j < n && found_count < k {
                    let e2 = &events[j];
                    if e2.chr != e.chr {
                        break;
                    }
                    if e2.is_start && !e2.first_set && e2.pos >= e.pos {
                        // This is a set2 start after (or at) e.pos
                        let dist_val = distance(e.pos, e2.pos);
                        let item = NearestInterval {
                            idx1: e.idx,
                            idx2: e2.idx,
                            distance: dist_val,
                        };
                        push_heap(&mut heaps1, e.idx, item, k);

                        // symmetrical push from that set2 interval's perspective
                        let item2 = NearestInterval {
                            idx1: e.idx,
                            idx2: e2.idx,
                            distance: dist_val,
                        };
                        push_heap(&mut heaps2, e2.idx, item2, k);

                        found_count += 1;
                    }
                    j += 1;
                }

                // *** CHANGED: Now that the set1 interval e.idx is ending (closed),
                // we finalize its results by draining its heap into the output vectors.
                if let Some(mut heap) = heaps1.remove(&e.idx) {
                    while let Some(NearestInterval { idx1, idx2, distance }) = heap.pop() {
                        overlaps_1.push(idx1);
                        overlaps_2.push(idx2);
                        distances.push(distance);
                    }
                }

            } else {
                // A set2 interval is ending
                active2.remove(&e.idx);

                // Add it to our queue of recently ended intervals in set2
                last_ended_set2.push_back(e.idx);
                if last_ended_set2.len() > k {
                    last_ended_set2.pop_front();
                }

                // *** CHANGED: finalize the set2 interval's heap as well.
                if let Some(mut heap) = heaps2.remove(&e.idx) {
                    while let Some(NearestInterval { idx1, idx2, distance }) = heap.pop() {
                        overlaps_1.push(idx1);
                        overlaps_2.push(idx2);
                        distances.push(distance);
                    }
                }
            }
        }
    }

    // Return the three parallel vectors
    (overlaps_1, overlaps_2, distances)
}
