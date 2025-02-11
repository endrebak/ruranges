use rustc_hash::FxHashSet;
use indexmap::IndexSet;

use crate::ruranges_structs::{MinEvent, OverlapPair};
use crate::sorts;

/// Perform a four-way merge sweep to find cross overlaps.

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum WhichList {
    StartSet1,
    EndSet1,
    StartSet2,
    EndSet2,
}

impl WhichList {
    #[inline]
    fn is_start(&self) -> bool {
        match self {
            WhichList::StartSet1 | WhichList::StartSet2 => true,
            WhichList::EndSet1 | WhichList::EndSet2 => false,
        }
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
pub fn sweep_line_overlaps(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> (Vec<usize>, Vec<usize>) {
    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() | chrs2.is_empty() {
        return (overlaps, overlaps2);
    };

    let events = sorts::build_sorted_events(
        chrs, starts, ends, chrs2, starts2, ends2, slack,
    );

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
                };
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

pub fn sweep_line_overlaps_set1(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> Vec<usize> {
    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() | chrs2.is_empty() {
        return overlaps;
    };

    let events = sorts::build_sorted_events(
        chrs, starts, ends, chrs2, starts2, ends2, slack,
    );

    // Active sets
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    let mut current_chr: i64 = events.first().unwrap().chr;

    // Process events in ascending order of position
    for e in events {
        if e.chr != current_chr {
            active1.clear();
            current_chr = e.chr;
        }

        if e.is_start {
            // Interval is starting
            if e.first_set {
                // Overlaps with all currently active intervals in set2
                for &idx2 in active2.iter() {
                    overlaps.push(e.idx);
                    overlaps2.push(idx2);
                };
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

    overlaps
}

pub fn sweep_line_overlaps_overlap_pair(
    sorted_starts: &[MinEvent],  // set 1 starts
    sorted_ends: &[MinEvent],    // set 1 ends
    sorted_starts2: &[MinEvent], // set 2 starts
    sorted_ends2: &[MinEvent],   // set 2 ends
) -> Vec<OverlapPair> {
    let mut out_idxs = Vec::new();
    // Quick check: if no starts exist in either set, no overlaps.
    if sorted_starts.is_empty() || sorted_starts2.is_empty() {
        return out_idxs;
    }
    // Active intervals for set1, set2
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();
    // Pointers into each list
    let mut i1 = 0usize; // pointer into sorted_starts  (set 1)
    let mut i2 = 0usize; // pointer into sorted_starts2 (set 2)
    let mut i3 = 0usize; // pointer into sorted_ends    (set 1)
    let mut i4 = 0usize; // pointer into sorted_ends2   (set 2)
    // Figure out the very first chromosome we encounter (if any):
    // We'll look at the heads of each list and pick the lexicographically smallest.
    let first_candidate = pick_winner_of_four(
        sorted_starts.get(i1).map(|e| (WhichList::StartSet1, e)),
        sorted_starts2.get(i2).map(|e| (WhichList::StartSet2, e)),
        sorted_ends.get(i3).map(|e| (WhichList::EndSet1, e)),
        sorted_ends2.get(i4).map(|e| (WhichList::EndSet2, e)),
    );
    // Unwrap the first candidate’s chromosome
    let mut current_chr = first_candidate.unwrap().1.chr;
    // Main sweep-line loop
    while i1 < sorted_starts.len()
        || i2 < sorted_starts2.len()
        || i3 < sorted_ends.len()
        || i4 < sorted_ends2.len()
    {
        let (which_list, event) = if let Some((which_list, event)) = pick_winner_of_four(
            sorted_starts.get(i1).map(|e| (WhichList::StartSet1, e)),
            sorted_starts2.get(i2).map(|e| (WhichList::StartSet2, e)),
            sorted_ends.get(i3).map(|e| (WhichList::EndSet1, e)),
            sorted_ends2.get(i4).map(|e| (WhichList::EndSet2, e)),
        ) {
            (which_list, event)
        } else {
            break;
        };
        // If we've moved to a new chromosome, reset active sets
        if event.chr != current_chr {
            active1.clear();
            active2.clear();
            current_chr = event.chr;
        }
        // Advance the pointer for whichever list we took an event from
        match which_list {
            WhichList::StartSet1 => {
                for &idx2 in active2.iter() {
                    out_idxs.push(OverlapPair {
                        idx: event.idx,
                        idx2: idx2,
                    })
                }
                // Now add it to active1
                active1.insert(event.idx);
                i1 += 1
            }
            WhichList::StartSet2 => {
                for &idx1 in active1.iter() {
                    out_idxs.push(OverlapPair {
                        idx: idx1,
                        idx2: event.idx,
                    })
                }
                // Now add it to active2
                active2.insert(event.idx);
                i2 += 1
            }
            WhichList::EndSet1 => {
                active1.remove(&event.idx);
                i3 += 1
            }
            WhichList::EndSet2 => {
                active2.remove(&event.idx);
                i4 += 1
            }
        }
    }
    out_idxs
}


pub fn sweep_line_first_overlaps(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> (Vec<usize>, Vec<usize>) {
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() || chrs2.is_empty() {
        return (overlaps, overlaps2);
    }

    // Build a sorted list of events. 
    // Each event has { chr, pos, idx, first_set, is_start }.
    // (Not shown here, but presumably something like:)
    let events = sorts::build_sorted_events(chrs, starts, ends, chrs2, starts2, ends2, slack);

    // Use IndexSet instead of Vec.
    let mut active1 = IndexSet::new();
    let mut active2 = IndexSet::new();

    let mut current_chr = events.first().unwrap().chr;

    for e in events {
        // Whenever chromosome changes, clear the active sets
        if e.chr != current_chr {
            active1.clear();
            active2.clear();
            current_chr = e.chr;
        }

        if e.is_start {
            // Interval is starting
            if e.first_set {
                // This is an interval from set1; try overlapping with the first active from set2
                if let Some(&idx2) = active2.get_index(0) {
                    // Record only this single overlap
                    overlaps.push(e.idx);
                    overlaps2.push(idx2);
                    // Do NOT add `e.idx` to active1, since we only want the first overlap
                } else {
                    // No active intervals in set2, so keep this one active for future
                    active1.insert(e.idx);
                }
            } else {
                // This is an interval from set2; try overlapping with the first active from set1
                if let Some(&idx1) = active1.get_index(0) {
                    overlaps.push(idx1);
                    overlaps2.push(e.idx);
                    // Do NOT add `e.idx` to active2
                } else {
                    active2.insert(e.idx);
                }
            }
        } else {
            // Interval is ending
            if e.first_set {
                // Remove from active1
                active1.remove(&e.idx);
            } else {
                // Remove from active2
                active2.remove(&e.idx);
            }
        }
    }

    (overlaps, overlaps2)
}


pub fn sweep_line_overlaps_containment(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> (Vec<usize>, Vec<usize>) {
    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() | chrs2.is_empty() {
        return (overlaps, overlaps2);
    };

    let events = sorts::build_sorted_events(
        chrs, starts, ends, chrs2, starts2, ends2, slack,
    );

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
                };
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

fn pick_winner_of_four<'a>(
    s1: Option<(WhichList, &'a MinEvent)>,
    s2: Option<(WhichList, &'a MinEvent)>,
    e1: Option<(WhichList, &'a MinEvent)>,
    e2: Option<(WhichList, &'a MinEvent)>,
) -> Option<(WhichList, &'a MinEvent)> {
    let starts_winner = pick_winner_of_two_choose_first_if_equal(s1, e1);
    let ends_winner = pick_winner_of_two_choose_first_if_equal(s2, e2);
    pick_winner_of_two_choose_first_if_equal(starts_winner, ends_winner)
}

fn pick_winner_of_two_choose_first_if_equal<'a>(
    a: Option<(WhichList, &'a MinEvent)>,
    b: Option<(WhichList, &'a MinEvent)>,
) -> Option<(WhichList, &'a MinEvent)> {
    match (a, b) {
        (None, None) => None,
        (Some(x), None) => Some(x),
        (None, Some(y)) => Some(y),
        (Some((wh_a, ev_a)), Some((wh_b, ev_b))) => {
            // Compare by chromosome
            if ev_a.chr < ev_b.chr {
                return Some((wh_a, ev_a));
            } else if ev_b.chr < ev_a.chr {
                return Some((wh_b, ev_b));
            }
            // Same chr => compare by pos
            if ev_a.pos < ev_b.pos {
                return Some((wh_a, ev_a));
            } else if ev_b.pos < ev_a.pos {
                return Some((wh_b, ev_b));
            }
            // Same (chr, pos) => tie break: end < start
            let a_is_end = !wh_a.is_start();
            let b_is_end = !wh_b.is_start();
            match (a_is_end, b_is_end) {
                // If both are ends or both are starts, just pick either. We'll pick `a`.
                (true, true) | (false, false) => Some((wh_a, ev_a)),
                // If only one is end, that one is “smaller”
                (true, false) => Some((wh_a, ev_a)),
                (false, true) => Some((wh_b, ev_b)),
            }
        }
    }
}
