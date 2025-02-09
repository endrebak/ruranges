use rustc_hash::FxHashSet;

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

pub fn sweep_line_overlaps_nearest(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> Vec<OverlapPair> {
    // We'll collect all cross overlaps here
    let (idx1s, idxs2) = sweep_line_overlaps(chrs, starts, ends, chrs2, starts2, ends2, slack);
    println!("after sweep line overlaps {}", idx1s.len());

    let mut overlaps = Vec::new();

    for i in 0..idx1s.len() {
        overlaps.push(
            OverlapPair {
                idx: idx1s[i],
                idx2: idxs2[i],
            }
        );

    }

    overlaps
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
