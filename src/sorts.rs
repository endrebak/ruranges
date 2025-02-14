use std::collections::HashMap;

use radsort::sort_by_key;

use crate::ruranges_structs::Event;
use crate::ruranges_structs::EventUsize;
use crate::ruranges_structs::Interval;
use crate::ruranges_structs::MaxEvent;
use crate::ruranges_structs::MinEvent;
use crate::ruranges_structs::SplicedSubsequenceInterval;
use crate::ruranges_structs::SubsequenceInterval;

pub fn build_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
) -> Vec<Interval> {
    let mut intervals: Vec<Interval> = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(Interval {
            group: chrs[i],
            start: starts[i] - slack,
            end: ends[i] + slack,
            idx: i,
        });
    }

    intervals
}

pub fn build_subsequence_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    strand_flags: &[bool],
) -> Vec<SplicedSubsequenceInterval> {
    let mut intervals: Vec<SplicedSubsequenceInterval> = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(SplicedSubsequenceInterval {
            chr: chrs[i],
            start: if strand_flags[i] {
                starts[i]
            } else {
                -starts[i]
            }, // so that negative strand intervals are sorted in the correct direction
            end: if strand_flags[i] { ends[i] } else { -ends[i] }, // we will find the absolute value when using them
            idx: i,
            forward_strand: strand_flags[i],
            temp_cumsum: 0,
            temp_length: 0,
        });
    }

    intervals
}

pub fn build_sequence_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    strand_flags: &[bool],
    force_plus_strand: bool,
) -> Vec<SubsequenceInterval> {
    let mut intervals: Vec<SubsequenceInterval> = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(SubsequenceInterval {
            group_id: chrs[i],
            start: if force_plus_strand || strand_flags[i] {
                starts[i]
            } else {
                -starts[i]
            }, // so that negative strand intervals are sorted in the correct direction
            end: if force_plus_strand || strand_flags[i] {
                ends[i]
            } else {
                -ends[i]
            }, // we will find the absolute value when using them
            idx: idxs[i],
            forward_strand: strand_flags[i],
        });
    }

    intervals
}

pub fn build_sorted_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
    sort_on_ends_too: bool,
) -> Vec<Interval> {
    let mut intervals = build_intervals(chrs, starts, ends, slack);

    if sort_on_ends_too {
        sort_by_key(&mut intervals, |i| i.end);
    };
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.group);

    intervals
}

pub fn build_sorted_subsequence_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    strand_flags: &[bool],
) -> Vec<SplicedSubsequenceInterval> {
    let mut intervals = build_subsequence_intervals(chrs, starts, ends, strand_flags);

    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.chr);

    intervals
}

pub fn build_sorted_sequence_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    strand_flags: &[bool],
    force_plus_strand: bool,
) -> Vec<SubsequenceInterval> {
    let mut intervals =
        build_sequence_intervals(chrs, starts, ends, idxs, strand_flags, force_plus_strand);

    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.group_id);

    intervals
}

pub fn sort_order_idx(chrs: &[i64], starts: &[i64], ends: &[i64]) -> Vec<usize> {
    build_sorted_intervals(chrs, starts, ends, 0, true)
        .iter()
        .map(|i| i.idx)
        .collect()
}

fn split_by_chromosome(mut intervals: Vec<Interval>) -> HashMap<i64, Vec<Interval>> {
    let mut result = HashMap::new();
    if intervals.len() == 0 {
        return result;
    }
    sort_by_key(&mut intervals, |i| i.group);

    let mut current_chr = intervals.first().unwrap().group;
    let mut current_group = Vec::new();

    for interval in intervals {
        if current_chr != interval.group {
            // We encountered a new chromosome, so store the old group
            result.insert(current_chr, std::mem::take(&mut current_group));
            current_chr = interval.group;
        } else {
            // First interval
            current_chr = interval.group;
        }
        current_group.push(interval);
    }

    // Push the last group if it exists
    if !current_group.is_empty() {
        result.insert(current_chr, current_group);
    }

    result
}

pub fn align_interval_collections_on_chromosome(
    intervals1: &mut [Interval],
    intervals2: &mut [Interval],
) -> HashMap<i64, (Vec<Interval>, Vec<Interval>)> {
    // Group each set of intervals by chromosome.
    let map1 = split_by_chromosome(intervals1.to_vec());
    let map2 = split_by_chromosome(intervals2.to_vec());

    let mut result = HashMap::new();

    // First pass: fill with everything from map1.
    // If the same chromosome also appears in map2, pair them; otherwise, second vector is empty.
    for (chr, group1) in map1 {
        // If `chr` doesn't exist in map2, we get an empty vec.
        let group2 = map2.get(&chr).cloned().unwrap_or_default();
        result.insert(chr, (group1, group2));
    }

    // Second pass: make sure we also include chromosomes that appear *only* in map2.
    // If we've already inserted them above, do nothing. Otherwise, create an entry
    // with an empty vector for the first collection, and `group2` for the second.
    for (chr, group2) in map2 {
        // Only insert if the key is missing
        result.entry(chr).or_insert_with(|| (Vec::new(), group2));
    }

    result
}

pub fn build_sorted_events_single_position(
    chrs: &[i64],
    pos: &[i64],
    start: bool,
    first_set: bool,
    negative_position: bool,
    slack: i64,
) -> Vec<Event> {
    let mut events: Vec<Event> = Vec::with_capacity(2 * (chrs.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        let pos = if start {
            pos[i] - slack
        } else {
            pos[i] + slack
        };
        events.push(Event {
            chr: chrs[i],
            pos: if negative_position { -pos } else { pos },
            is_start: start,
            first_set: first_set,
            idx: i,
        });
    }

    sort_by_key(&mut events, |e| (e.chr, e.pos, e.is_start));

    events
}

pub fn build_sorted_events_single_collection(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
) -> Vec<Event> {
    let mut events: Vec<Event> = Vec::with_capacity(2 * (chrs.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(Event {
            chr: chrs[i],
            pos: starts[i],
            is_start: true,
            first_set: true,
            idx: i,
        });
        events.push(Event {
            chr: chrs[i],
            pos: ends[i] + slack,
            is_start: false,
            first_set: true,
            idx: i,
        });
    }

    // Sort events by:
    // 1. pos (ascending)
    // 2. is_start before is_end (if pos ties)
    // (We don't strictly need to tie-break by set_id or idx, but we can.)

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_events_single_collection_separate_outputs(
    chrs: &[i64],
    pos: &[i64],
    slack: i64,
) -> Vec<MinEvent> {
    let mut out_pos: Vec<MinEvent> = Vec::with_capacity(chrs.len());

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        out_pos.push(MinEvent {
            chr: chrs[i],
            pos: pos[i] - slack,
            idx: i,
        });
    }

    sort_by_key(&mut out_pos, |e| e.pos);
    sort_by_key(&mut out_pos, |e| e.chr);

    out_pos
}

pub fn build_sorted_events_with_starts_ends(
    chrs: &[i64],
    pos: &[i64],
    slack: i64,
) -> Vec<MinEvent> {
    let mut out_pos: Vec<MinEvent> = Vec::with_capacity(chrs.len());

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        out_pos.push(MinEvent {
            chr: chrs[i],
            pos: pos[i] - slack,
            idx: i,
        });
    }

    sort_by_key(&mut out_pos, |e| e.pos);
    sort_by_key(&mut out_pos, |e| e.chr);

    out_pos
}

pub fn build_sorted_events(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> Vec<EventUsize> {
    let mut events: Vec<EventUsize> = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(EventUsize {
            chr: chrs[i],
            pos: if slack < starts[i] {
                starts[i] - slack
            } else {
                0
            },
            is_start: true,
            first_set: true,
            idx: i,
        });
        events.push(EventUsize {
            chr: chrs[i],
            pos: ends[i].saturating_add(slack),
            is_start: false,
            first_set: true,
            idx: i,
        });
    }

    for j in 0..chrs2.len() {
        events.push(EventUsize {
            chr: chrs2[j],
            pos: starts2[j],
            is_start: true,
            first_set: false,
            idx: j,
        });
        events.push(EventUsize {
            chr: chrs2[j],
            pos: ends2[j],
            is_start: false,
            first_set: false,
            idx: j,
        });
    }

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_maxevents_with_starts_ends(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> Vec<MaxEvent> {
    let mut events: Vec<MaxEvent> = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(MaxEvent {
            chr: chrs[i],
            pos: starts[i] - slack,
            start: starts[i] - slack,
            end: ends[i] + slack,
            is_start: true,
            first_set: true,
            idx: i,
        });
        events.push(MaxEvent {
            chr: chrs[i],
            pos: ends[i] + slack,
            end: ends[i] + slack,
            start: starts[i] - slack,
            is_start: false,
            first_set: true,
            idx: i,
        });
    }

    for i in 0..chrs2.len() {
        events.push(MaxEvent {
            chr: chrs2[i],
            pos: starts2[i],
            start: starts2[i],
            end: ends2[i],
            is_start: true,
            first_set: false,
            idx: i,
        });
        events.push(MaxEvent {
            chr: chrs2[i],
            pos: ends2[i],
            start: starts2[i],
            end: ends2[i],
            is_start: false,
            first_set: false,
            idx: i,
        });
    }

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_events_idxs(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> Vec<Event> {
    let mut events: Vec<Event> = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(Event {
            chr: chrs[i],
            pos: starts[i] - slack,
            is_start: true,
            first_set: true,
            idx: i,
        });
        events.push(Event {
            chr: chrs[i],
            pos: ends[i] + slack,
            is_start: false,
            first_set: true,
            idx: i,
        });
    }

    for j in 0..chrs2.len() {
        events.push(Event {
            chr: chrs2[j],
            pos: starts2[j],
            is_start: true,
            first_set: false,
            idx: j,
        });
        events.push(Event {
            chr: chrs2[j],
            pos: ends2[j],
            is_start: false,
            first_set: false,
            idx: j,
        });
    }

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_events_from_intervals(
    intervals1: &mut [Interval],
    intervals2: &mut [Interval],
) -> Vec<Event> {
    let mut events: Vec<Event> = Vec::with_capacity(2 * (intervals1.len() + intervals2.len()));

    // Convert set1 intervals into events
    for interval in intervals1 {
        events.push(Event {
            chr: interval.group,
            pos: interval.start,
            is_start: true,
            first_set: true,
            idx: interval.idx,
        });
        events.push(Event {
            chr: interval.group,
            pos: interval.end,
            is_start: false,
            first_set: true,
            idx: interval.idx,
        });
    }

    for interval in intervals2 {
        events.push(Event {
            chr: interval.group,
            pos: interval.start,
            is_start: true,
            first_set: false,
            idx: interval.idx,
        });
        events.push(Event {
            chr: interval.group,
            pos: interval.end,
            is_start: false,
            first_set: false,
            idx: interval.idx,
        });
    }

    // Sort events by:
    // 1. pos (ascending)
    // 2. is_start before is_end (if pos ties)
    // (We don't strictly need to tie-break by set_id or idx, but we can.)
    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);

    events
}
