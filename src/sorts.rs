use std::collections::HashMap;
use std::time::Instant;

use radsort::sort_by_key;

use crate::ruranges_structs::Interval;
use crate::ruranges_structs::Event;


pub fn build_intervals(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
) -> Vec<Interval> {
    let mut intervals: Vec<Interval> = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(Interval {
            chr:   chrs[i],
            start: starts[i],
            end:   ends[i],
            idx:   idxs[i],
            _idx: 0,
        });
    }

    intervals
}


pub fn build_sorted_intervals(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
) -> Vec<Interval> {
    let mut intervals = build_intervals(chrs, starts, ends, idxs);

    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.chr);

    intervals
}


pub fn sort_order_idx(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
) -> Vec<i64> {
    build_sorted_intervals(chrs, starts, ends, idxs).iter().map(|i| i.idx).collect()
}


fn split_by_chromosome(mut intervals: Vec<Interval>) -> HashMap<i32, Vec<Interval>> {
    let mut result = HashMap::new();
    if intervals.len() == 0 {
        return result;
    }
    sort_by_key(&mut intervals, |i| i.chr);

    let mut k = 0;
    let mut current_chr = intervals.first().unwrap().chr;
    let mut current_group = Vec::new();

    for mut interval in intervals {
        if current_chr != interval.chr {
            // We encountered a new chromosome, so store the old group
            result.insert(current_chr, std::mem::take(&mut current_group));
            current_chr = interval.chr;
            k = 0;
        } else {
            // First interval
            current_chr = interval.chr;
        }
        interval._idx = k;
        current_group.push(interval);
        k += 1;
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
) -> HashMap<i32, (Vec<Interval>, Vec<Interval>)> {
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


pub fn build_sorted_events(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i32],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
) -> Vec<Event> {
    let mut events: Vec<Event> = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(Event {
            chr: chrs[i],
            pos: starts[i],
            is_start: true,
            first_set: true,
            idx: idxs[i],
        });
        events.push(Event {
            chr: chrs[i],
            pos: ends[i],
            is_start: false,
            first_set: true,
            idx: idxs[i],
        });
    }

    for j in 0..chrs2.len() {
        events.push(Event {
            chr: chrs2[j],
            pos: starts2[j],
            is_start: true,
            first_set: false,
            idx: idxs2[j],
        });
        events.push(Event {
            chr: chrs2[j],
            pos: ends2[j],
            is_start: false,
            first_set: false,
            idx: idxs2[j],
        });
    }

    // Sort events by:
    // 1. pos (ascending)
    // 2. is_start before is_end (if pos ties)
    // (We don't strictly need to tie-break by set_id or idx, but we can.)

    let start = Instant::now();
    sort_by_key(&mut events, |e| e.is_start);
    let duration = start.elapsed();
    println!("Time elapsed sorting events1: {:?}", duration);
    sort_by_key(&mut events, |e| e.pos);
    let duration = start.elapsed();
    println!("Time elapsed sorting events2: {:?}", duration);
    sort_by_key(&mut events, |e| e.chr);
    let duration = start.elapsed();
    println!("Time elapsed sorting events: {:?}", duration);

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
            chr: interval.chr,
            pos: interval.start,
            is_start: true,
            first_set: true,
            idx: interval.idx,
        });
        events.push(Event {
            chr: interval.chr,
            pos: interval.end,
            is_start: false,
            first_set: true,
            idx: interval.idx,
        });
    }

    for interval in intervals2 {
        events.push(Event {
            chr: interval.chr,
            pos: interval.start,
            is_start: true,
            first_set: false,
            idx: interval.idx,
        });
        events.push(Event {
            chr: interval.chr,
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
