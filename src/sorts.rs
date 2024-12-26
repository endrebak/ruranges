use radsort::sort_by_key;

use crate::ruranges_structs::Interval;
use crate::ruranges_structs::Event;

pub fn sort_intervals(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
) -> Vec<i64> {

    let mut intervals: Vec<Interval> = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(Interval {
            chr:   chrs[i],
            start: starts[i],
            end:   ends[i],
            idx:   idxs[i],
        });
    }

    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.chr);

    intervals.iter().map(|i| i.idx).collect::<Vec<_>>()
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
    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}
