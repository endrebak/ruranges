
use std::time::Instant;

use rustc_hash::FxHashSet;

use crate::sorts;

pub fn sweep_line_merge(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let start = Instant::now();

    let mut out_chrs = Vec::with_capacity(chrs.len());
    let mut out_starts = Vec::with_capacity(chrs.len());
    let mut out_ends = Vec::with_capacity(chrs.len());

    if chrs.is_empty() {
        return (out_chrs, out_starts, out_ends);
    };

    let events = sorts::build_sorted_events_single_collection(chrs, starts, ends, idxs);

    let mut active = FxHashSet::default();

    let mut current_chr: i64 = events.first().unwrap().chr;
    let mut current_start: i64 = 0;

    for e in events {
        if e.chr != current_chr {
            active.clear();
            current_chr = e.chr;
        }

        if active.is_empty() {
            current_start = e.pos;
        }

        if e.is_start {
            active.insert(e.idx);
        } else {
            active.remove(&e.idx);
            if active.is_empty() {
                out_chrs.push(e.chr);
                out_starts.push(current_start);
                out_ends.push(e.pos);
            }
        }
    }

    let duration = start.elapsed();
    println!("Time elapsed merging: {:?}", duration);

    (out_chrs, out_starts, out_ends)
}
