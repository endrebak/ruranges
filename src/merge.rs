
use std::time::Instant;

use crate::sorts;

pub fn sweep_line_merge(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    slack: i64,
) -> (Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>) {
    let start = Instant::now();

    let mut out_indices = Vec::with_capacity(chrs.len());
    let mut out_starts = Vec::with_capacity(chrs.len());
    let mut out_ends = Vec::with_capacity(chrs.len());
    let mut counts = Vec::with_capacity(chrs.len());

    if chrs.is_empty() {
        return (out_indices, out_starts, out_ends, counts);
    };

    let events = sorts::build_sorted_events_single_collection(chrs, starts, ends, idxs, slack);

    let mut current_chr: i64 = events.first().unwrap().chr;
    let mut current_start: i64 = 0;
    let mut active_count: i64 = 0;
    let mut current_cluster_count: i64 = 0;

    for e in events {
        if e.chr != current_chr {
            active_count = 0;
            current_cluster_count = 0;
            current_chr = e.chr;
        }

        if active_count == 0 {
            current_start = e.pos;
            current_cluster_count = 0;
        }

        if e.is_start {
            active_count += 1;
            current_cluster_count += 1;
        } else {
            active_count -= 1;
            if active_count == 0 {
                out_indices.push(e.idx);
                out_starts.push(current_start);
                out_ends.push(e.pos);
                counts.push(current_cluster_count);
            }
        }
    }

    let duration = start.elapsed();
    println!("Time elapsed merging: {:?}", duration);

    (out_indices, out_starts, out_ends, counts)
}
