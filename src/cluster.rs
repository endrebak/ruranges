use std::time::Instant;

use rustc_hash::FxHashSet;

use crate::sorts;

pub fn sweep_line_cluster(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
) -> (Vec<i64>, Vec<i64>) {
    let start = Instant::now();

    let mut indices = Vec::with_capacity(chrs.len());
    let mut cluster_ids = Vec::with_capacity(chrs.len());

    if chrs.is_empty() {
        return (cluster_ids, indices);
    };

    let events = sorts::build_sorted_events_single_collection(chrs, starts, ends, idxs);
    let duration = start.elapsed();

    let mut active1 = FxHashSet::default();

    let mut current_chr: i32 = events.first().unwrap().chr;
    let mut current_cluster: i64 = 0;

    for e in events {
        if e.chr != current_chr {
            active1.clear();
            current_chr = e.chr;
        }

        if e.is_start {
            indices.push(e.idx);
            cluster_ids.push(current_cluster);
            // Interval is starting
            // Now add it to active1
            active1.insert(e.idx);
        } else {
            // Interval is ending
            active1.remove(&e.idx);
            if active1.is_empty() {
                current_cluster += 1;
            }
        }
    }

    let duration = start.elapsed();
    println!("Time elapsed finding overlaps: {:?}", duration);

    (cluster_ids, indices)
}
