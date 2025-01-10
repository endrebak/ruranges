use std::time::Instant;


use crate::sorts;

pub fn sweep_line_cluster(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    slack: i64,
) -> (Vec<i64>, Vec<i64>) {
    let start = Instant::now();

    let mut indices = Vec::with_capacity(chrs.len());
    let mut cluster_ids = Vec::with_capacity(chrs.len());

    if chrs.is_empty() {
        return (cluster_ids, indices);
    };

    let events = sorts::build_sorted_events_single_collection(chrs, starts, ends, idxs, slack);

    let mut current_chr: i64 = events.first().unwrap().chr;
    let mut current_cluster: i64 = 0;
    let mut active_intervals: i64 = 0;

    for e in events {
        if e.chr != current_chr {
            current_cluster += 1;
            active_intervals = 0;
            current_chr = e.chr;
        }

        if e.is_start {
            indices.push(e.idx);
            cluster_ids.push(current_cluster);
            active_intervals += 1;
        } else {
            active_intervals -= 1;
            if active_intervals == 0 {
                current_cluster += 1;
            }
        }
    }

    let duration = start.elapsed();
    println!("Time elapsed finding clusters: {:?}", duration);

    (cluster_ids, indices)
}
