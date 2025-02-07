use std::time::Instant;

use crate::{
    nearest_unique_k::sweep_line_k_nearest,
    overlaps::{self, sweep_line_overlaps_nearest},
    ruranges_structs::{MinEvent, Nearest, OverlapPair},
    sorts::build_sorted_events_single_collection_separate_outputs,
};

/// For each MinEvent in `sorted_ends`, find up to `k` *unique positions*
/// in `sorted_starts2` that lie to the right (including equal position on the
/// same chromosome). If multiple entries in `sorted_starts2` share the same
/// position, they all get reported, but they count as one unique position.
pub fn nearest_intervals_to_the_right(
    sorted_ends: Vec<MinEvent>,
    sorted_starts2: Vec<MinEvent>,
    k: usize,
) -> Vec<Nearest> {
    // We might need more than `sorted_ends.len()` because each end could
    // contribute up to `k` *unique positions* (potentially multiplied by the
    // number of intervals sharing those positions). So we set capacity
    // accordingly.
    // This is not strictly required, but it helps performance to reserve enough space.
    let mut output = Vec::with_capacity(sorted_ends.len().saturating_mul(k));

    let n_starts = sorted_starts2.len();

    // `j` will track our position in sorted_starts2 as we move through sorted_ends.
    let mut j = 0usize;

    // Iterate over each 'end' event
    for end in &sorted_ends {
        let end_chr = end.chr;
        let end_pos = end.pos;

        // Advance `j` so that sorted_starts2[j] is the first start
        // that is >= end_pos on the same chrom (or beyond).
        // Because both arrays are sorted, we never need to move `j` backward.
        while j < n_starts {
            let start = &sorted_starts2[j];
            if start.chr < end_chr {
                // still on a smaller chromosome; move j forward
                j += 1;
            } else if start.chr == end_chr && start.pos < end_pos {
                // same chrom but still to the left; move j forward
                j += 1;
            } else {
                // now start.chr > end_chr (i.e. next chromosome) OR
                // start.chr == end_chr && start.pos >= end_pos
                // -> we've reached a region that is "to the right" or next chrom
                break;
            }
        }

        // Now collect up to k unique positions (on the same chromosome).
        let mut unique_count = 0;
        let mut last_pos: Option<i64> = None;

        // We'll scan from `j` onward, but we do NOT move `j` itself
        // because the next 'end' might need a similar or slightly advanced position.
        // Instead, we use `local_idx` to look ahead for this specific end.
        let mut local_idx = j;
        while local_idx < n_starts {
            let start = &sorted_starts2[local_idx];

            // If we've passed beyond the chromosome of this end, we won't find
            // any more right-side intervals for this end.
            if start.chr != end_chr {
                break;
            }

            // Check if we're at a new unique position
            if last_pos.map_or(true, |lp| start.pos != lp) {
                unique_count += 1;
                if unique_count > k {
                    // we've reached the limit of k unique positions
                    break;
                }
                last_pos = Some(start.pos);
            }

            // This start is included in the results
            let distance = start.pos - end_pos + 1; // can be 0 or positive
            output.push(Nearest {
                distance,
                idx: end.idx,
                idx2: start.idx,
            });

            local_idx += 1;
        }
    }

    output
}

/// For each MinEvent in `sorted_ends`, find up to `k` *unique positions*
/// in `sorted_starts2` that lie to the left (strictly smaller position on
/// the same chromosome). If multiple entries in `sorted_starts2` share
/// the same position, they all get reported, but they count as one
/// unique position in the limit `k`.
pub fn nearest_intervals_to_the_left(
    sorted_ends: Vec<MinEvent>,
    sorted_starts2: Vec<MinEvent>,
    k: usize,
) -> Vec<Nearest> {
    // The max possible size is (number of ends) * (k + duplicates at each of those k positions).
    // We reserve a rough upper bound for efficiency.
    let mut output = Vec::with_capacity(sorted_ends.len().saturating_mul(k));

    let n_starts = sorted_starts2.len();
    let mut j = 0_usize; // Points into sorted_starts2

    for end in &sorted_ends {
        let end_chr = end.chr;
        let end_pos = end.pos;

        // Move `j` forward so that:
        // - All start events at indices < j have start.chr < end_chr
        //   OR (start.chr == end_chr && start.pos < end_pos).
        // - Equivalently, sorted_starts2[j] is the *first* event that is NOT
        //   strictly to the left of `end`.
        while j < n_starts {
            let start = &sorted_starts2[j];
            if start.chr < end_chr {
                // still a smaller chromosome => definitely to the left
                j += 1;
            } else if start.chr == end_chr && start.pos < end_pos {
                // same chrom, smaller position => to the left
                j += 1;
            } else {
                // we've reached a start that is not to the left
                break;
            }
        }

        // Now, everything in [0..j) is strictly to the left of `end`.
        // We'll look backwards from j-1 to gather up to k unique positions
        // on the same chromosome.
        if j == 0 {
            // No intervals to the left; skip
            continue;
        }

        let mut local_idx = j - 1;
        let mut unique_count = 0;
        let mut last_pos: Option<i64> = None;

        // Descend from j-1 down to 0 (or until we break).
        loop {
            let start = &sorted_starts2[local_idx];

            // Must match the same chromosome
            if start.chr != end_chr {
                break;
            }

            // Check if we have a new (unique) position
            if last_pos.map_or(true, |lp| start.pos != lp) {
                unique_count += 1;
                if unique_count > k {
                    break;
                }
                last_pos = Some(start.pos);
            }

            // Calculate the distance (end.pos - start.pos)
            // Here, start.pos < end.pos by definition if we get here.
            let distance = end_pos - start.pos + 1;
            output.push(Nearest {
                distance,
                idx: end.idx,    // the 'end' event's idx
                idx2: start.idx, // the 'start' event's idx
            });

            if local_idx == 0 {
                break;
            }
            local_idx -= 1;
        }
    }

    output
}

/// Merges three lists:
///  1) Overlaps (distance = 0),
///  2) Nearest-left (distance >= 1 in ascending order),
///  3) Nearest-right (distance >= 1 in ascending order).
///
/// We collect up to `k` *distinct* distances.  If multiple items share
/// a distance, we keep them all.  Overlaps count as distance=0 (the
/// first distinct distance), and the "nearest" entries are offset by +1
/// when emitted so that distance=1 => 2, distance=2 => 3, etc.
///
/// Returns three parallel vectors:
///  - `idx1_vec`: i64
///  - `idx2_vec`: i64
///  - `dist_vec`: i64
///
/// No extra sorting or temporary vectors are required.
pub fn merge_three_way(
    overlaps: &[OverlapPair],  // all distance=0
    nearest_left: &[Nearest],  // sorted ascending by .distance
    nearest_right: &[Nearest], // sorted ascending by .distance
    k: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // Reserve enough for all possible results (upper bound).
    let capacity = overlaps.len() + nearest_left.len() + nearest_right.len();
    let mut idx1_vec = Vec::with_capacity(capacity);
    let mut idx2_vec = Vec::with_capacity(capacity);
    let mut dist_vec = Vec::with_capacity(capacity);

    // Indices for each of the three inputs
    let (mut i, mut j, mut r) = (0_usize, 0_usize, 0_usize);

    // We'll track how many *distinct* distances we've emitted.
    let mut distinct_count = 0_usize;
    let mut last_dist: Option<i64> = None;

    loop {
        // 1) Distance from the "overlaps" side is always 0, if any remain.
        let dist_o = if i < overlaps.len() { 0 } else { i64::MAX };

        // 2) Distance from the "nearest_left" side is `nearest_left[j].distance + 1`.
        let dist_l = if j < nearest_left.len() {
            nearest_left[j].distance + 1
        } else {
            i64::MAX
        };

        // 3) Distance from the "nearest_right" side is `nearest_right[r].distance + 1`.
        let dist_r = if r < nearest_right.len() {
            nearest_right[r].distance + 1
        } else {
            i64::MAX
        };

        // Pick the minimum distance among the three
        let current_dist = dist_o.min(dist_l.min(dist_r));
        if current_dist == i64::MAX {
            // All lists are exhausted
            break;
        }

        // Check if this is a new distinct distance
        if last_dist.map_or(true, |prev| prev != current_dist) {
            distinct_count += 1;
            if distinct_count > k {
                // We've already reached k distinct distances
                break;
            }
            last_dist = Some(current_dist);
        }

        // Pull all Overlaps that match `current_dist` (which can be 0 only)
        while i < overlaps.len() && current_dist == 0 {
            idx1_vec.push(overlaps[i].idx);
            idx2_vec.push(overlaps[i].idx2);
            dist_vec.push(0);
            i += 1;
        }

        // Pull all Nearest-Left that match `current_dist`
        while j < nearest_left.len() && (nearest_left[j].distance + 1) == current_dist {
            idx1_vec.push(nearest_left[j].idx as i64);
            idx2_vec.push(nearest_left[j].idx2 as i64);
            dist_vec.push(current_dist);
            j += 1;
        }

        // Pull all Nearest-Right that match `current_dist`
        while r < nearest_right.len() && (nearest_right[r].distance + 1) == current_dist {
            idx1_vec.push(nearest_right[r].idx as i64);
            idx2_vec.push(nearest_right[r].idx2 as i64);
            dist_vec.push(current_dist);
            r += 1;
        }
    }

    (idx1_vec, idx2_vec, dist_vec)
}

pub fn nearest(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    slack: i64,
    k: usize,
    include_overlaps: bool,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let start = Instant::now();
    let (sorted_starts, sorted_ends) =
        build_sorted_events_single_collection_separate_outputs(chrs, starts, ends, idxs, slack);
    let (sorted_starts2, sorted_ends2) =
        build_sorted_events_single_collection_separate_outputs(chrs2, starts2, ends2, idxs2, 0);
    println!("sorts {:.2?}", start.elapsed());

    let overlaps = if include_overlaps {
        sweep_line_overlaps_nearest(&sorted_starts, &sorted_ends, &sorted_starts2, &sorted_ends2)
    } else {
        Vec::new()
    };
    println!("overlaps {:.2?}", start.elapsed());
    println!("overlaps {}", overlaps.len());
    let nearest_left = nearest_intervals_to_the_left(sorted_ends, sorted_starts2, k);
    println!("left {:.2?}", start.elapsed());
    println!("left {}", nearest_left.len());
    let nearest_right = nearest_intervals_to_the_right(sorted_starts, sorted_ends2, k);
    println!("right {:.2?}", start.elapsed());
    println!("right {}", nearest_right.len());

    let merged = merge_three_way(&overlaps, &nearest_left, &nearest_right, k);
    println!("merge {:.2?}", start.elapsed());
    merged
}
