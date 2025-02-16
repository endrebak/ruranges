use std::{str::FromStr, time::Instant};

use crate::{
    overlaps::{self, sweep_line_overlaps, sweep_line_overlaps_overlap_pair},
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

/// Merges th
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Backward,
    Any,
}

impl FromStr for Direction {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "forward" => Ok(Direction::Forward),
            "backward" => Ok(Direction::Backward),
            "any" => Ok(Direction::Any),
            _ => Err("Invalid direction string"),
        }
    }
}

pub fn nearest(
    chrs: &[u32],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[u32],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
    k: usize,
    include_overlaps: bool,
    direction: &str,
) -> (Vec<u32>, Vec<u32>, Vec<i64>) {
    let dir = Direction::from_str(direction).unwrap();

    let sorted_starts = build_sorted_events_single_collection_separate_outputs(chrs, starts, slack);
    let sorted_ends = build_sorted_events_single_collection_separate_outputs(chrs, ends, slack);

    let sorted_starts2 = build_sorted_events_single_collection_separate_outputs(chrs2, starts2, 0);
    let sorted_ends2 = build_sorted_events_single_collection_separate_outputs(chrs2, ends2, 0);

    let overlaps = if include_overlaps {
        sweep_line_overlaps_overlap_pair(
            &sorted_starts,
            &sorted_ends,
            &sorted_starts2,
            &sorted_ends2,
        )
    } else {
        Vec::new()
    };
    let nearest_left = if dir == Direction::Backward || dir == Direction::Any {
        let mut tmp = nearest_intervals_to_the_left(sorted_starts, sorted_ends2, k);
        radsort::sort_by_key(&mut tmp, |n| (n.idx, n.distance));
        tmp
    } else {
        Vec::new()
    };
    let nearest_right = if dir == Direction::Forward || dir == Direction::Any {
        let mut tmp = nearest_intervals_to_the_right(sorted_ends, sorted_starts2, k);
        radsort::sort_by_key(&mut tmp, |n| (n.idx, n.distance));
        tmp
    } else {
        Vec::new()
    };

    let merged = merge_three_way_by_index_distance(&overlaps, &nearest_left, &nearest_right, k);
    merged
}

/// Merges three sources of intervals, grouped by `idx` (i.e. `idx1` in overlaps).
/// For each unique `idx`, it returns up to `k` *distinct* distances (including
/// all intervals at those distances). Overlaps are treated as distance=0 (or 1).
///
/// The data is assumed to be sorted in ascending order by `(idx, distance)`.
pub fn merge_three_way_by_index_distance(
    overlaps: &[OverlapPair],  // sorted by idx1
    nearest_left: &[Nearest],  // sorted by (idx, distance)
    nearest_right: &[Nearest], // sorted by (idx, distance)
    k: usize,
) -> (Vec<u32>, Vec<u32>, Vec<i64>) {
    // We'll return tuples: (idx, idx2, distance).
    // You can adapt if you want a custom struct instead.
    let mut idxs1 = Vec::new();
    let mut idxs2 = Vec::new();
    let mut distance = Vec::new();

    // Pointers over each input
    let (mut i, mut j, mut r) = (0_usize, 0_usize, 0_usize);

    // Outer loop: pick the smallest index among the three lists
    while i < overlaps.len() || j < nearest_left.len() || r < nearest_right.len() {
        // Current index (None if that list is exhausted)
        let idx_o = overlaps.get(i).map(|o| o.idx);
        let idx_l = nearest_left.get(j).map(|n| n.idx);
        let idx_r = nearest_right.get(r).map(|n| n.idx);

        // If all three are None, we're done
        let current_idx = match (idx_o, idx_l, idx_r) {
            (None, None, None) => break,
            (Some(a), Some(b), Some(c)) => a.min(b.min(c)),
            (Some(a), Some(b), None) => a.min(b),
            (Some(a), None, Some(c)) => a.min(c),
            (None, Some(b), Some(c)) => b.min(c),
            (Some(a), None, None) => a,
            (None, Some(b), None) => b,
            (None, None, Some(c)) => c,
        };

        // Gather all overlaps for current_idx
        let i_start = i;
        while i < overlaps.len() && overlaps[i].idx == current_idx {
            i += 1;
        }
        let overlaps_slice = &overlaps[i_start..i];

        // Gather all nearest_left for current_idx
        let j_start = j;
        while j < nearest_left.len() && nearest_left[j].idx == current_idx {
            j += 1;
        }
        let left_slice = &nearest_left[j_start..j];

        // Gather all nearest_right for current_idx
        let r_start = r;
        while r < nearest_right.len() && nearest_right[r].idx == current_idx {
            r += 1;
        }
        let right_slice = &nearest_right[r_start..r];

        // Now we have three *already-sorted* slices (by distance) for this index:
        //  1) overlaps_slice (distance=0 or 1, or if you store it in OverlapPair, read it)
        //  2) left_slice (sorted ascending by distance)
        //  3) right_slice (sorted ascending by distance)
        //
        // We'll do a 3-way merge *by distance*, collecting up to k *distinct* distances.
        // If you store overlap distances in OverlapPair, you can read them;
        // otherwise, assume overlap distance=0.

        let mut used_distances = std::collections::HashSet::new();
        let mut distinct_count = 0;

        let (mut oi, mut lj, mut rr) = (0, 0, 0);

        // Helper closures to peek distance from each slice
        let mut overlap_dist = |ix: usize| -> i64 {
            // If you store distance in OverlapPair, return that. Otherwise 0 or 1.
            // For the example, let's assume actual Overlap distance=0:
            0
        };
        let mut left_dist = |ix: usize| -> i64 { left_slice[ix].distance };
        let mut right_dist = |ix: usize| -> i64 { right_slice[ix].distance };

        // Inner loop: pick the next *smallest* distance among the three slices
        while oi < overlaps_slice.len() || lj < left_slice.len() || rr < right_slice.len() {
            // Peek next distance (or i64::MAX if none)
            let d_o = if oi < overlaps_slice.len() {
                overlap_dist(oi)
            } else {
                i64::MAX
            };
            let d_l = if lj < left_slice.len() {
                left_dist(lj)
            } else {
                i64::MAX
            };
            let d_r = if rr < right_slice.len() {
                right_dist(rr)
            } else {
                i64::MAX
            };

            let smallest = d_o.min(d_l.min(d_r));
            if smallest == i64::MAX {
                // no more items
                break;
            }

            // We'll pull everything from Overlaps that has distance == smallest
            while oi < overlaps_slice.len() {
                let dcur = overlap_dist(oi);
                if dcur == smallest {
                    // If this is a *new* distance (not in used_distances),
                    // we check if it would exceed k distinct distances
                    if !used_distances.contains(&dcur) {
                        distinct_count += 1;
                        if distinct_count > k {
                            // no new distances allowed
                            break;
                        }
                        used_distances.insert(dcur);
                    }
                    // Add to result
                    let OverlapPair { idx, idx2 } = overlaps_slice[oi];
                    idxs1.push(idx);
                    idxs2.push(idx2);
                    distance.push(dcur);
                    oi += 1;
                } else {
                    break;
                }
            }
            if distinct_count > k {
                break;
            }

            // Pull everything from Left that has distance == smallest
            while lj < left_slice.len() {
                let dcur = left_dist(lj);
                if dcur == smallest {
                    if !used_distances.contains(&dcur) {
                        distinct_count += 1;
                        if distinct_count > k {
                            break;
                        }
                        used_distances.insert(dcur);
                    }
                    let item = &left_slice[lj];
                    idxs1.push(item.idx);
                    idxs2.push(item.idx2);
                    distance.push(dcur);
                    lj += 1;
                } else {
                    break;
                }
            }
            if distinct_count > k {
                break;
            }

            // Pull everything from Right that has distance == smallest
            while rr < right_slice.len() {
                let dcur = right_dist(rr);
                if dcur == smallest {
                    if !used_distances.contains(&dcur) {
                        distinct_count += 1;
                        if distinct_count > k {
                            break;
                        }
                        used_distances.insert(dcur);
                    }
                    let item = &right_slice[rr];
                    idxs1.push(item.idx);
                    idxs2.push(item.idx2);
                    distance.push(dcur);
                    rr += 1;
                } else {
                    break;
                }
            }
            if distinct_count > k {
                break;
            }
        }
        // done collecting up to k distinct distances for this index
    }

    (idxs1, idxs2, distance)
}
