/// A small struct to hold an interval plus its original index.
#[derive(Debug)]
struct Interval {
    // Chromosome or group
    group: i64,
    start: i64,
    end: i64,
    idx: i64,
}

/// Returns all cross‐overlaps between two sets of intervals,
/// enumerating every pair of overlapping intervals (iv1 in set1, iv2 in set2).
///
/// This uses an incremental two‐pointer approach:
///  1. Sort both sets by (chr, start).
///  2. For each interval in set1, advance pointer `j` in set2 so that
///     set2[j] is not completely to the left of set1[i].
///  3. From there, iterate forward to capture all overlaps on the same chromosome
///     until set2[k] starts after set1[i] ends.
///  4. Never move `j` backward; only move it forward (this avoids re-checking old intervals).
///
pub fn two_pointer_overlaps(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
    slack: i64,
) -> (Vec<i64>, Vec<i64>) {
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    if chrs.is_empty() || chrs2.is_empty() {
        return (overlaps, overlaps2);
    }

    // Build and sort intervals
    let set1 = crate::sorts::build_sorted_intervals(chrs, starts, ends, idxs, slack, false);
    let set2 = crate::sorts::build_sorted_intervals(chrs2, starts2, ends2, idxs2, 0, false);

    let mut j = 0;

    // For each interval in set1
    for i in 0..set1.len() {
        let iv1 = &set1[i];

        // 1) Advance j so that set2[j] is not completely to the left of iv1
        //    or on a smaller chromosome
        while j < set2.len() {
            let iv2 = &set2[j];
            if iv2.group < iv1.group {
                // set2[j] is on an earlier chromosome
                j += 1;
            } else if iv2.group > iv1.group {
                // set2[j] is on a later chromosome, so we won't find any overlap here
                break;
            } else {
                // same chromosome => check if set2[j] is still left of iv1
                if iv2.end < iv1.start {
                    j += 1;
                } else {
                    // iv2 is now possibly overlapping or beyond the start of iv1
                    break;
                }
            }
        }

        // 2) From j forward, collect all overlaps on the same chromosome
        let mut k = j;
        while k < set2.len() {
            let iv2 = &set2[k];
            // If we've moved on to another chromosome or iv2 starts after iv1 ends, break
            if iv2.group != iv1.group || iv2.start > iv1.end {
                break;
            }
            // If they overlap, record
            if iv1.end >= iv2.start && iv2.end >= iv1.start {
                overlaps.push(iv1.idx);
                overlaps2.push(iv2.idx);
            }
            k += 1;
        }
    }

    (overlaps, overlaps2)
}
