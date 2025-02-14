use rustc_hash::FxHashMap;

use crate::{ruranges_structs::Event, sorts};

pub fn sweep_line_complement(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
    chrom_lens: &FxHashMap<i64, i64>,
    include_first_interval: bool, // <-- new parameter
) -> (Vec<i64>, Vec<i64>, Vec<i64>, Vec<usize>) {
    let mut out_chrs = Vec::with_capacity(chrs.len());
    let mut out_starts = Vec::with_capacity(chrs.len());
    let mut out_ends = Vec::with_capacity(chrs.len());
    let mut out_idxs = Vec::with_capacity(chrs.len());

    // Early return if no input
    if chrs.is_empty() {
        return (out_chrs, out_starts, out_ends, out_idxs);
    }

    // Build your events array, sorted by chr and pos
    let events: Vec<Event> =
        sorts::build_sorted_events_single_collection(chrs, starts, ends, slack);

    // Initialize
    let mut current_chr = events[0].chr;
    let mut active_count = 0_i64;
    // Whether we start "in a hole" (i.e., complement) depends on `include_first_interval`
    let mut in_complement = include_first_interval;
    // Start the first hole at position 0 of the chromosome (only matters if `in_complement == true`)
    let mut current_start = 0_i64;
    let mut current_index = 0_usize;

    for e in events {
        // If we hit a new chromosome
        if e.chr != current_chr {
            // If we ended the previous chromosome still in a hole,
            // optionally close it out at the chromosome’s end
            if let Some(chlen) = chrom_lens.get(&current_chr) {
                if in_complement {
                    out_chrs.push(current_chr);
                    out_starts.push(current_start);
                    out_ends.push(*chlen);
                    out_idxs.push(current_index);
                }
            }

            // Reset for new chromosome
            current_chr = e.chr;
            active_count = 0;
            in_complement = include_first_interval;
            current_start = 0;
            current_index = e.idx;
        }

        // Process this event
        if e.is_start {
            // coverage X → X + 1
            active_count += 1;
            // If coverage was zero, we just ended a hole
            if active_count == 1 && in_complement && current_start != e.pos {
                // That hole ends at e.pos
                out_chrs.push(current_chr);
                out_starts.push(current_start);
                out_ends.push(e.pos);
                out_idxs.push(current_index);

                // We're no longer in a hole
                in_complement = false;
            }
        } else {
            // coverage X → X - 1
            active_count -= 1;
            // If coverage has just dropped back to zero,
            // we start a new hole here
            if active_count == 0 {
                in_complement = true;
                current_start = e.pos;
            }
        }
    }

    // End of all events: if we finished in a hole and have chromosome lengths
    if let Some(chlen) = chrom_lens.get(&current_chr) {
        if in_complement {
            out_chrs.push(current_chr);
            out_starts.push(current_start);
            out_ends.push(*chlen);
            out_idxs.push(current_index);
        }
    }

    (out_chrs, out_starts, out_ends, out_idxs)
}
