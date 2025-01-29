use rustc_hash::FxHashMap;

use crate::sorts;

pub fn sweep_line_subtract(
    chrs1: &[i64],
    starts1: &[i64],
    ends1: &[i64],
    idxs1: &[i64],
    chrs2: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // If either set is empty, set1 is unchanged (or trivially subtracted).
    if chrs1.is_empty() || chrs2.is_empty() {
        return (
            idxs1.to_vec(),
            starts1.to_vec(),
            ends1.to_vec(),
        );
    }

    // Build sorted events
    let events = sorts::build_sorted_events(chrs1, starts1, ends1, idxs1, chrs2, starts2, ends2, idxs2, 0);

    // Output buffers
    let mut result_idxs = Vec::new();
    let mut result_starts = Vec::new();
    let mut result_ends = Vec::new();

    // Track how many set2 intervals are active
    let mut active2_count: i64 = 0;

    // For each active interval in set1, store the position at which
    // we last started a "valid" sub-interval (when active2_count == 0).
    // i.e. active1[idx] = Some(position) means we are currently capturing
    // a sub-interval for that idx that started at `position`.
    let mut active1: FxHashMap<i64, Option<i64>> = FxHashMap::default();

    let mut current_chr = events.first().unwrap().chr;

    // We'll sweep in ascending order
    for e in events.iter() {
        // If we jumped to a new chromosome, close out everything
        // because intervals do not cross chromosome boundaries.
        if e.chr != current_chr {
            // for any active sub-interval in the old chromosome, we close them at the last event pos
            // but in typical coordinate intervals, they should already be ended by the end event. 
            // We'll do a final cleanup if you want. Usually, each interval on the old chr
            // has presumably ended with an event, but if not, you can decide to finalize them here.

            // Clear everything
            active1.clear();
            active2_count = 0;
            current_chr = e.chr;
        }

        let pos = e.pos;

        // --- 1. If we have *just arrived* at a new position, and `active2_count == 0`,
        // we are "continuing" sub-intervals for all active1. 
        //
        // But typically, the actual writing out of intervals
        // occurs at the event boundaries (start or end). 
        // We'll handle that logic around the transitions.

        // --- 2. Now handle the event itself:

        if e.first_set {
            // This event is from set1
            if e.is_start {
                // A set1 interval starts
                // If we are outside set2 (active2_count==0), 
                // that means we can immediately start capturing a sub-interval.
                if active2_count == 0 {
                    active1.insert(e.idx, Some(pos));
                } else {
                    // set2 is active, so we do not start capturing yet
                    active1.insert(e.idx, None);
                }
            } else {
                // A set1 interval ends
                // If we have been capturing a sub-interval for this idx, close it
                if let Some(start_pos) = active1.get(&e.idx).cloned().unwrap_or(None) {
                    // We are capturing. End the sub-interval at e.pos
                    if start_pos < pos {
                        result_idxs.push(e.idx);
                        result_starts.push(start_pos);
                        result_ends.push(pos);
                    }
                }
                // Remove it from active1
                active1.remove(&e.idx);
            }
        } else {
            // This event is from set2
            if e.is_start {
                // set2 interval starts
                active2_count += 1;

                // If we just went from 0 -> 1, that means we need to close
                // *all currently capturing intervals in set1* right at this boundary.
                if active2_count == 1 {
                    // close everyone
                    for (&idx1, &maybe_start) in active1.iter() {
                        if let Some(start_pos) = maybe_start {
                            // Close at current event pos (exclusive or inclusive depends on your semantics)
                            if start_pos < pos {
                                result_idxs.push(idx1);
                                result_starts.push(start_pos);
                                result_ends.push(pos);
                            }
                        }
                    }
                    // Now, set them all to None, since we cannot capture while set2 is active
                    for v in active1.values_mut() {
                        *v = None;
                    }
                }
            } else {
                // set2 interval ends
                active2_count -= 1;

                // If we just went from 1 -> 0, that means we can *resume capturing*
                // for all set1 intervals that are still active.
                if active2_count == 0 {
                    // For every set1 interval that is active, we set the start to the boundary
                    // so we resume capturing at e.pos
                    for (_idx1, v) in active1.iter_mut() {
                        if v.is_none() {
                            *v = Some(pos);
                        }
                    }
                }
            }
        }

        // Optionally, you can look ahead to the next event's position 
        // to handle the "between events" region if needed. 
        // But typically, the creation of sub-intervals at boundaries is enough.

        // 3. Move on to the next event
    }

    // No final cleanup is strictly necessary if every set1 interval has a corresponding end event.

    (result_idxs, result_starts, result_ends)
}