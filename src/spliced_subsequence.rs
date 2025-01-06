use crate::{
    ruranges_structs::SplicedSubsequenceInterval,
    sorts::build_sorted_subsequence_intervals,
};

/// Replicates the "spliced_subseq" logic in one pass for intervals sorted by (chrom, start, end).
///
/// - `chrs`: chromosome (encoded) array
/// - `starts`: genomic start coordinates
/// - `ends`: genomic end coordinates
/// - `idxs`: some index array
/// - `strand_flags`: array of bool indicating forward strand (true) or reverse (false)
/// - `start`: spliced start (can be negative, to count from the 3' end)
/// - `end`: spliced end (can be None = unlimited, or negative => 3' offset)
/// - `force_plus_strand`: if `true`, treat **all** intervals as if forward strand
///
/// Returns tuple of (out_idxs, out_starts, out_ends).
pub fn spliced_subseq(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    strand_flags: &[bool],
    start: i64,
    end: Option<i64>,
    force_plus_strand: bool,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    // If no intervals, just return.
    if chrs.is_empty() {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    // Build the vector of intervals, which is already sorted by (chr, start, end) in your code.
    let intervals: Vec<SplicedSubsequenceInterval> =
        build_sorted_subsequence_intervals(chrs, starts, ends, idxs, strand_flags);

    // We'll accumulate the results here.
    let mut out_idxs: Vec<i64> = Vec::with_capacity(intervals.len());
    let mut out_starts: Vec<i64> = Vec::with_capacity(intervals.len());
    let mut out_ends: Vec<i64> = Vec::with_capacity(intervals.len());

    // A small buffer for intervals belonging to the "current chrom."
    let mut group_buf: Vec<SplicedSubsequenceInterval> = Vec::new();

    // Keep track of the current chrom and running cumsum across intervals with that chrom.
    let mut current_chrom = intervals[0].chr;
    let mut running_sum = 0_i64;

    // This closure finalizes one chrom-group: it applies negative indexing, forward/reverse logic,
    // filters out intervals with `start >= end`, and pushes results into output vectors.
    let finalize_group = |group: &mut [SplicedSubsequenceInterval],
                          start: i64,
                          end: Option<i64>,
                          force_plus: bool,
                          out_idxs: &mut Vec<i64>,
                          out_starts: &mut Vec<i64>,
                          out_ends: &mut Vec<i64>| {
        if group.is_empty() {
            return;
        }

        // The total length of this chrom group is the cumsum of the last interval in the group.
        let total_length = group.last().unwrap().temp_cumsum;

        // If end is None, use `total_length`.
        let end_val = end.unwrap_or(total_length);

        // Convert negative start/end to their positive equivalents from the 3' end:
        let global_start = if start < 0 {
            total_length + start
        } else {
            start
        };
        let global_end = if end_val < 0 {
            total_length + end_val
        } else {
            end_val
        };

        // Adjust each interval according to the spliced subsequence logic.
        for iv in group.iter_mut() {
            let cumsum_start = iv.temp_cumsum - iv.temp_length; // spliced start of this exon
            let cumsum_end = iv.temp_cumsum;                    // spliced end of this exon

            // Determine if we use forward or reverse logic:
            // if `force_plus == true`, always do forward logic
            // otherwise, check iv.forward_strand
            let is_forward = force_plus || iv.forward_strand;

            if is_forward {
                // Forward logic:
                //   start_adjust = global_start - cumsum_start
                //   if start_adjust > 0 => shift iv.start right
                //   end_adjust   = cumsum_end - global_end
                //   if end_adjust > 0 => shift iv.end left
                let start_adjust = global_start - cumsum_start;
                if start_adjust > 0 {
                    iv.start += start_adjust;
                }

                let end_adjust = cumsum_end - global_end;
                if end_adjust > 0 {
                    iv.end -= end_adjust;
                }
            } else {
                // Reverse logic (mirroring your Python snippet for strand == '-'):
                //
                //   start_adjust = global_start - cumsum_start
                //   if start_adjust > 0 => shift iv.end left
                //   end_adjust   = cumsum_end - global_end
                //   if end_adjust > 0 => shift iv.start right
                let start_adjust = global_start - cumsum_start;
                if start_adjust > 0 {
                    iv.end -= start_adjust;
                }

                let end_adjust = cumsum_end - global_end;
                if end_adjust > 0 {
                    iv.start += end_adjust;
                }
            }
        }

        let strand = group.first().unwrap().forward_strand;
        // Filter out intervals where start >= end, then push to results
        if !strand { group.reverse() };
        for iv in group.iter() {
            if iv.start < iv.end {
                out_idxs.push(iv.idx);
                out_starts.push(iv.start);
                out_ends.push(iv.end);
            }
        }
    };

    // Single pass over all intervals
    for mut interval in intervals {
        interval.start = interval.start.abs();
        interval.end = interval.end.abs();
        // If we've moved to a new chrom, finalize the previous group, then start fresh.
        if interval.chr != current_chrom {
            finalize_group(
                &mut group_buf,
                start,
                end,
                force_plus_strand,
                &mut out_idxs,
                &mut out_starts,
                &mut out_ends,
            );
            group_buf.clear();
            running_sum = 0;
            current_chrom = interval.chr;
        }

        // Prepare a copy of the interval with temp_length + cumsum.
        let mut iv_copy = interval.clone();
        iv_copy.temp_length = iv_copy.end - iv_copy.start;
        iv_copy.temp_cumsum = running_sum + iv_copy.temp_length;
        running_sum = iv_copy.temp_cumsum;

        group_buf.push(iv_copy);
    }

    // Finalize the last group
    finalize_group(
        &mut group_buf,
        start,
        end,
        force_plus_strand,
        &mut out_idxs,
        &mut out_starts,
        &mut out_ends,
    );

    (out_idxs, out_starts, out_ends)
}
