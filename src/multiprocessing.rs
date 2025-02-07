use std::fmt;

use crate::ruranges_structs::MinEvent;


pub fn find_chr_boundaries_minevents(data: &[MinEvent]) -> Vec<usize> {
    let mut boundaries = Vec::new();

    // Start boundary (beginning of first chromosome group)
    boundaries.push(0);

    // Identify every index `i` where the chromosome changes
    for i in 1..data.len() {
        if data[i].chr != data[i - 1].chr {
            boundaries.push(i);
        }
    }

    // Final boundary (end of the last chromosome group)
    boundaries.push(data.len());

    boundaries
}

/// Holds combined boundaries for a single chromosome across two vectors.
#[derive(Debug, Clone)]
pub struct ChrBound {
    pub chr: i64,
    pub start1: usize,
    pub end1: usize,
    pub start2: usize,
    pub end2: usize,
}

impl fmt::Display for ChrBound {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Customize the output format as desired.
        write!(f, "ChrBound {{ chr: {}, start1: {}, end1: {}, start2: {}, end2: {}, len1: {}, len2: {}, }}",
            self.chr, self.start1, self.end1, self.start2, self.end2, self.end1 - self.start1, self.end2 - self.start2)
    }
}

/// Returns boundary indices [0, ..., data.len()] whenever `chr` changes.
/// E.g. if `data` has chr=1 for indices [0..2), chr=2 for [2..5), etc.,
/// then you might get [0, 2, 5] (and finally data.len()).
fn find_chr_boundaries(data: &[MinEvent]) -> Vec<usize> {
    let mut boundaries = Vec::new();
    if data.is_empty() {
        return boundaries;
    }

    // Always push the start index 0
    boundaries.push(0);

    // Mark where the chromosome changes
    for i in 1..data.len() {
        if data[i].chr != data[i - 1].chr {
            boundaries.push(i);
        }
    }

    // Push the final end
    boundaries.push(data.len());
    boundaries
}

/// Converts boundary indices into a list of (chr, start_index, end_index) blocks.
/// Each block covers all MinEvents for a single chromosome in `data`.
fn build_chr_blocks(data: &[MinEvent], boundaries: &[usize]) -> Vec<(i64, usize, usize)> {
    let mut blocks = Vec::new();
    for w in boundaries.windows(2) {
        let start = w[0];
        let end = w[1];
        // If start < end, we have a real block
        if start < end {
            let chrom = data[start].chr;
            blocks.push((chrom, start, end));
        }
    }
    blocks
}


/// A small helper struct for the final results.
/// Each partition covers [start1..end1) in `sorted_starts`
/// and [start2..end2) in `sorted_starts2`.
#[derive(Debug)]
pub struct PartitionIndex {
    pub start1: usize,
    pub end1: usize,
    pub start2: usize,
    pub end2: usize,
}

/// A helper struct to store the range of indices for a contiguous
/// set of events on a single chromosome.
#[derive(Debug)]
struct ChromRange {
    chr: i64,
    start_idx: usize,
    end_idx: usize, // end_idx is exclusive
}

/// Given a sorted slice of MinEvents, group them by chromosome
/// and return a Vec of (chr, start_idx, end_idx).
fn group_by_chromosome(events: &[MinEvent]) -> Vec<ChromRange> {
    if events.is_empty() {
        return vec![];
    }

    let mut ranges = Vec::new();

    let mut current_chr = events[0].chr;
    let mut current_start = 0usize;

    for i in 1..events.len() {
        if events[i].chr != current_chr {
            // We've hit a new chromosome, close out the old range
            ranges.push(ChromRange {
                chr: current_chr,
                start_idx: current_start,
                end_idx: i,
            });
            // start a new range
            current_chr = events[i].chr;
            current_start = i;
        }
    }
    // close the final range
    ranges.push(ChromRange {
        chr: current_chr,
        start_idx: current_start,
        end_idx: events.len(),
    });

    ranges
}

/// Partition a single sorted slice (grouped by chromosome) into N partitions.
/// Each partition is represented as (start_index, end_index) into the original slice.
fn partition_chrom_ranges(
    events: &[MinEvent],
    num_partitions: usize,
) -> Vec<(usize, usize)> {
    if events.is_empty() {
        // If no events, return num_partitions empty partitions
        return (0..num_partitions).map(|_| (0, 0)).collect();
    }
    if num_partitions == 0 {
        return vec![];
    }

    let chrom_ranges = group_by_chromosome(events);

    // total events
    let total_len = events.len();
    let target_chunk_size = (total_len as f64 / num_partitions as f64).ceil() as usize;

    let mut partitions = Vec::with_capacity(num_partitions);
    let mut current_start = chrom_ranges[0].start_idx;
    let mut accumulated = 0; // count of events in the current partition
    let mut partition_count = 1;

    for (i, chr_range) in chrom_ranges.iter().enumerate() {
        let chr_range_size = chr_range.end_idx - chr_range.start_idx;
        let tentative_new_size = accumulated + chr_range_size;

        // If adding this chromosome range exceeds target_chunk_size
        // and we still have space for more partitions, then we close
        // the current partition before adding this chromosome.
        if partition_count < num_partitions && // we can only "cut" if we still have partitions left to form
           accumulated > 0 && // avoid zero-length partition in normal logic
           (tentative_new_size > target_chunk_size)
        {
            // close out the previous partition
            partitions.push((current_start, chr_range.start_idx));
            partition_count += 1;

            // start a new partition with this chromosome
            current_start = chr_range.start_idx;
            accumulated = 0;
        }

        accumulated += chr_range_size;

        // if this is the last chromosome or if we have formed
        // the last partition, we close it automatically
        if i == chrom_ranges.len() - 1 {
            // close out final partition
            partitions.push((current_start, chr_range.end_idx));
        }
    }

    // If we still don't have exactly num_partitions, we can pad out
    // or merge the last ones. This naive approach just merges any
    // extras at the end if we created more than needed (which can
    // happen if lots of single-chr partitions blow up).
    // In many real-world scenarios, you might do a more sophisticated
    // balancing, but here we keep it simple.
    if partitions.len() > num_partitions {
        // Merge the extra partitions into the last one
        let mut merged = partitions[0..(num_partitions - 1)].to_vec();
        // The last partition we merge everything from the leftover
        let last_start = partitions[num_partitions - 1].0;
        let last_end = partitions.last().unwrap().1;
        merged.push((last_start, last_end));
        partitions = merged;
    } else if partitions.len() < num_partitions {
        // If fewer partitions, we can just duplicate the last range as no-ops
        while partitions.len() < num_partitions {
            let last = *partitions.last().unwrap();
            partitions.push(last);
        }
    }

    partitions
}

/// Create `num_partitions` partitions for *both* slices, ensuring no chromosome boundaries
/// are crossed in either slice. Each returned element describes the start/end in slice1
/// and the start/end in slice2.
pub fn partition_two_arrays(
    sorted_starts: &[MinEvent],
    sorted_starts2: &[MinEvent],
    num_partitions: usize,
) -> Vec<PartitionIndex> {
    let parts1 = partition_chrom_ranges(sorted_starts, num_partitions);
    let parts2 = partition_chrom_ranges(sorted_starts2, num_partitions);

    // Zip them into a single vector of PartitionIndex
    parts1
        .into_iter()
        .zip(parts2.into_iter())
        .map(|((start1, end1), (start2, end2))| PartitionIndex {
            start1,
            end1,
            start2,
            end2,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_two_arrays() {
        // A small mock dataset with two chromosomes, 5 events on chr1,
        // then 4 events on chr2, for each slice.
        let ev1 = vec![
            MinEvent { chr: 1, pos: 10, idx: 0 },
            MinEvent { chr: 1, pos: 20, idx: 1 },
            MinEvent { chr: 1, pos: 30, idx: 2 },
            MinEvent { chr: 1, pos: 40, idx: 3 },
            MinEvent { chr: 1, pos: 50, idx: 4 },
            MinEvent { chr: 2, pos: 10, idx: 5 },
            MinEvent { chr: 2, pos: 20, idx: 6 },
            MinEvent { chr: 2, pos: 30, idx: 7 },
            MinEvent { chr: 2, pos: 40, idx: 8 },
        ];
        let ev2 = vec![
            MinEvent { chr: 1, pos: 15, idx: 0 },
            MinEvent { chr: 1, pos: 25, idx: 1 },
            MinEvent { chr: 1, pos: 35, idx: 2 },
            MinEvent { chr: 2, pos: 5,  idx: 3 },
            MinEvent { chr: 2, pos: 15, idx: 4 },
            MinEvent { chr: 2, pos: 25, idx: 5 },
        ];

        // Let's request 3 partitions
        let results = partition_two_arrays(&ev1, &ev2, 3);

        for (i, part) in results.iter().enumerate() {
            println!("Partition {}: {:?}", i, part);
        }

        // Here we only check that we got exactly 3 partitions:
        assert_eq!(results.len(), 3);

        // Additional checks or asserts can verify the boundaries do not cross chrs, etc.
    }
}
