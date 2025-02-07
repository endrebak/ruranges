#[derive(Debug, Clone)]
pub struct Interval {
    pub group: i64,
    pub start: i64,
    pub end: i64,
    pub idx: i64,
}

/// An "event" in the sweep line:
/// - `pos`: the coordinate (start or end of an interval)
/// - `is_start`: true if it's a start event, false if it's an end event
/// - `set_id`: which set does this interval belong to? (1 or 2)
/// - `idx`: the interval's ID/index
#[derive(Debug, Clone, Hash)]
pub struct Event {
    pub chr: i64,
    pub pos: i64,
    pub is_start: bool,
    pub first_set: bool,
    pub idx: i64,
}

#[derive(Debug, Clone, Hash)]
pub struct MinEvent {
    pub chr: i64,
    pub pos: i64,
    pub idx: i64,
}

#[derive(Debug, Clone, Hash)]
pub struct OverlapPair {
    pub idx: i64,
    pub idx2: i64,
}

#[derive(Debug, Clone, Hash)]
pub struct Nearest {
    pub distance: i64,
    pub idx: i64,
    pub idx2: i64,
}

#[derive(Debug, Clone)]
pub struct SplicedSubsequenceInterval {
    /// Encoded chromosome (or chrom+strand+gene) ID.
    pub chr: i64,

    /// The genomic start coordinate.
    pub start: i64,

    /// The genomic end coordinate.
    pub end: i64,

    pub idx: i64,

    pub forward_strand: bool,

    /// Temporary: length = (end - start).
    pub temp_length: i64,

    /// Temporary: cumulative sum of lengths within this chrom group.
    pub temp_cumsum: i64,
}

/// A simple struct to hold each interval's data for "subsequence" logic.
#[derive(Clone)]
pub struct SubsequenceInterval {
    pub group_id: i64,        // grouping ID
    pub start: i64,           // genomic start
    pub end: i64,             // genomic end
    pub idx: i64,             // e.g. row index or something else
    pub forward_strand: bool, // true => + strand, false => - strand
}
