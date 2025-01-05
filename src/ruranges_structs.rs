#[derive(Debug, Clone)]
pub struct Interval {
    pub chr: i32,
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
    pub chr: i32,
    pub pos: i64,
    pub is_start: bool,
    pub first_set: bool,
    pub idx: i64,
}

#[derive(Debug, Clone, Hash)]
pub struct Nearest {
    pub distance: i64,
    pub idx: i64,
}
