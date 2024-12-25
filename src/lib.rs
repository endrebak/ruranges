use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use radsort::sort_by_key;
use rustc_hash::FxHashSet;
use std::time::Instant;
use numpy::{PyReadonlyArray1, PyArray1, IntoPyArray};


#[derive(Debug, Clone)]
struct Interval {
    chr: i32,
    start: i64,
    end: i64,
    idx: i64,
}

#[pyfunction]
pub fn sort_intervals(
    chrs: PyReadonlyArray1<i32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    py: Python,
) -> PyResult<Py<PyArray1<i64>>> {
    let chrs_slice   = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice   = ends.as_slice()?;
    let idxs_slice   = idxs.as_slice()?;

    let mut intervals: Vec<Interval> = Vec::with_capacity(chrs_slice.len());
    let start = Instant::now();
    for i in 0..chrs_slice.len() {
        intervals.push(Interval {
            chr:   chrs_slice[i],
            start: starts_slice[i],
            end:   ends_slice[i],
            idx:   idxs_slice[i],
        });
    }
    let duration_first = start.elapsed();
    println!("`first_function()` took: {:?}", duration_first);

    let start = Instant::now();
    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.chr);
    let duration_first = start.elapsed();
    println!("`second_function()` took: {:?}", duration_first);

    let start = Instant::now();
    let indexes = intervals.iter().map(|i| i.idx).collect::<Vec<_>>();
    let duration_first = start.elapsed();
    println!("`third_function()` took: {:?}", duration_first);

    // `into_pyarray(py)` creates a borrowed &PyArray1<i64> internally,
    // but we need an owned Py<PyArray1<i64>> to return from the function safely.
    let start = Instant::now();
    let res = indexes.into_pyarray(py).to_owned().into();
    let duration_first = start.elapsed();
    println!("`fourth_function()` took: {:?}", duration_first);
    Ok(res)
}

/// An "event" in the sweep line:
/// - `pos`: the coordinate (start or end of an interval)
/// - `is_start`: true if it's a start event, false if it's an end event
/// - `set_id`: which set does this interval belong to? (1 or 2)
/// - `idx`: the interval's ID/index
#[derive(Debug, Clone)]
struct Event {
    chr: i32,
    pos: i64,
    is_start: bool,
    first_set: bool,
    idx: i64,
}

fn build_sorted_events(
    chrs: &[i32],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    chrs2: &[i32],
    starts2: &[i64],
    ends2: &[i64],
    idxs2: &[i64],
) -> Vec<Event> {
    let mut events: Vec<Event> = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(Event {
            chr: chrs[i],
            pos: starts[i],
            is_start: true,
            first_set: true,
            idx: idxs[i],
        });
        events.push(Event {
            chr: chrs[i],
            pos: ends[i],
            is_start: false,
            first_set: true,
            idx: idxs[i],
        });
    }

    for j in 0..chrs2.len() {
        events.push(Event {
            chr: chrs2[j],
            pos: starts2[j],
            is_start: true,
            first_set: false,
            idx: idxs2[j],
        });
        events.push(Event {
            chr: chrs2[j],
            pos: ends2[j],
            is_start: false,
            first_set: false,
            idx: idxs2[j],
        });
    }

    // Sort events by:
    // 1. pos (ascending)
    // 2. is_start before is_end (if pos ties)
    // (We don't strictly need to tie-break by set_id or idx, but we can.)
    let start = Instant::now();
    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);
    let duration_first = start.elapsed();
    println!("`first_function()` took: {:?}", duration_first);

    events
}

// use std::collections::HashSet;

/// Returns all overlapping pairs (idx1, idx2) between intervals in set1 and set2.
/// This uses a line-sweep / active-set approach.
///
/// Algorithm steps:
///   1. Build a list of events (start & end) for each interval in both sets.
///   2. Sort events by coordinate. Where coordinates tie, put start before end.
///   3. Maintain active sets (for set1 and set2). For a start event in set1,
///      record overlap with all active in set2, then insert into active1. Etc.
///   4. Return the list of all cross-set overlaps.
fn sweep_line_overlaps(
    events: Vec<Event>,
) -> (Vec<i64>, Vec<i64>) {

    // Active sets
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    // We'll collect all cross overlaps here
    let mut overlaps = Vec::new();
    let mut overlaps2 = Vec::new();

    // Process events in ascending order of position
    for e in events {
        if e.is_start {
            // Interval is starting
            if e.first_set {
                    // Overlaps with all currently active intervals in set2
                    for &idx2 in active2.iter() {
                        overlaps.push(e.idx);
                        overlaps2.push(idx2);
                    }
                    // Now add it to active1
                    active1.insert(e.idx);
                } else
                {
                    // Overlaps with all currently active intervals in set1
                    for &idx1 in active1.iter() {
                        overlaps.push(idx1);
                        overlaps2.push(e.idx);
                    }
                    // Now add it to active2
                    active2.insert(e.idx);
                }
        } else {
            // Interval is ending
            if e.first_set {
                    active1.remove(&e.idx);
                } else {
                    active2.remove(&e.idx);
                }
            }
        }

    (overlaps, overlaps2)
}


#[pyfunction]
pub fn chromsweep(
    py: Python,
    chrs: PyReadonlyArray1<i32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice   = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice   = ends.as_slice()?;
    let idxs_slice   = idxs.as_slice()?;
    let chrs_slice2   = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2   = ends2.as_slice()?;
    let idxs_slice2   = idxs2.as_slice()?;

    let events = build_sorted_events(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        idxs_slice2,
    );
    let result = sweep_line_overlaps(events);
    Ok(
        (
        result.0.into_pyarray(py).to_owned().into(),
        result.1.into_pyarray(py).to_owned().into(),
        )
    )
}

/// The Python module definition
#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Notice we pass `py` here, not `m`
    // m.add_function(wrap_pyfunction!(build_sorted_intervals, m)?)?;
    m.add_function(wrap_pyfunction!(chromsweep, m)?)?;
    m.add_function(wrap_pyfunction!(sort_intervals, m)?)?;
    Ok(())
}