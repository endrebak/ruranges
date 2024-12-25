use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use radsort::sort_by_key;

/// The struct we want to create and sort
#[derive(Clone, Debug)]
struct Interval {
    chr: i32,
    start: i64,
    end: i64,
}

#[pyfunction]
fn build_sorted_intervals(
    chrs: Vec<i32>,
    starts: Vec<i64>,
    ends: Vec<i64>,
) -> PyResult<Vec<(i32, i64, i64)>> {
    if chrs.len() != starts.len() || starts.len() != ends.len() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "All input vectors must have the same length",
        ));
    }

    // Build a list of Interval structs
    let mut intervals = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(Interval {
            chr: chrs[i],
            start: starts[i],
            end: ends[i],
        });
    }

    // Sort by chr, then by start, then by end
    sort_by_key(&mut intervals, |ival| ival.end);
    sort_by_key(&mut intervals, |ival| ival.start);
    sort_by_key(&mut intervals, |ival| ival.chr);

    // Convert sorted Interval structs to list of tuples
    let result = intervals
        .into_iter()
        .map(|ival| (ival.chr, ival.start, ival.end))
        .collect();

    Ok(result)
}

/// The Python module definition
#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Notice we pass `py` here, not `m`
    m.add_function(wrap_pyfunction!(build_sorted_intervals, m)?)?;
    Ok(())
}