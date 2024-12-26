use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use numpy::{PyReadonlyArray1, PyArray1, IntoPyArray};

use crate::sorts;
use crate::overlaps;

#[pyfunction]
pub fn chromsweep_numpy(
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

    let events = sorts::build_sorted_events(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        idxs_slice2,
    );
    let result = overlaps::sweep_line_overlaps(events);
    Ok(
        (
        result.0.into_pyarray(py).to_owned().into(),
        result.1.into_pyarray(py).to_owned().into(),
        )
    )
}

#[pyfunction]
pub fn sort_intervals_numpy(
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

    let indexes = sorts::sort_intervals(chrs_slice, starts_slice, ends_slice, idxs_slice);
    Ok(indexes.into_pyarray(py).to_owned().into())
}

#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Notice we pass `py` here, not `m`
    // m.add_function(wrap_pyfunction!(build_sorted_intervals, m)?)?;
    m.add_function(wrap_pyfunction!(chromsweep_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(sort_intervals_numpy, m)?)?;
    Ok(())
}