use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use numpy::{PyReadonlyArray1, PyArray1, IntoPyArray};

use crate::nearest;
use crate::nearest::nearest_k;
use crate::nearest::nearest_nonoverlapping_k;
use crate::sorts;
use crate::overlaps;
use crate::sorts::build_intervals;

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

    let result = overlaps::sweep_line_overlaps(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        idxs_slice2,
    );
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

    let indexes = sorts::sort_order_idx(chrs_slice, starts_slice, ends_slice, idxs_slice);
    Ok(indexes.into_pyarray(py).to_owned().into())
}

#[pyfunction]
pub fn nearest_next_intervals_numpy(
    chrs: PyReadonlyArray1<i32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice   = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice   = ends.as_slice()?;
    let idxs_slice   = idxs.as_slice()?;
    let chrs_slice2   = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2   = ends2.as_slice()?;
    let idxs_slice2   = idxs2.as_slice()?;

    let mut is1 = build_intervals(chrs_slice, starts_slice, ends_slice, idxs_slice);
    let mut is2 = build_intervals(chrs_slice2, starts_slice2, ends_slice2, idxs_slice2);

    let indexes_and_dists = nearest::nearest_next_nonoverlapping_k(is1.as_mut_slice(), is2.as_mut_slice(), 2);
    Ok(
        (
        indexes_and_dists.0.into_pyarray(py).to_owned().into(),
        indexes_and_dists.1.into_pyarray(py).to_owned().into(),
        indexes_and_dists.2.into_pyarray(py).to_owned().into(),
        )
    )
}

#[pyfunction]
pub fn nearest_previous_intervals_numpy(
    chrs: PyReadonlyArray1<i32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice   = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice   = ends.as_slice()?;
    let idxs_slice   = idxs.as_slice()?;
    let chrs_slice2   = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2   = ends2.as_slice()?;
    let idxs_slice2   = idxs2.as_slice()?;

    let mut is1 = build_intervals(chrs_slice, starts_slice, ends_slice, idxs_slice);
    let mut is2 = build_intervals(chrs_slice2, starts_slice2, ends_slice2, idxs_slice2);

    let indexes_and_dists = nearest::nearest_previous_nonoverlapping_k(is1.as_mut_slice(), is2.as_mut_slice(), 2);
    Ok(
        (
        indexes_and_dists.0.into_pyarray(py).to_owned().into(),
        indexes_and_dists.1.into_pyarray(py).to_owned().into(),
        indexes_and_dists.2.into_pyarray(py).to_owned().into(),
        )
    )
}

#[pyfunction]
pub fn nearest_nonoverlapping_intervals_numpy(
    chrs: PyReadonlyArray1<i32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
    k: usize,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice   = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice   = ends.as_slice()?;
    let idxs_slice   = idxs.as_slice()?;
    let chrs_slice2   = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2   = ends2.as_slice()?;
    let idxs_slice2   = idxs2.as_slice()?;

    let indexes_and_dists = nearest_nonoverlapping_k(chrs_slice, starts_slice, ends_slice, idxs_slice, chrs_slice2, starts_slice2, ends_slice2, idxs_slice2, k);
    Ok(
        (
        indexes_and_dists.0.into_pyarray(py).to_owned().into(),
        indexes_and_dists.1.into_pyarray(py).to_owned().into(),
        indexes_and_dists.2.into_pyarray(py).to_owned().into(),
        )
    )
}

#[pyfunction]
pub fn nearest_intervals_numpy(
    chrs: PyReadonlyArray1<i32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
    k: usize,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice   = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice   = ends.as_slice()?;
    let idxs_slice   = idxs.as_slice()?;
    let chrs_slice2   = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2   = ends2.as_slice()?;
    let idxs_slice2   = idxs2.as_slice()?;

    let (idx1, idx2, d) = nearest_k(chrs_slice, starts_slice, ends_slice, idxs_slice, chrs_slice2, starts_slice2, ends_slice2, idxs_slice2, k);
    Ok(
        (
            idx1.into_pyarray(py).to_owned().into(),
            idx2.into_pyarray(py).to_owned().into(),
            d.into_pyarray(py).to_owned().into(),
        )
    )
}


#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(chromsweep_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(sort_intervals_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_next_intervals_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_previous_intervals_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_intervals_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_nonoverlapping_intervals_numpy, m)?)?;
    Ok(())
}