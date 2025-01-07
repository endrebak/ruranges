use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::cluster::sweep_line_cluster;
use crate::merge::sweep_line_merge;
use crate::nearest;
use crate::nearest::sweep_line_k_nearest;
use crate::overlaps;
use crate::sorts;
use crate::spliced_subsequence::spliced_subseq;

#[pyfunction]
pub fn chromsweep_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i64>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;
    let idxs_slice2 = idxs2.as_slice()?;

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
    Ok((
        result.0.into_pyarray(py).to_owned().into(),
        result.1.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
pub fn sort_intervals_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    py: Python,
) -> PyResult<Py<PyArray1<i64>>> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;

    let indexes = sorts::sort_order_idx(chrs_slice, starts_slice, ends_slice, idxs_slice);
    Ok(indexes.into_pyarray(py).to_owned().into())
}

#[pyfunction]
pub fn nearest_intervals_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i64>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
    k: usize,
    overlaps: bool,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;
    let idxs_slice2 = idxs2.as_slice()?;

    let (idx1, idx2, d) = sweep_line_k_nearest(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        idxs_slice2,
        k,
        overlaps,
    );
    Ok((
        idx1.into_pyarray(py).to_owned().into(),
        idx2.into_pyarray(py).to_owned().into(),
        d.into_pyarray(py).to_owned().into(),
    ))
}


#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, idxs, slack=None))]
pub fn cluster_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    slack: Option<i64>,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (cluster_ids, indices) = sweep_line_cluster(chrs.as_slice()?, starts.as_slice()?, ends.as_slice()?, idxs.as_slice()?, slack);
    Ok((cluster_ids.into_pyarray(py).to_owned().into(), indices.into_pyarray(py).to_owned().into()))
}

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, idxs, slack=None))]
pub fn merge_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    slack: Option<i64>,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (indices, starts, ends, counts) = sweep_line_merge(chrs.as_slice()?, starts.as_slice()?, ends.as_slice()?, idxs.as_slice()?, slack);
    Ok((indices.into_pyarray(py).to_owned().into(), starts.into_pyarray(py).to_owned().into(), ends.into_pyarray(py).to_owned().into(), counts.into_pyarray(py).to_owned().into()))
}

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, idxs, strand_flags, start, end = None, force_plus_strand = false))]
pub fn spliced_subsequence_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    strand_flags: PyReadonlyArray1<bool>,
    start: i64,
    end: Option<i64>,
    force_plus_strand: bool,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (outidx, outstarts, outends) = spliced_subseq(chrs.as_slice()?, starts.as_slice()?, ends.as_slice()?, idxs.as_slice()?, strand_flags.as_slice()?, start, end, force_plus_strand);
    Ok(
        (
            outidx.into_pyarray(py).to_owned().into(),
            outstarts.into_pyarray(py).to_owned().into(),
            outends.into_pyarray(py).to_owned().into(),
        )
    )
}

#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(chromsweep_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(sort_intervals_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_intervals_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(spliced_subsequence_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(merge_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_next_intervals_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_previous_intervals_numpy, m)?)?;
    Ok(())
}
