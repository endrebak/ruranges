use rayon::scope;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::types::PyTuple;
use pyo3::wrap_pyfunction;
use rustc_hash::FxHashMap;

use crate::boundary::sweep_line_boundary;
use crate::cluster::sweep_line_cluster;
use crate::complement::sweep_line_non_overlaps;
use crate::complement_single::sweep_line_complement;
use crate::merge::sweep_line_merge;
use crate::nearest::{sweep_line_k_nearest, pick_k_distances_combined};
use crate::overlaps;
use crate::overlaps::sweep_line_overlaps;
use crate::sorts;
use crate::spliced_subsequence::spliced_subseq;
use crate::subtract::sweep_line_subtract;

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
    slack: i64,
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
        slack,
    );
    Ok((
        result.0.into_pyarray(py).to_owned().into(),
        result.1.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
pub fn subtract_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i64>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;
    let idxs_slice2 = idxs2.as_slice()?;

    let result = sweep_line_subtract(
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
        result.2.into_pyarray(py).to_owned().into(),
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

    let mut right_result = None;
    let mut left_result: Option<(Vec<i64>, Vec<i64>, Vec<i64>)> = None;
    let mut overlap_result: Option<(Vec<i64>, Vec<i64>)> = None;

    rayon::scope(|s| {
        let right_ref = &mut right_result;
        let left_ref = &mut left_result;
        let overlap_ref = &mut overlap_result;

        s.spawn(|_| {
            let tmp = sweep_line_k_nearest(
                chrs_slice,
                ends_slice,
                idxs_slice,
                chrs_slice2,
                starts_slice2,
                idxs_slice2,
                false,
                k,
            );
            *left_ref = Some(tmp);
        });
        s.spawn(|_| {
            let tmp = sweep_line_k_nearest(
                chrs_slice,
                starts_slice,
                idxs_slice,
                chrs_slice2,
                ends_slice2,
                idxs_slice2,
                true,
                k,
            );
            *right_ref = Some(tmp);
        });
        if overlaps {
            s.spawn(|_| {
                        let tmp_overlap = sweep_line_overlaps(
                            chrs_slice,
                            starts_slice,
                            ends_slice,
                            idxs_slice,
                            chrs_slice2,
                            starts_slice2,
                            ends_slice2,
                            idxs_slice2,
                            0,
                        );
                        *overlap_ref = Some(tmp_overlap);
                    });
            }
            });
    let (r_idxs1, r_idxs2, r_dists) = right_result.unwrap();
    let (l_idxs1, l_idxs2, l_dists) = left_result.unwrap();

    let (_out_idx1, _out_idx2, _out_dists) = pick_k_distances_combined(
        &l_idxs1,
        &l_idxs2,
        &l_dists,
        &r_idxs1,
        &r_idxs2,
        &r_dists,
        k,
    );

    let (out_idx1, out_idx2, out_dists) = if overlaps {
        let (o_idxs1, o_idxs2) = overlap_result.unwrap();
        let dists = vec![0; o_idxs1.len()];
        print_vec_info("_out_idx1", &_out_idx1);
        print_vec_info("_out_idx2", &_out_idx2);
        print_vec_info("_out_dists", &_out_dists);
        print_vec_info("o_idxs1", &o_idxs1);
        print_vec_info("o_idxs2", &o_idxs2);
        print_vec_info("dists", &dists);
        pick_k_distances_combined(
            &_out_idx1,
            &_out_idx2,
            &_out_dists,
            &o_idxs1,
            &o_idxs2,
            &dists,
            k + 1,
        )
    } else {
        (_out_idx1, _out_idx2, _out_dists)
    };

    Ok((
        out_idx1.into_pyarray(py).to_owned().into(),
        out_idx2.into_pyarray(py).to_owned().into(),
        out_dists.into_pyarray(py).to_owned().into(),
    ))
}
fn print_vec_info(name: &str, data: &Vec<i64>) {
    println!("{} length: {}", name, data.len());
    println!("{} first 10: {:?}", name, &data[..data.len().min(10)]);
}


#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, idxs, slack=0))]
pub fn cluster_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    slack: i64,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (cluster_ids, indices) = sweep_line_cluster(chrs.as_slice()?, starts.as_slice()?, ends.as_slice()?, idxs.as_slice()?, slack);
    Ok((cluster_ids.into_pyarray(py).to_owned().into(), indices.into_pyarray(py).to_owned().into()))
}

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, idxs, slack=0))]
pub fn merge_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    slack: i64,
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

// #[pyfunction]
// #[pyo3(signature = (chrs, starts, ends, idxs, strand_flags, start, end = None, force_plus_strand = false))]
// pub fn subsequence_numpy(
//     chrs: PyReadonlyArray1<i64>,
//     starts: PyReadonlyArray1<i64>,
//     ends: PyReadonlyArray1<i64>,
//     idxs: PyReadonlyArray1<i64>,
//     strand_flags: PyReadonlyArray1<bool>,
//     start: i64,
//     end: Option<i64>,
//     force_plus_strand: bool,
//     py: Python,
// ) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
//     let (outidx, outstarts, outends) = crate::subsequence::subseq(chrs.as_slice()?, starts.as_slice()?, ends.as_slice()?, idxs.as_slice()?, strand_flags.as_slice()?, start, end, force_plus_strand);
//     Ok(
//         (
//             outidx.into_pyarray(py).to_owned().into(),
//             outstarts.into_pyarray(py).to_owned().into(),
//             outends.into_pyarray(py).to_owned().into(),
//         )
//     )
// }

#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(chromsweep_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(complement_overlaps_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(sort_intervals_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_intervals_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(complement_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(boundary_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(subtract_numpy, m)?)?;
//     m.add_function(wrap_pyfunction!(subsequence_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(spliced_subsequence_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(merge_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_next_intervals_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_previous_intervals_numpy, m)?)?;
    Ok(())
}


#[pyfunction]
pub fn complement_overlaps_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i64>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    idxs2: PyReadonlyArray1<i64>,
    slack: i64,
) -> PyResult<Py<PyArray1<i64>>> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;
    let idxs_slice2 = idxs2.as_slice()?;

    let result = sweep_line_non_overlaps(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        idxs_slice2,
        slack,
    );
    Ok(
        result.into_pyarray(py).to_owned().into()
    )
}

#[pyfunction]
pub fn complement_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
    slack: i64,
    chrom_len_ids: PyReadonlyArray1<i64>,
    chrom_lens: PyReadonlyArray1<i64>,
    include_first_interval: bool,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;

    let keys = chrom_len_ids.as_slice()?;
    let vals = chrom_lens.as_slice()?;

    if keys.len() != vals.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "keys array and values array must have the same length",
        ));
    }
    let mut lens_map = FxHashMap::default();
    for (&k, &v) in keys.iter().zip(vals.iter()) {
        lens_map.insert(k, v);
    }

    let (outchrs, outstarts, outends, outidxs) = sweep_line_complement(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
        slack,
        &lens_map,
        include_first_interval,
    );
    Ok(
        (
            outchrs.into_pyarray(py).to_owned().into(),
            outstarts.into_pyarray(py).to_owned().into(),
            outends.into_pyarray(py).to_owned().into(),
            outidxs.into_pyarray(py).to_owned().into()
        )
    )
}

#[pyfunction]
pub fn boundary_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    idxs: PyReadonlyArray1<i64>,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let idxs_slice = idxs.as_slice()?;

    let (outidxs, outstarts, outends, counts) = sweep_line_boundary(
        chrs_slice,
        starts_slice,
        ends_slice,
        idxs_slice,
    );
    Ok(
        (
            outidxs.into_pyarray(py).to_owned().into(),
            outstarts.into_pyarray(py).to_owned().into(),
            outends.into_pyarray(py).to_owned().into(),
            counts.into_pyarray(py).to_owned().into()
        )
    )
}