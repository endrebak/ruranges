use std::str::FromStr;
use std::time::Instant;

use numpy::PyArrayMethods;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::types::PyTuple;
use pyo3::wrap_pyfunction;
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;

use crate::boundary::sweep_line_boundary;
use crate::cluster::sweep_line_cluster;
use crate::complement::sweep_line_non_overlaps;
use crate::complement_single::sweep_line_complement;
use crate::merge::sweep_line_merge;
use crate::nearest::nearest;
// use crate::nearest::nearest;
use crate::overlaps::{self, sweep_line_overlaps_containment};
use crate::ruranges_structs::OverlapPair;
use crate::sorts;
use crate::sorts::build_sorted_events_single_collection_separate_outputs;
use crate::spliced_subsequence::spliced_subseq;
use crate::split::sweep_line_split;
use crate::subtract::sweep_line_subtract;
use crate::tile::{tile, window};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum OverlapType {
    First,
    Last,
    All,
}

impl FromStr for OverlapType {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "all" => Ok(OverlapType::All),
            "first" => Ok(OverlapType::First),
            "last" => Ok(OverlapType::Last),
            _ => Err("Invalid direction string"),
        }
    }
}

#[pyfunction]
pub fn chromsweep_numpy(
    py: Python,
    chrs: PyReadonlyArray1<u32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<u32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    slack: i64,
    overlap_type: &str,
    contained: bool,
) -> PyResult<(Py<PyArray1<u32>>, Py<PyArray1<u32>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;

    let overlap_type = OverlapType::from_str(overlap_type).unwrap();
    let invert = overlap_type == OverlapType::Last;

    let result = if overlap_type == OverlapType::All && !contained {
        // The common, super-optimized case
        overlaps::sweep_line_overlaps(
            chrs_slice,
            starts_slice,
            ends_slice,
            chrs_slice2,
            starts_slice2,
            ends_slice2,
            slack,
        )
    } else {
        if !contained {
            let (sorted_starts, sorted_ends) = overlaps::compute_sorted_events(
                chrs_slice,
                starts_slice,
                ends_slice,
                slack,
                invert,
            );
            let (sorted_starts2, sorted_ends2) =
                overlaps::compute_sorted_events(chrs_slice2, starts_slice2, ends_slice2, 0, invert);

            let mut pairs = overlaps::sweep_line_overlaps_overlap_pair(
                &sorted_starts,
                &sorted_ends,
                &sorted_starts2,
                &sorted_ends2,
            );
            keep_first_by_idx(&mut pairs);
            pairs.into_iter().map(|pair| (pair.idx, pair.idx2)).unzip()
        } else {
            let maxevents = overlaps::compute_sorted_maxevents(
                chrs_slice,
                starts_slice,
                ends_slice,
                chrs_slice2,
                starts_slice2,
                ends_slice2,
                slack,
                invert,
            );
            let mut pairs = overlaps::sweep_line_overlaps_containment(maxevents);
            if overlap_type == OverlapType::All {
                pairs.into_iter().map(|pair| (pair.idx, pair.idx2)).unzip()
            } else {
                keep_first_by_idx(&mut pairs);
                pairs.into_iter().map(|pair| (pair.idx, pair.idx2)).unzip()
            }
        }
    };

    let res = Ok((
        result.0.into_pyarray(py).to_owned().into(),
        result.1.into_pyarray(py).to_owned().into(),
    ));
    res
}

fn keep_first_by_idx(pairs: &mut Vec<OverlapPair>) {
    let mut seen_idx = FxHashSet::default();
    pairs.retain(|pair| seen_idx.insert(pair.idx));
}

#[pyfunction]
#[pyo3(signature = (*, chrs, starts, ends, chrs2, starts2, ends2, slack=0, k=1, include_overlaps=true, direction="any"))]
pub fn nearest_numpy(
    py: Python,
    chrs: PyReadonlyArray1<u32>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<u32>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    slack: i64,
    k: usize,
    include_overlaps: bool,
    direction: &str,
) -> PyResult<(Py<PyArray1<u32>>, Py<PyArray1<u32>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;

    let result = nearest(
        chrs_slice,
        starts_slice,
        ends_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        slack,
        k,
        include_overlaps,
        direction,
    );
    let res = Ok((
        result.0.into_pyarray(py).to_owned().into(),
        result.1.into_pyarray(py).to_owned().into(),
        result.2.into_pyarray(py).to_owned().into(),
    ));
    res
}

#[pyfunction]
pub fn subtract_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i64>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
) -> PyResult<(Py<PyArray1<usize>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;

    let result = sweep_line_subtract(
        chrs_slice,
        starts_slice,
        ends_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
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
    py: Python,
) -> PyResult<Py<PyArray1<usize>>> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;

    let indexes = sorts::sort_order_idx(chrs_slice, starts_slice, ends_slice);
    Ok(indexes.into_pyarray(py).to_owned().into())
}

// #[pyfunction]
// pub fn nearest_intervals_unique_k_numpy(
//     chrs: PyReadonlyArray1<i64>,
//     starts: PyReadonlyArray1<i64>,
//     ends: PyReadonlyArray1<i64>,
//     idxs: PyReadonlyArray1<i64>,
//     chrs2: PyReadonlyArray1<i64>,
//     starts2: PyReadonlyArray1<i64>,
//     ends2: PyReadonlyArray1<i64>,
//     idxs2: PyReadonlyArray1<i64>,
//     k: usize,
//     overlaps: bool,
//     direction: &str,
//     py: Python,
// ) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
//     let dir: Direction = direction.parse().expect("Invalid direction, must be forward, backwards, or any.");
//
//     let chrs_slice = chrs.as_slice()?;
//     let starts_slice = starts.as_slice()?;
//     let ends_slice = ends.as_slice()?;
//     let idxs_slice = idxs.as_slice()?;
//     let chrs_slice2 = chrs2.as_slice()?;
//     let starts_slice2 = starts2.as_slice()?;
//     let ends_slice2 = ends2.as_slice()?;
//     let idxs_slice2 = idxs2.as_slice()?;
//
//     let mut right_result = Some((Vec::new(), Vec::new(), Vec::new()));
//     let mut left_result = Some((Vec::new(), Vec::new(), Vec::new()));
//     let mut overlap_result: Option<(Vec<i64>, Vec<i64>)> = None;
//
//     rayon::scope(|s| {
//         let right_ref = &mut right_result;
//         let left_ref = &mut left_result;
//         let overlap_ref = &mut overlap_result;
//         if dir != Direction::Forward {
//             s.spawn(|_| {
//                 let tmp = sweep_line_k_nearest(
//                     chrs_slice,
//                     ends_slice,
//                     idxs_slice,
//                     chrs_slice2,
//                     starts_slice2,
//                     idxs_slice2,
//                     false,
//                     k,
//                 );
//                 *left_ref = Some(tmp);
//             });
//         }
//         if dir != Direction::Backward {
//             s.spawn(|_| {
//                 let tmp = sweep_line_k_nearest(
//                     chrs_slice,
//                     starts_slice,
//                     idxs_slice,
//                     chrs_slice2,
//                     ends_slice2,
//                     idxs_slice2,
//                     true,
//                     k,
//                 );
//                 *right_ref = Some(tmp);
//             });
//         }
//         if overlaps {
//             s.spawn(|_| {
//                         let tmp_overlap = sweep_line_overlaps(
//                             chrs_slice,
//                             starts_slice,
//                             ends_slice,
//                             idxs_slice,
//                             chrs_slice2,
//                             starts_slice2,
//                             ends_slice2,
//                             idxs_slice2,
//                             0,
//                         );
//                         *overlap_ref = Some(tmp_overlap);
//                     });
//             }
//             });
//     let (r_idxs1, r_idxs2, r_dists) = right_result.unwrap();
//     let (l_idxs1, l_idxs2, l_dists) = left_result.unwrap();
//
//     let (_out_idx1, _out_idx2, _out_dists) = pick_k_distances_combined(
//         &l_idxs1,
//         &l_idxs2,
//         &l_dists,
//         &r_idxs1,
//         &r_idxs2,
//         &r_dists,
//         k,
//     );
//
//     let (out_idx1, out_idx2, out_dists) = if overlaps {
//         let (o_idxs1, o_idxs2) = overlap_result.unwrap();
//         let dists = vec![0; o_idxs1.len()];
//         pick_k_distances_combined(
//             &_out_idx1,
//             &_out_idx2,
//             &_out_dists,
//             &o_idxs1,
//             &o_idxs2,
//             &dists,
//             k,
//         )
//     } else {
//         (_out_idx1, _out_idx2, _out_dists)
//     };
//
//     Ok((
//         out_idx1.into_pyarray(py).to_owned().into(),
//         out_idx2.into_pyarray(py).to_owned().into(),
//         out_dists.into_pyarray(py).to_owned().into(),
//     ))
// }

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, slack=0))]
pub fn cluster_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    slack: i64,
    py: Python,
) -> PyResult<(Py<PyArray1<i64>>, Py<PyArray1<usize>>)> {
    let (cluster_ids, indices) = sweep_line_cluster(
        chrs.as_slice()?,
        starts.as_slice()?,
        ends.as_slice()?,
        slack,
    );
    Ok((
        cluster_ids.into_pyarray(py).to_owned().into(),
        indices.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
#[pyo3(signature = (starts, ends, negative_strand, tile_size))]
pub fn tile_numpy(
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    negative_strand: PyReadonlyArray1<bool>,
    tile_size: i64,
    py: Python,
) -> PyResult<(
    Py<PyArray1<usize>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<f64>>,
)> {
    let (starts, ends, indices, overlap_fraction) =
        tile(starts.as_slice()?, ends.as_slice()?, negative_strand.as_slice()?, tile_size);
    Ok((
        indices.into_pyarray(py).to_owned().into(),
        starts.into_pyarray(py).to_owned().into(),
        ends.into_pyarray(py).to_owned().into(),
        overlap_fraction.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
#[pyo3(signature = (starts, ends, negative_strand, window_size))]
pub fn window_numpy(
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    negative_strand: PyReadonlyArray1<bool>,
    window_size: i64,
    py: Python,
) -> PyResult<(Py<PyArray1<usize>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (starts, ends, indices) = window(starts.as_slice()?, ends.as_slice()?, negative_strand.as_slice()?, window_size);
    Ok((
        indices.into_pyarray(py).to_owned().into(),
        starts.into_pyarray(py).to_owned().into(),
        ends.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, slack=0))]
pub fn merge_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    slack: i64,
    py: Python,
) -> PyResult<(
    Py<PyArray1<usize>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
)> {
    let (indices, starts, ends, counts) = sweep_line_merge(
        chrs.as_slice()?,
        starts.as_slice()?,
        ends.as_slice()?,
        slack,
    );
    Ok((
        indices.into_pyarray(py).to_owned().into(),
        starts.into_pyarray(py).to_owned().into(),
        ends.into_pyarray(py).to_owned().into(),
        counts.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, slack=0, between=false))]
pub fn split_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    slack: i64,
    between: bool,
    py: Python,
) -> PyResult<(Py<PyArray1<usize>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (indices, starts, ends) = sweep_line_split(
        chrs.as_slice()?,
        starts.as_slice()?,
        ends.as_slice()?,
        slack,
        between,
    );
    Ok((
        indices.into_pyarray(py).to_owned().into(),
        starts.into_pyarray(py).to_owned().into(),
        ends.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
#[pyo3(signature = (chrs, starts, ends, strand_flags, start, end = None, force_plus_strand = false))]
pub fn spliced_subsequence_numpy(
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    strand_flags: PyReadonlyArray1<bool>,
    start: i64,
    end: Option<i64>,
    force_plus_strand: bool,
    py: Python,
) -> PyResult<(Py<PyArray1<usize>>, Py<PyArray1<i64>>, Py<PyArray1<i64>>)> {
    let (outidx, outstarts, outends) = spliced_subseq(
        chrs.as_slice()?,
        starts.as_slice()?,
        ends.as_slice()?,
        strand_flags.as_slice()?,
        start,
        end,
        force_plus_strand,
    );
    Ok((
        outidx.into_pyarray(py).to_owned().into(),
        outstarts.into_pyarray(py).to_owned().into(),
        outends.into_pyarray(py).to_owned().into(),
    ))
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

#[pyfunction]
pub fn complement_overlaps_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    chrs2: PyReadonlyArray1<i64>,
    starts2: PyReadonlyArray1<i64>,
    ends2: PyReadonlyArray1<i64>,
    slack: i64,
) -> PyResult<Py<PyArray1<usize>>> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;
    let chrs_slice2 = chrs2.as_slice()?;
    let starts_slice2 = starts2.as_slice()?;
    let ends_slice2 = ends2.as_slice()?;

    let result = sweep_line_non_overlaps(
        chrs_slice,
        starts_slice,
        ends_slice,
        chrs_slice2,
        starts_slice2,
        ends_slice2,
        slack,
    );
    Ok(result.into_pyarray(py).to_owned().into())
}

#[pyfunction]
pub fn complement_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
    slack: i64,
    chrom_len_ids: PyReadonlyArray1<i64>,
    chrom_lens: PyReadonlyArray1<i64>,
    include_first_interval: bool,
) -> PyResult<(
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<usize>>,
)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;

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
        slack,
        &lens_map,
        include_first_interval,
    );
    Ok((
        outchrs.into_pyarray(py).to_owned().into(),
        outstarts.into_pyarray(py).to_owned().into(),
        outends.into_pyarray(py).to_owned().into(),
        outidxs.into_pyarray(py).to_owned().into(),
    ))
}

#[pyfunction]
pub fn boundary_numpy(
    py: Python,
    chrs: PyReadonlyArray1<i64>,
    starts: PyReadonlyArray1<i64>,
    ends: PyReadonlyArray1<i64>,
) -> PyResult<(
    Py<PyArray1<usize>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
    Py<PyArray1<i64>>,
)> {
    let chrs_slice = chrs.as_slice()?;
    let starts_slice = starts.as_slice()?;
    let ends_slice = ends.as_slice()?;

    let (outidxs, outstarts, outends, counts) =
        sweep_line_boundary(chrs_slice, starts_slice, ends_slice);
    Ok((
        outidxs.into_pyarray(py).to_owned().into(),
        outstarts.into_pyarray(py).to_owned().into(),
        outends.into_pyarray(py).to_owned().into(),
        counts.into_pyarray(py).to_owned().into(),
    ))
}


#[derive(Debug, PartialEq)]
enum Direction {
    Forward,
    Backward,
    Any,
}

impl FromStr for Direction {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "forward" => Ok(Direction::Forward),
            "backward" => Ok(Direction::Backward),
            "any" => Ok(Direction::Any),
            _ => Err(format!("Invalid direction: {}", s)),
        }
    }
}

#[pymodule]
fn ruranges(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(chromsweep_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(complement_overlaps_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(window_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(tile_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(sort_intervals_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_intervals_unique_k_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(complement_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(boundary_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(subtract_numpy, m)?)?;
    //     m.add_function(wrap_pyfunction!(subsequence_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(spliced_subsequence_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(merge_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(split_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_next_intervals_numpy, m)?)?;
    // m.add_function(wrap_pyfunction!(nearest_previous_intervals_numpy, m)?)?;
    Ok(())
}
