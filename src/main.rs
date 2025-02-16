use std::hash::Hash;
use clap::Parser;
use num_traits::{PrimInt, Signed, Zero};
use polars::prelude::*;
use polars::datatypes::DataType;


use ruranges::overlaps;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::time::Instant;

/// Simple program to read two CSV files as DataFrames using Polars.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the first CSV file
    input1: PathBuf,

    /// Path to the second CSV file
    input2: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments
    let args = Args::parse();

    let columns = Arc::new(vec![0, 1, 2]);

    // Open the first CSV file and create a CsvReader
    let parse_options: CsvParseOptions = CsvParseOptions::default().with_separator(b'\t');
    let fields = vec![
        Field::new("column_1".into(), DataType::Categorical(None, CategoricalOrdering::Physical)),
        Field::new("column_2".into(), DataType::Int32),
        Field::new("column_3".into(), DataType::Int32),
        // Add more fields as needed.
    ];
    let schema = Schema::from_iter(fields);

    
    // Build a schema from the fields.
    let file1 = args.input1;
    let csv_reader = CsvReadOptions::default()
    .with_has_header(false)
    .with_schema_overwrite(Some(std::sync::Arc::new(schema.clone())))
    .with_projection(Some(columns.clone()))  // Add this line
    .with_rechunk(true)
    .with_parse_options(parse_options);
    let csv = csv_reader.clone().try_into_reader_with_file_path(Some(file1.clone()))?
    .finish()?;

    let parse_options: CsvParseOptions = CsvParseOptions::default().with_separator(b'\t');
    let file2 = args.input2;
    let csv2 = CsvReadOptions::default()
    .with_has_header(false)
    .with_schema_overwrite(Some(std::sync::Arc::new(schema)))
    .with_projection(Some(columns))  // Add this line
    .with_rechunk(true)
    .with_parse_options(parse_options)
    .try_into_reader_with_file_path(Some(file2))?
    .finish()?;

    // Print the DataFrames to the terminal
    println!("DataFrame 1:");
    println!("{:?}", csv);

    println!("\nDataFrame 2:");
    println!("{:?}", csv2);

    // Force contiguous memory layout


    let start = Instant::now();
    let col1 = csv.column("column_1")?.as_series().ok_or("Failed to convert column_1 to Series")?;
    let col2 = csv2.column("column_1")?.as_series().ok_or("Failed to convert column_1 to Series")?;
    let (c1, c2) = factorize_binary(col1, col2)?;
    let (c1, c2) = process_columns::<UInt32Type>(&csv, &csv2, "column_1".into())?;
    println!("{:?}", start.elapsed());
    let starts = csv.column("column_2")?.i32()?;
    println!("{:?}", start.elapsed());
    let ends = csv.column("column_3")?.i32()?;
    println!("{:?}", start.elapsed());
    println!("{:?}", start.elapsed());
    let starts2 = csv2.column("column_2")?.i32()?;
    println!("{:?}", start.elapsed());
    let ends2 = csv2.column("column_3")?.i32()?;
    println!("{:?}", start.elapsed());
    let (mut idx, _) = overlaps::sweep_line_overlaps(
        c1.cont_slice()?,
        starts.cont_slice()?,
        ends.cont_slice()?,
        c2.cont_slice()?,
        starts2.cont_slice()?,
        ends2.cont_slice()?,
        0_i32,
    );
    println!("{:?}", idx.len());
    radsort::sort(&mut idx);
    println!("{:?}", idx.len());

    // let idx_ca = UInt32Chunked::from_vec("idx".into(), idx);
    // let mut res = csv.take(&idx_ca)?;
    // let stdout = std::io::stdout();
    // let mut handle = stdout.lock();
    // CsvWriter::new(&mut handle)
    //     .include_header(false)
    //     .with_separator(b'\t')
    //     .finish(&mut res)?;
    // let lf = LazyCsvReader::new(file1)
    // .with_has_header(false)
    // .finish()?;
    // process_chunks(lf, idx, 1_000_000u32)?;
    
    Ok(())
}

fn process_columns<T>(
    csv1: &DataFrame,
    csv2: &DataFrame,
    column: &str,
) -> PolarsResult<(ChunkedArray<T>, ChunkedArray<T>)>
where
    T: PolarsNumericType,
    // The *native* integer type must be prim-int-like etc.
    T::Native: PrimInt + Signed + Hash + Copy + radsort::Key + Zero + TryFrom<usize>, {
    // Create a local map that will be dropped when this function ends.
    let mut global_map: FxHashMap<String, T::Native> = FxHashMap::default();
    
    let c1 = encode_strings_to_codes(csv1.column(column)?, &mut global_map)?;
    let c2 = encode_strings_to_codes(csv2.column(column)?, &mut global_map)?;
    
    Ok((c1, c2))
}

fn encode_strings_to_codes<T>(
    s: &Column,
    global_map: &mut FxHashMap<String, T::Native>,
) -> PolarsResult<ChunkedArray<T>> 
where
    T: PolarsNumericType,
    // Require T::Native to be convertible from usize.
    T::Native: PrimInt + Signed + Hash + Copy + radsort::Key + Zero + TryFrom<usize>,
{
    if s.dtype() != &DataType::String {
        return Err(PolarsError::ComputeError(
            format!("Series '{}' is not String", s.name()).into(),
        ));
    }
    
    let str_chunked = s.str()?;
    
    if global_map.capacity() - global_map.len() < str_chunked.len() {
        global_map.reserve(str_chunked.len());
    }

    // Map each string (or null) to a code.
    let codes: Vec<Option<T::Native>> = str_chunked
        .into_iter()
        .map(|opt_str| -> PolarsResult<Option<T::Native>> {
            if let Some(val) = opt_str {
                let len_as_native = T::Native::try_from(global_map.len())
                    .map_err(|_| {
                        PolarsError::ComputeError("Conversion from usize to T::Native failed".into())
                    })?;
                Ok(Some(*global_map.entry(val.to_string()).or_insert(len_as_native)))
            } else {
                Ok(None)
            }
        })
        .collect::<PolarsResult<Vec<Option<T::Native>>>>()?;
    
    // Create the ChunkedArray from the iterator over Option<T::Native>.
    Ok(ChunkedArray::<T>::from_iter_options(s.name().clone(), codes.into_iter()))
}

fn factorize_binary(
    series1: &Series,
    series2: &Series,
) -> PolarsResult<(UInt32Chunked, UInt32Chunked)> {
    // First create a temporary categorical Series that contains all unique values
    let mut tmp = vec![];
    tmp.extend_from_slice(series1.iter().collect::<Vec<_>>().as_slice());
    tmp.extend_from_slice(series2.iter().collect::<Vec<_>>().as_slice());
    let tmp_series = Series::new("tmp".into(), tmp);
    
    // Cast to categorical to get a unified mapping
    let cat_type = DataType::Categorical(None, CategoricalOrdering::Physical);
    let unified_cat = tmp_series.cast(&cat_type)?;
    
    // Now cast both original series using this unified mapping
    let cat1 = series1.cast(&unified_cat.dtype())?;
    let cat2 = series2.cast(&unified_cat.dtype())?;
    
    // Extract codes
    let codes1 = cat1.to_physical_repr().u32()?.clone();
    let codes2 = cat2.to_physical_repr().u32()?.clone();
    
    Ok((codes1, codes2))
}


fn process_chunks(
    lf: LazyFrame,
    sorted_indices: Vec<u32>,
    chunk_size: u32,
) -> PolarsResult<()> {
    use polars::prelude::*;
    let stdout = std::io::stdout();
    let mut handle = stdout.lock();
    let total_rows = lf.clone().collect()?.height() as u32;
    let mut indices_iter = sorted_indices.into_iter().peekable();
    let mut chunk_start = 0u32;
    
    while chunk_start < total_rows && indices_iter.peek().is_some() {
        let chunk_end = (chunk_start + chunk_size).min(total_rows);
        let mut chunk_indices = Vec::new();
        
        while let Some(&idx) = indices_iter.peek() {
            if idx >= chunk_end {
                break;
            }
            chunk_indices.push(indices_iter.next().unwrap());
        }
        
        if chunk_indices.is_empty() {
            chunk_start = chunk_end;
            continue;
        }
        
        let partial_df = lf
            .clone()
            .slice(chunk_start as i64, (chunk_end - chunk_start) as u32)
            .collect()?;
            
        // Create a Series of local indices that preserves duplicates
        let local_offsets: Vec<u32> = chunk_indices
            .iter()
            .map(|&idx| idx - chunk_start)
            .collect();
            
        // Create a new DataFrame with just the rows we need, potentially duplicated
        let idx_ca = UInt32Chunked::from_vec("idx".into(), local_offsets);

        let mut subset = partial_df.take(&idx_ca)?;
        
        // Write all rows at once
        CsvWriter::new(&mut handle)
            .include_header(false)
            .with_separator(b'\t')
            .finish(&mut subset)?;
            
        chunk_start = chunk_end;
    }
    Ok(())
}