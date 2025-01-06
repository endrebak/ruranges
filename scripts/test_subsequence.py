from time import time
import numpy as np
import pandas as pd
import polars as pl

import pyranges as pr
import ruranges

def factorize_across_columns(df, columns):
    return df.groupby(columns).ngroup()

p  = pr.PyRanges({"Chromosome": [1, 1, 2, 2, 3],
                  "Strand": ["+", "+", "-", "-", "+"],
                  "Start": [1, 40, 10, 70, 140],
                  "End": [11, 60, 25, 80, 152],
                  "transcript_id":["t1", "t1", "t2", "t2", "t3"] })

print(p)

def ss(df, col, start, end):

    gids = factorize_across_columns(df, ["Chromosome", col])
    strands = df.Strand.replace({"+": True, "-": False}).astype(bool).values
    print(strands)
    idx, starts, ends = ruranges.spliced_subsequence_numpy(
        chrs=gids.to_numpy(),
        starts=df.Start.to_numpy(),
        ends=df.End.to_numpy(),
        idxs=df.index.to_numpy(),
        strand_flags=strands,
        start=start,
        end=end,
    )
    
    df_temp = df.reindex(idx)
    df_temp.loc[:, "Start"] = starts
    df_temp.loc[:, "End"] = ends
    return df_temp

res = ss(p, "transcript_id", 3, -3)
print(res)

start = time()
df = pl.read_csv("/Users/endrebakkenstovner/benchmark_pyranges/downloads/gencode/annotation.bed", separator="\t", has_header=True, infer_schema_length=10000).to_pandas()
df = df[["Chromosome", "Start", "End", "Feature", "Strand", "gene_id", "transcript_id"]]

df = df[df.Feature == "exon"]

df = pr.PyRanges(df)

print(df)

sf = time()
gids = factorize_across_columns(df, ["Chromosome", "Strand", "gene_id"])
ef = time()
print("custom factorize", ef - sf)

sf = time()
# factorized, unique_vals = pd.factorize(df.Chromosome)
ef = time()
print(ef - sf)

df_len = len(df)

idx, starts, ends = ruranges.spliced_subsequence_numpy(
    chrs=gids.to_numpy(),
    starts=df.Start.to_numpy(),
    ends=df.End.to_numpy(),
    idxs=df.index.to_numpy(),
    strand_flags=(df.Strand == "+").values,
    start=0,
    end=100,
)
end = time()

print(idx)

df_temp = df.reindex(idx)
df_temp.loc[:, "Start"] = starts
df_temp.loc[:, "End"] = ends
print(df_temp)
print(df_temp[df_temp.transcript_id.duplicated(keep=False)])

print("time spliced_subsequence", end - start)