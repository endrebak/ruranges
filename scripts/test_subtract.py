import pyranges as pr
import numpy as np
import ruranges
from ruranges import chromsweep_numpy, nearest_intervals_numpy
from time import time

import pandas as pd
import polars as pl


start = time()
# df = pl.read_csv("/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/annotation/hg38/10_000_000/1000.bed", separator="\t", has_header=False).to_pandas()
# df2 = pl.read_csv("/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/reads/hg38/10_000_000/1000.bed", separator="\t", has_header=False).to_pandas()
df = pl.read_csv(
    "/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/annotation/hg38/10_000_000/1000.bed",
    separator="\t",
    has_header=False,
).to_pandas()
df2 = pl.read_csv(
    "/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/reads/hg38/10_000_000/1000.bed",
    separator="\t",
    has_header=False,
).to_pandas()
# df = pd.read_csv("/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/annotation/hg38/10_000_000/1000.bed", sep="\t")
# df2 = pd.read_csv("/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/reads/hg38/10_000_000/1000.bed", sep="\t")
end = time()
print(end - start)
print(df)

df.columns = ["Chromosome", "Start", "End"]
df2.columns = ["Chromosome", "Start", "End"]

df = pr.PyRanges(df).sort_ranges().reset_index(drop=True)
df2 = pr.PyRanges(df2).sort_ranges().reset_index(drop=True)

# df = pd.DataFrame(df.loci["chr19", 5992080:5992475])
# df2 = pd.DataFrame(df2.loci["chr19", 5000000:6000000])

print(df.head(10))
print(df2.head(10))

start = time()
combined = pd.concat([df["Chromosome"], df2["Chromosome"]], ignore_index=True)

factorized, unique_vals = pd.factorize(combined)

df_len = len(df)
df["Chromosome2"] = factorized[:df_len]
df2["Chromosome2"] = factorized[df_len:]

# Factorize it => factorized is a NumPy array of the codes; unique_vals is an Index of unique chroms
factorized, unique_vals = pd.factorize(combined)

print(df.Chromosome)

idx1, starts, ends = ruranges.subtract_numpy(
    chrs=df.Chromosome2.to_numpy(),
    starts=df.Start.to_numpy(),
    ends=df.End.to_numpy(),
    idxs=df.index.to_numpy(),
    chrs2=df2.Chromosome2.to_numpy(),
    starts2=df2.Start.to_numpy(),
    ends2=df2.End.to_numpy(),
    idxs2=df2.index.to_numpy(),
)
end = time()
print("time nearest", end - start)

temp_df1 = df.take(idx1)

final_df = pd.DataFrame(
    {
        "chr_left": temp_df1["Chromosome"].values,
        "start_left": starts,
        "end_left": ends,
    },
    index=idx1,
)

print(final_df)
print("time total nearest", time() - start)
raise

# final_df.to_csv("nearest.csv", index=False)

# res1, res2, dist = df.loc[idx1[ignore_mask]].reset_index(drop=True), df2.loc[idx2[ignore_mask]].reset_index(drop=True), pd.Series(dist[ignore_mask])
# res2.columns = ["C", "S", "E", "C2"]
# print(idx1)
# print(idx2)
# print(dist)
# res = pd.concat([res1, res2, dist], axis=1)
# print(res)
print(final_df)

print("total time nearest", time() - start)

start = time()
idx = ruranges.sort_intervals_numpy(
    chrs=res.Chromosome2.to_numpy().astype(np.int32),
    starts=res.Start.to_numpy(),
    ends=res.End.to_numpy(),
    idxs=res.index.to_numpy(),
)
print("sorted", res.loc[idx])
print("time rurange sort", time() - start)

start = time()
print(res.sort_ranges())
print("regular sort ranges", time() - start)
print(res[res[0] == 0])

# bd = pr.read_bed("result.bed")[["Chromosome", "Start", "End"]]

raise

start = time()
print(pr.PyRanges(res).sort_ranges())
print("time sort ranges", time() - start)
print(len(dist))

# start = time()
# d = ruranges.build_sorted_events(
#     chrs=df.Chromosome.to_numpy().astype(np.int32),
#     starts=df.Start.to_numpy(),
#     ends=df.End.to_numpy(),
#     idxs=df.index.to_numpy(),
#     chrs2=df2.Chromosome.to_numpy().astype(np.int32),
#     starts2=df2.Start.to_numpy(),
#     ends2=df2.End.to_numpy(),
#     idxs2=df2.index.to_numpy(),
# )
#
# end = time()
# print(end - start)
# print(d)
raise

df = pr.PyRanges(df).sort_ranges().reset_index(drop=True)
df2 = pr.PyRanges(df2).sort_ranges().reset_index(drop=True)
df["Chromosome"] = 1
df2["Chromosome"] = 1

# df = pd.read_table("/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/annotation/proteome/10_000_000/100.bed", header=None)
# df2 = pd.read_table("/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/reads/proteome/10_000_000/100.bed", header=None)

# df = pr.read_gtf("/Users/endrebakkenstovner/benchmark_pyranges/downloads/gencode/annotation.gtf")
# df2 = pr.read_bed("/Users/endrebakkenstovner/benchmark_pyranges/downloads/remc/H1_cell_line.bed")

# df2, df = df, df2

start = time()


combined = pd.concat([df["Chromosome"], df2["Chromosome"]], ignore_index=True)

# Factorize it => factorized is a NumPy array of the codes; unique_vals is an Index of unique chroms
factorized, unique_vals = pd.factorize(combined)

# Split those codes back into df and df2
df_len = len(df)
df["Chromosome"] = factorized[:df_len]
df2["Chromosome"] = factorized[df_len:]
print(df)
print(df2)

idx1, idx2, dist = ruranges.nearest_intervals_numpy(
    chrs=df.Chromosome.to_numpy().astype(np.int32),
    starts=df.Start.to_numpy(),
    ends=df.End.to_numpy(),
    idxs=df.index.to_numpy(),
    chrs2=df2.Chromosome.to_numpy().astype(np.int32),
    starts2=df2.Start.to_numpy(),
    ends2=df2.End.to_numpy(),
    idxs2=df2.index.to_numpy(),
    k=1,
)

print(dist[:10], dist[-10:])

ignore_mask = ~((idx1 == -1) | (idx2 == -1))

res1, res2, dist = (
    df.loc[idx1[ignore_mask]].reset_index(drop=True),
    df2.loc[idx2[ignore_mask]].reset_index(drop=True),
    pd.Series(dist[ignore_mask]),
)
print(res1)
print(res2)
print(dist)
res = pd.concat([res1, res2, dist], axis=1)

print(time() - start)

print(res)
print(len(dist))

start = time()
res2 = df.nearest(df2, exclude_overlaps=True).reset_index(drop=True)
print(time() - start)
print(res2)

difference_mask = res[0] != res2.Distance
print(res[difference_mask])
print(pd.DataFrame(res2[difference_mask]))
