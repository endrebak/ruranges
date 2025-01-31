from time import time
import numpy as np
import pandas as pd
import polars as pl

import pyranges as pr
import ruranges


def factorize_across_columns(df, columns):
    return df.groupby(columns).ngroup()


start = time()
df = pl.read_csv(
    "/Users/endrebakkenstovner/benchmark_pyranges/downloads/generated/annotation/hg38/10_000_000/1000.bed",
    separator="\t",
    has_header=False,
).to_pandas()
df.columns = ["Chromosome", "Start", "End"]

df = pr.PyRanges(df)
df["Number"] = np.random.randint(0, 5, size=len(df))

print(df)

sf = time()
gids = factorize_across_columns(df, ["Chromosome", "Number"])
ef = time()
print("custom factorize", ef - sf)

sf = time()
# factorized, unique_vals = pd.factorize(df.Chromosome)
ef = time()
print(ef - sf)

df_len = len(df)

cluster, idx = ruranges.cluster_numpy(
    chrs=gids.to_numpy(),
    starts=df.Start.to_numpy(),
    ends=df.End.to_numpy(),
    idxs=df.index.to_numpy(),
)
end = time()

print(idx)

df_temp = df.iloc[idx]
df_temp["Cluster"] = cluster
print(df_temp)

print("time cluster", end - start)
