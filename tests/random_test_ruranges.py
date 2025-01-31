import os
import pandas as pd
import ruranges
import pyranges as pr
import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st
import pandas as pd
import subprocess
import tempfile
import os

# Suppose you have your own interval library, my_interval_lib,
# with a function `overlap_intervals(intervals_a, intervals_b)`
# that returns something akin to a list of overlapping intervals or
# a custom object. We’ll just make a placeholder import:
# from my_interval_lib import overlap_intervals

# We’ll define a strategy to generate intervals in the form [(start, end, chrom), ...]
# or [(start, end), ...]. We'll keep it simpler and ignore the chromosome/strand for now.


def interval_strategy():
    """
    Hypothesis strategy to create a list of intervals.
    Each interval is (start, end), with start <= end.
    """

    # We generate a tuple of two integers between 0 and 1000.
    # We'll ensure start <= end by sorting them after generation.
    def fix_interval(t):
        start, end = t
        start, end = (min(start, end), max(start, end))
        return start, end + 1

    return st.lists(
        st.tuples(
            st.integers(min_value=0, max_value=10000),
            st.integers(min_value=0, max_value=10000),
        ).map(fix_interval),
        min_size=1,  # at least one interval
        max_size=50,  # limit for performance
    )


@given(intervals_a=interval_strategy(), intervals_b=interval_strategy())
@settings(
    max_examples=1000,  # Increase or decrease based on desired test thoroughness
    deadline=None,  # Sometimes dealing with intervals can be slow; remove Hypothesis deadline
)
def test_overlap_algorithms(intervals_a, intervals_b):
    """
    Compare overlap results from pyranges (gold standard) with your library's result.
    """

    # 1. Convert intervals_a and intervals_b to PyRanges
    # PyRanges typically expects a DataFrame with columns: Chromosome, Start, End.
    # For simplicity, let's assume everything is on the same chromosome "chr1".
    df_a = {
        "Chromosome": [1] * len(intervals_a),
        "Start": [i[0] for i in intervals_a],
        "End": [i[1] for i in intervals_a],
    }
    df_b = {
        "Chromosome": [1] * len(intervals_b),
        "Start": [i[0] for i in intervals_b],
        "End": [i[1] for i in intervals_b],
    }

    pr_a = pr.PyRanges(df_a).sort_ranges()
    pr_b = pr.PyRanges(df_b).sort_ranges()
    print("- " * 100)
    print(pr_a)
    print(pr_b)

    # PyRanges overlap
    gold_standard_overlap = pr.PyRanges(call_with_bedtools(pr_a, pr_b)).sort_ranges()
    print("pyranges " * 10, gold_standard_overlap)

    idx1, idx2, dist = ruranges.nearest_intervals_numpy(
        chrs=pr_a.Chromosome.to_numpy().astype(np.int32),
        starts=pr_a.Start.to_numpy(),
        ends=pr_a.End.to_numpy(),
        idxs=pr_a.index.to_numpy(),
        chrs2=pr_b.Chromosome.to_numpy().astype(np.int32),
        starts2=pr_b.Start.to_numpy(),
        ends2=pr_b.End.to_numpy(),
        idxs2=pr_b.index.to_numpy(),
        k=1,
    )
    ru1 = pd.DataFrame(pr_a.loc[idx1]).reset_index(drop=True)
    ru2 = pd.DataFrame(pr_b.loc[idx2]).reset_index(drop=True)
    ru2.columns = ["C", "S", "E"]
    ru = pr.PyRanges(pd.concat([ru1, ru2, pd.DataFrame(dist)], axis=1)).sort_ranges()
    print("ruranges " * 10, ru)

    assert len(gold_standard_overlap) == len(ru)
    if len(ru) > 0:
        print(gold_standard_overlap)
        print(ru)
        assert all(gold_standard_overlap.Start.to_numpy() == ru.Start.to_numpy())
        assert all(gold_standard_overlap.End.to_numpy() == ru.End.to_numpy())
        assert all(gold_standard_overlap.Distance.to_numpy() == ru[0].to_numpy())

    # 2. Convert intervals_a and intervals_b to the format needed by your library
    # (Adjust based on how your library represents intervals.)
    # my_intervals_a = ...
    # my_intervals_b = ...
    #
    # # Then compute overlap
    # my_overlap = overlap_intervals(my_intervals_a, my_intervals_b)

    # For illustration, let's just pretend the result we get from your library is similar
    # We will convert gold_standard_overlap to a list of (start, end) so we can compare
    # directly, but normally you’d do the reverse: convert your library’s result to
    # the same shape or representation and compare to PyRanges.


def call_with_bedtools(gr1, gr2):
    # Create two temporary BED files from the input PyRanges objects
    with tempfile.NamedTemporaryFile(
        suffix=".bed", delete=False
    ) as tmp1, tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as tmp2:
        tmp1_name = tmp1.name
        tmp2_name = tmp2.name
    # Write PyRanges to these BED files
    gr1.to_bed(tmp1_name)
    gr2.to_bed(tmp2_name)
    # Construct your bedtools command. For example, bedtools closest:
    command = f"bedtools closest -t 'first' -d -a {tmp1_name} -b {tmp2_name}"
    try:
        # Run the command and capture the output
        output = subprocess.check_output(command, shell=True)
        # Decode the byte-string output to text
        output_str = output.decode("utf-8").strip()
        # Split into lines and then columns
        lines = output_str.split("\n")
        if not lines or lines == [""]:
            # If there's no result
            return pd.DataFrame(
                {"Chromsome": [], "Start": [], "End": [], "Distance": []}
            )
        rows = [
            (int(el) for i, el in enumerate(line.split("\t")) if i in {0, 1, 2, 12})
            for line in lines
        ]
        df_result = pd.DataFrame(rows)
        df_result.columns = ["Chromosome", "Start", "End", "Distance"]
        return df_result
    finally:
        # Clean up temporary files
        if os.path.exists(tmp1_name):
            os.remove(tmp1_name)
        if os.path.exists(tmp2_name):
            os.remove(tmp2_name)
