import pandas as pd
import ruranges
import pyranges as pr
import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

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
            st.integers(min_value=0, max_value=1000),
            st.integers(min_value=0, max_value=1000)
        ).map(fix_interval),
        min_size=1,  # at least one interval
        max_size=50  # limit for performance
    )

@given(
    intervals_a=interval_strategy(),
    intervals_b=interval_strategy()
)
@settings(
    max_examples=1000,  # Increase or decrease based on desired test thoroughness
    deadline=None       # Sometimes dealing with intervals can be slow; remove Hypothesis deadline
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

    pr_a = pr.PyRanges(df_a)
    pr_b = pr.PyRanges(df_b)
    print("- " * 100)
    print(pr_a)
    print(pr_b)

    # PyRanges overlap
    gold_standard_overlap = pr_a.overlap(pr_b).sort_ranges()
    print("pyranges " * 10, gold_standard_overlap)

    idx1, idx2 = ruranges.chromsweep(
        pr_a.Chromosome,
        pr_a.Start,
        pr_a.End,
        pr_a.index,
        pr_b.Chromosome,
        pr_b.Start,
        pr_b.End,
        pr_b.index,
    )
    ru = pr_a.loc[idx1].sort_ranges()
    print("ruranges " * 10, ru)

    
    assert gold_standard_overlap.Start == ru.Start
    assert gold_standard_overlap.End == ru.End

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

