import pandas as pd
import numpy as np

import ruranges


def test_nearest_1():
    intervals_a = [(0, 1), (1, 2)]
    intervals_b = [(0, 1), (2, 3)]

    df, df2 = create_dataframes(intervals_a, intervals_b=intervals_b)

    ru = call_nearest(df, df2)
    print(ru)

    assert [*ru.Distance] == [1]
    assert [*ru.Start] == [1]
    assert [*ru.End] == [2]
    assert [*ru.S] == [0]
    assert [*ru.E] == [1]


# def test_nearest2():
#     df, df2 = create_dataframes(intervals_a=[(0, 3), (1, 2)], intervals_b=[(3, 4)])
#     print(df)
#     print(df2)
#     ru = call_nearest(df, df2)
#     print(ru)
#     assert [*ru.Distance] == [1, 2]
#     assert [*ru.Start] == [0, 1]
#     assert [*ru.End] == [3, 2]
#     assert [*ru.S] == [3, 3]
#     assert [*ru.E] == [4, 4]
#
# def test_nearest3():
#     df, df2 = create_dataframes(intervals_a=[(0, 3), (1, 2)], intervals_b=[(0, 1), (2, 3)])
#     print(df)
#     print(df2)
#     ru = call_nearest(df, df2)
#     print(ru)
#     assert [*ru.Distance] == [1]
#     assert [*ru.Start] == [1]
#     assert [*ru.End] == [2]
#     assert [*ru.S] == [0]
#     assert [*ru.E] == [1]
#
# def test_nearest4():
#     df, df2 = create_dataframes(
#         intervals_a=[(1, 33), (3, 5)],
#         intervals_b=[(0, 1), (34, 35)],
#     )
#     print(df)
#     print(df2)
#
#     ru = call_nearest(df, df2)
#     print(ru)
#     assert [*ru.Distance] == [1, 3]
#
#
# def test_nearest5():
#     df, df2 = create_dataframes(
#         intervals_a=[(0, 1), (1, 2)],
#         intervals_b=[(0, 1), (3, 4)],
#     )
#     print(df)
#     print(df2)
#
#     ru = call_nearest(df, df2)
#     print(ru)
#     assert [*ru.Distance] == [3, 1]


def call_nearest(df, df2):
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
    dist = pd.DataFrame(dist)
    dist.columns = ["Distance"]
    ru1 = pd.DataFrame(df.loc[idx1]).reset_index(drop=True)
    ru2 = pd.DataFrame(df2.loc[idx2]).reset_index(drop=True)
    ru2.columns = ["C", "S", "E"]
    return pd.concat([ru1, ru2, pd.DataFrame(dist)], axis=1)


def create_dataframes(intervals_a, intervals_b):
    starts_a = [el[0] for el in intervals_a]
    ends_a = [el[1] for el in intervals_a]

    starts_b = [el[0] for el in intervals_b]
    ends_b = [el[1] for el in intervals_b]

    df = pd.DataFrame({"Chromosome": 1, "Start": starts_a, "End": ends_a})
    df2 = pd.DataFrame({"Chromosome": 1, "Start": starts_b, "End": ends_b})
    return df, df2
