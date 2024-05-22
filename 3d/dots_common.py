import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm, Normalize

from coolpuppy import plotpup
import numpy as np



def createRange(
    loc: int,
    accepted_range: int,
) -> list[int]:
    """
    returns a range of values around the given location.
    """
    return list(range(loc - accepted_range, loc + accepted_range + 1, accepted_range))


def findUniqueCommonLoops(
    dots_df1: pd.DataFrame,
    dots_df2: pd.DataFrame,
    merge_common: bool = False,
    accepted_range: int = 10000,
) -> dict[str, pd.DataFrame]:
    """
    returns a dictionary of 3 dataframes:
        - `common_loops`: a list of loops in the first sample that are also in the second sample
        - `specific_loops1`: a list of loops in the first sample that are not in the second sample
        - `specific_loops2`: a list of loops in the second sample that are not in the first sample

    the dataframes have, and are sorted by, the following columns:
        - `chrom1`
        - `start1`
        - `chrom2`
        - `start2`
    """

    df1_comp = (
        dots_df1.sort_values(by=["chrom1", "start1", "chrom2", "start2"])
        .reset_index(drop=True)
        .iloc[:, [0, 1, 3, 4]]
    )
    df2_comp = (
        dots_df2.sort_values(by=["chrom1", "start1", "chrom2", "start2"])
        .reset_index(drop=True)
        .iloc[:, [0, 1, 3, 4]]
    )
    df1_merged = [
        zipped
        for zipped in zip(
            df1_comp["chrom1"],
            df1_comp["start1"],
            df1_comp["chrom2"],
            df1_comp["start2"],
        )
    ]
    df2_merged = [
        zipped
        for zipped in zip(
            df2_comp["chrom1"],
            df2_comp["start1"],
            df2_comp["chrom2"],
            df2_comp["start2"],
        )
    ]

    # find common loops. common loops should have the same start and end in both samples, with a tolerance of 10kb
    common_1 = list()
    common_2 = list()

    for loop in df1_merged:
        range1 = createRange(loop[1], accepted_range)
        range2 = createRange(loop[3], accepted_range)
        accepted_loops = list(
            filter(
                lambda x: x[0] == loop[0]
                and x[1] in range1
                and x[2] == loop[2]
                and x[3] in range2,
                df2_merged,
            )
        )

        if len(accepted_loops) > 0:
            common_1.append(loop)
    
    for loop in df2_merged:
        range1 = createRange(loop[1], accepted_range)
        range2 = createRange(loop[3], accepted_range)
        accepted_loops = list(
            filter(
                lambda x: x[0] == loop[0]
                and x[1] in range1
                and x[2] == loop[2]
                and x[3] in range2,
                df1_merged,
            )
        )

        if len(accepted_loops) > 0:
            common_2.append(loop)

    # find specific loops
    specific_loops1 = list(
        filter(lambda x: x not in common_1 
               and x not in common_2, 
               df1_merged)
    )  # loops in sample 1 that are not in sample 2

    specific_loops2 = list(
        filter(lambda x: x not in common_1 
               and x not in common_2, 
               df2_merged)
    )  # loops in sample 2 that are not in sample 1

    assert len(common_1) + len(specific_loops1) == len(df1_merged)
    assert len(common_2) + len(specific_loops2) == len(df2_merged)

    loops: dict[str, pd.DataFrame] = {}

    loops["common_loops1"] = pd.DataFrame(
        common_1, columns=["chrom1", "start1", "chrom2", "start2"]
    )
    loops["common_loops2"] = pd.DataFrame(
        common_2, columns=["chrom1", "start1", "chrom2", "start2"]
    )
    loops["specific_loops1"] = pd.DataFrame(
        specific_loops1, columns=["chrom1", "start1", "chrom2", "start2"]
    )
    loops["specific_loops2"] = pd.DataFrame(
        specific_loops2, columns=["chrom1", "start1", "chrom2", "start2"]
    )

    if merge_common:
        loops["common_loops"] = pd.concat(
            [loops["common_loops1"], loops["common_loops2"]]
        ).drop_duplicates()
        del loops["common_loops1"]
        del loops["common_loops2"]
        return loops
    else:
        return loops


def findUniqueCommonAnchors(
    dots_df1: pd.DataFrame,
    dots_df2: pd.DataFrame,
    merge_common: bool = False,
    accepted_range: int = 10000,
) -> dict[str, pd.DataFrame]:
    """
    returns a dictionary of 4 dataframes:
        - `common_anchors1`: a list of anchors in the first sample that are also anchors in the second sample
        - `common_anchors2`: a list of anchors in the second sample that are also anchors in the first sample
        - `specific_anchors1`: a list of anchors in the first sample that are not in the second sample
        - `specific_anchors2`: a list of anchors in the second sample that are not in the first sample

    the dataframes have, and are sorted by, the following columns:
        - `chrom`
        - `start`
    """

    anchors: dict[str, pd.DataFrame] = {}

    df1_comp = (
        dots_df1.sort_values(by=["chrom1", "start1", "chrom2", "start2"])
        .reset_index(drop=True)
        .iloc[:, [0, 1, 3, 4]]
    )
    df2_comp = (
        dots_df2.sort_values(by=["chrom1", "start1", "chrom2", "start2"])
        .reset_index(drop=True)
        .iloc[:, [0, 1, 3, 4]]
    )
    df1_anch1_merged = [
        zipped for zipped in zip(df1_comp["chrom1"], df1_comp["start1"])
    ]
    df1_anch2_merged = [
        zipped for zipped in zip(df1_comp["chrom2"], df1_comp["start2"])
    ]
    df2_anch1_merged = [
        zipped for zipped in zip(df2_comp["chrom1"], df2_comp["start1"])
    ]
    df2_anch2_merged = [
        zipped for zipped in zip(df2_comp["chrom2"], df2_comp["start2"])
    ]

    # find common anchors. common anchors should have the same start and end in both samples, with a tolerance of 10kb
    common_anchors1 = list()

    for anch1 in df1_anch1_merged + df1_anch2_merged:
        range1 = createRange(anch1[1], accepted_range)
        accepted_anchors = list(
            filter(
                lambda x: x[0] == anch1[0] and x[1] in range1,
                df2_anch1_merged + df2_anch2_merged,
            )
        )

        if len(accepted_anchors) > 0:
            common_anchors1.append(anch1)

    common_anchors2 = list()
    
    for anch2 in df2_anch1_merged + df2_anch2_merged:
        range2 = createRange(anch2[1], accepted_range)
        accepted_anchors = list(
            filter(
                lambda x: x[0] == anch2[0] and x[1] in range2,
                df2_anch2_merged + df2_anch2_merged,
            )
        )

        if len(accepted_anchors) > 0:
            common_anchors2.append(anch2)

    # find specific loops
    specific_anchors1 = list(
        filter(lambda x: x not in common_anchors1, df1_anch1_merged + df1_anch2_merged)
    )  # loops in sample 1 that are not in sample 2

    specific_anchors2 = list(
        filter(lambda x: x not in common_anchors2, df2_anch1_merged + df2_anch2_merged)
    )  # loops in sample 2 that are not in sample 1

    assert len(specific_anchors1) + len(common_anchors1) == len(dots_df1) * 2
    assert len(specific_anchors2) + len(common_anchors2) == len(dots_df2) * 2

    anchors["common_anchors1"] = pd.DataFrame(
        common_anchors1, columns=["chrom", "start"]
    )
    anchors["common_anchors2"] = pd.DataFrame(
        common_anchors2, columns=["chrom", "start"]
    )
    anchors["specific_anchors1"] = pd.DataFrame(
        specific_anchors1, columns=["chrom", "start"]
    )
    anchors["specific_anchors2"] = pd.DataFrame(
        specific_anchors2, columns=["chrom", "start"]
    )

    if merge_common:
        anchors["common_anchors"] = pd.concat(
            [anchors["common_anchors1"], anchors["common_anchors2"]]
        ).drop_duplicates()
        del anchors["common_anchors1"]
        del anchors["common_anchors2"]
        return anchors
    else:
        return anchors
