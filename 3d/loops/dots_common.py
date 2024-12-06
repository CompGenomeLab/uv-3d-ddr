import pandas as pd
import numpy as np
import cooltools
import bbi
from scipy.ndimage import gaussian_filter1d

def createRange(
    loc: int,
    accepted_range: int,
    resolution: int = 10000,
) -> list[int]:
    """
    returns a range of values around the given location.
    """
    return list(range(loc - accepted_range, loc + resolution + accepted_range + 1, resolution))

def merge_left_right_anchors(df): # from loop df 
    left_anchors_df = df.iloc[:, :3].copy()
    right_anchors_df = df.iloc[:, 3:].copy()
    right_anchors_df.columns = left_anchors_df.columns
    return pd.concat([left_anchors_df, right_anchors_df], ignore_index=True).rename(
        columns={'chrom1': 'chrom', 'start1': 'start', 'end1': 'end'}
    ).drop_duplicates()

def findUniqueCommonLoops(
    dots_df1: pd.DataFrame,
    dots_df2: pd.DataFrame,
    merge_common: bool = False,
    accepted_range: int = 10000,
    resolution: int = 10000,
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
        range1 = createRange(loop[1], accepted_range, resolution=resolution)
        range2 = createRange(loop[3], accepted_range, resolution=resolution)
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
        range1 = createRange(loop[1], accepted_range, resolution=resolution)
        range2 = createRange(loop[3], accepted_range, resolution=resolution)
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
    resolution: int = 10000,
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

    df1_anchors = merge_left_right_anchors(dots_df1.copy())
    df2_anchors = merge_left_right_anchors(dots_df2.copy())

    df1_anch_zipped = [zipped for zipped in zip(df1_anchors.iloc[:, :2]["chrom"], df1_anchors.iloc[:, :2]["start"])]
    df2_anch_zipped = [zipped for zipped in zip(df2_anchors.iloc[:, :2]["chrom"], df2_anchors.iloc[:, :2]["start"])]

    # find common anchors. common anchors should have the same start and end in both samples, with a tolerance of 10kb
    common_anchors1 = list()

    for anch1 in df1_anch_zipped:
        range1 = createRange(anch1[1], accepted_range, resolution=resolution)
        accepted_anchors = list(
            filter(
                lambda x: x[0] == anch1[0] and x[1] in range1,
                df2_anch_zipped,
            )
        )

        if len(accepted_anchors) > 0:
            common_anchors1.append(anch1)

    common_anchors2 = list()
    
    for anch2 in df2_anch_zipped:
        range2 = createRange(anch2[1], accepted_range, resolution=resolution)
        accepted_anchors = list(
            filter(
                lambda x: x[0] == anch2[0] and x[1] in range2,
                df1_anch_zipped,
            )
        )

        if len(accepted_anchors) > 0:
            common_anchors2.append(anch2)

    # find specific loops
    specific_anchors1 = list(
        filter(lambda x: x not in common_anchors1, df1_anch_zipped)
    )  # loops in sample 1 that are not in sample 2

    specific_anchors2 = list(
        filter(lambda x: x not in common_anchors2, df2_anch_zipped)
    )  # loops in sample 2 that are not in sample 1

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
        anchors_dedup = {}
        for key, value in anchors.items():
            anchors_dedup[key] = value.drop_duplicates()
        return anchors_dedup

def _2_step_loop_separation(
        df_1, df_2, 
        clr_i, clr_j, 
        expected_i, expected_j, 
        view_df, 
        threshold = 1, 
        merge_common = True, 
        accepted_range = 50_000, 
        resolution = 10_000,
        nproc = 4,
        weight_name = "weight",):
    
    # step 1
    # t0 vs t12
    loops = findUniqueCommonLoops(df_1, df_2, merge_common=False, accepted_range=50_000, resolution=resolution)

    # step 2
    # threshold is log2 fold change

    apa_flank = 0 # we are only interested in the central pixel (o/e interaction frequency of anchor1 and anchor2)

    new_specific_dict = {}
    ratio_data = []
    for i, (k, v) in enumerate(loops.items()):
        if "specific" in k:

            df = v.copy()
            df['end1'] = df['start1'] + resolution
            df['end2'] = df['start2'] + resolution
            #df = df.head(10)
            pup_0 = cooltools.pileup(clr_i, df, view_df=view_df, expected_df=expected_i, flank=apa_flank, nproc=nproc, clr_weight_name=weight_name)
            pup_1 = cooltools.pileup(clr_j, df, view_df=view_df, expected_df=expected_j, flank=apa_flank, nproc=nproc, clr_weight_name=weight_name)

            pup_0_cs = pup_0.flatten()
            pup_1_cs = pup_1.flatten()

            def _fix_nans(diff_arr, base_arr):
                diff_arr[np.isnan(diff_arr)] = base_arr[np.isnan(diff_arr)]
                return diff_arr

            if k == "specific_loops1":
                diff = pup_0_cs / pup_1_cs
                diff = _fix_nans(diff, pup_0_cs)
            elif k == "specific_loops2":
                diff = pup_1_cs / pup_0_cs
                diff = _fix_nans(diff, pup_1_cs)
            else:
                print("WARNING: not a specific loop")

            # omit nan
            diff = diff[~np.isnan(diff)]

            cs_l2fc = np.log2(diff)
            # l2fc threshold
            ratio_data.append((cs_l2fc, threshold))
            idx1 = np.where(cs_l2fc >= threshold)[0]
            
            new_df = df.iloc[idx1].copy()
            new_df['ratio'] = diff[idx1]
            new_specific_dict[k] = new_df.copy()

    # merge common loops1 to specific loops1 from previous step
    # merge common loops2 to specific loops2 from previous step

    dots_df_reiter_ = []
    for i, (k, v) in enumerate(new_specific_dict.items()):
        if 'specific' not in k:
            continue
        
        specific_df = v.iloc[:, :6].copy()
        specific_df = specific_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]

        current_common_loops_key = f"common_loops{i+1}"

        print(f"End of step1: Merging {current_common_loops_key} with {k}")

        common_df = loops[current_common_loops_key].copy()
        common_df['end1'] = common_df['start1'] + resolution
        common_df['end2'] = common_df['start2'] + resolution
        common_df = common_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
        frames = [specific_df, common_df]
        df_out = pd.concat(frames)
        df_out.drop_duplicates(inplace=True)
        dots_df_reiter_.append(df_out)

    loops = findUniqueCommonLoops(dots_df_reiter_[0], dots_df_reiter_[1], merge_common=merge_common, accepted_range=accepted_range, resolution=resolution)
    anchors = findUniqueCommonAnchors(dots_df_reiter_[0], dots_df_reiter_[1], merge_common=merge_common, accepted_range=accepted_range, resolution=resolution)

    for i, (k, v) in enumerate(loops.items()):
        loops[k] = loops[k].drop_duplicates(subset=['chrom1', 'start1', 'chrom2', 'start2'])

    for i, (k, v) in enumerate(anchors.items()):
        anchors[k] = anchors[k].drop_duplicates(subset=['chrom', 'start'])

    return loops, anchors#, dots_df_reiter_[0], dots_df_reiter_[1]

def make_loops_df(loops_dict, resolution=10000):
    dfs = []
    for k,v in loops_dict.items():
        df = v.copy()
        df["type"] = [k] * df.shape[0]
        df['end1'] = df['start1'] + resolution
        df['end2'] = df['start2'] + resolution
        df['length'] = df['start2'] - df['end1']
        dfs.append(df)

    df = pd.concat(dfs)[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'type', 'length']]
    return df

def make_anchors_df(anchors_dict, resolution=10000):
    dfs = []
    for k,v in anchors_dict.items():
        df = v.copy()
        df['end'] = df['start'] + resolution
        df["type"] = [k] * df.shape[0]
        dfs.append(df)

    df = pd.concat(dfs)[['chrom', 'start', 'end', 'type']]
    return df

def _bw_data_pe(
    curr_DF, 
    flank, 
    n_bins_start, 
    n_bins_end, 
    n_bins_main, 
    BWS,
    MEANS_BW=None,
    ):

    left_flank_data = {k: bbi.stackup(v,
                            curr_DF.chrom1,
                            curr_DF.start1 - flank,
                            curr_DF.start1,
                            bins=n_bins_start) 
                            for k, v in BWS.items()
                            }
    loop_body_data = {k: bbi.stackup(v,
                            curr_DF.chrom1,
                            curr_DF.start1,
                            curr_DF.end2,
                            bins=n_bins_main) 
                            for k, v in BWS.items()
                            }
    right_flank_data = {k: bbi.stackup(v,
                            curr_DF.chrom1,
                            curr_DF.end2,
                            curr_DF.end2 + flank,
                            bins=n_bins_end)
                            for k, v in BWS.items()
                            }

    if MEANS_BW is not None:
        for k in left_flank_data:
            left_flank_data[k] = left_flank_data[k] / MEANS_BW[k]
        for k in loop_body_data:
            loop_body_data[k] = loop_body_data[k] / MEANS_BW[k]
        for k in right_flank_data:
            right_flank_data[k] = right_flank_data[k] / MEANS_BW[k]
    return left_flank_data, loop_body_data, right_flank_data


def _bw_data_se(
    curr_DF, 
    flank, 
    n_bins_start, 
    n_bins_end, 
    n_bins_main, 
    BWS,
    MEANS_BW=None,
    ):

    left_flank_data = {k: bbi.stackup(v,
                            curr_DF.chrom,
                            curr_DF.start - flank,
                            curr_DF.start,
                            bins=n_bins_start) 
                            for k, v in BWS.items()
                            }
    loop_body_data = {k: bbi.stackup(v,
                            curr_DF.chrom,
                            curr_DF.start,
                            curr_DF.end,
                            bins=n_bins_main) 
                            for k, v in BWS.items()
                            }
    right_flank_data = {k: bbi.stackup(v,
                            curr_DF.chrom,
                            curr_DF.end,
                            curr_DF.end + flank,
                            bins=n_bins_end)
                            for k, v in BWS.items()
                            }

    if MEANS_BW is not None:
        for k in left_flank_data:
            left_flank_data[k] = left_flank_data[k] / MEANS_BW[k]
        for k in loop_body_data:
            loop_body_data[k] = loop_body_data[k] / MEANS_BW[k]
        for k in right_flank_data:
            right_flank_data[k] = right_flank_data[k] / MEANS_BW[k]
    return left_flank_data, loop_body_data, right_flank_data

def get_mean_bw(bw, chrom_sizes):
    for chrom in chrom_sizes.index:
        start, end = 0, chrom_sizes.loc[chrom] // 10000 * 10000
        yield np.nanmean(bbi.fetch(bw, chrom, start, end, bins=chrom_sizes.loc[chrom] // 1000))

def _make_bw_dict_to_df_avg(bw_dict, sigma=None):
    dfs = []
    for k, v in bw_dict.items():
        df_template = {
            "n_range": [],
            "data": [],
            "type": [],
        }

        vec = np.nanmean(v, axis=0)
        total_bins = vec.shape[0]
        df_template["n_range"] += np.arange(total_bins).tolist()
        df_template["data"] += vec.tolist() if sigma == None else gaussian_filter1d(vec, sigma=sigma).flatten().tolist()
        df_template["type"] += [k] * total_bins
        dfs.append(pd.DataFrame(df_template))
    return pd.concat(dfs)

from matplotlib import patches

def add_scores(scores_type, score_, ax, mtx_size, corner_size, font_scale = 2, height = 6):
    
    if scores_type == "central_score":

        ax.text(
                s=f"{score_:.3g}",
                y=0.95,
                x=0.05,
                ha="left",
                va="top",
                size=font_scale * (4.94 + height),
                transform=ax.transAxes,
        )

    elif scores_type == "corner_scores":
        
        # top left
        ax.text(
                s=f"{score_[0]:.3g}",
                y=0.95,
                x=0.05,
                ha="left",
                va="top",
                size=font_scale * (4.94 + height),
                transform=ax.transAxes,
        )
        
        # add patch
        rect = patches.Rectangle((0, 0), corner_size, corner_size, linewidth=2, edgecolor='#465775', facecolor='none')
        ax.add_patch(rect)

        # top right
        ax.text(
                s=f"{score_[1]:.3g}",
                y=0.95,
                x=0.95,
                ha="right",
                va="top",
                size=font_scale * (4.94 + height),
                transform=ax.transAxes,
        )

        # add patch
        rect = patches.Rectangle((mtx_size - corner_size, 0), corner_size, corner_size, linewidth=2, edgecolor='#465775', facecolor='none')
        ax.add_patch(rect)


        # bottom left
        ax.text(
                s=f"{score_[2]:.3g}",
                y=0.05,
                x=0.05,
                ha="left",
                va="bottom",
                size=font_scale * (4.94 + height),
                transform=ax.transAxes,
        )
        
        # add patch
        rect = patches.Rectangle((0, mtx_size - corner_size), corner_size, corner_size, linewidth=2, edgecolor='#465775', facecolor='none')
        ax.add_patch(rect)

        # bottom right
        ax.text(
                s=f"{score_[3]:.3g}",
                y=0.05,
                x=0.95,
                ha="right",
                va="bottom",
                size=font_scale * (4.94 + height),
                transform=ax.transAxes,
        )

        # add patch
        rect = patches.Rectangle((mtx_size - corner_size, mtx_size - corner_size), corner_size, corner_size, linewidth=2, edgecolor='#465775', facecolor='none')
        ax.add_patch(rect)