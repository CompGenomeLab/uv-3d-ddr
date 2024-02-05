import pandas as pd
import numpy as np
from itertools import combinations

from sklearn.preprocessing import MinMaxScaler
from matplotlib.colors import LogNorm, Normalize
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import seaborn as sns
from coolpuppy import plotpup


def merge(dicts):
    keys = set().union(*dicts)
    return {
        k: [d[k] for d in dicts if k in d]
        for k in keys
    }

def boundary_switch(ins_df_list, sample_names, window=500_000, accepted_range=2, only_strong_base=True):
    
    for i in range(len(ins_df_list)):
        ins_df_list[i]['is_insulating'] = ins_df_list[i][f'boundary_strength_{window}'] > 0
    
    combs_data = []
    mid_point = (accepted_range * 2 + 1) // 2

    combs = set(list(combinations(sample_names+sample_names[::-1] , 2)))
    combs = [comb for comb in combs if comb[0]!=comb[1]]

    for comb in combs:
        
        sample1_index = sample_names.index(comb[0])
        sample2_index = sample_names.index(comb[1])

        print(f"Processing {comb[0]} and {comb[1]}")
        print(f"Sample1 index: {sample1_index}, Sample2 index: {sample2_index}")

        df1, df2 = ins_df_list[sample1_index].copy(), ins_df_list[sample2_index].copy()

        if only_strong_base:
            boundries_1 = df1.loc[df1[f"is_boundary_{window}"] == True]
        else:
            boundries_1 = df1.loc[df1['is_insulating'] == True]
        
        boundries_1 = boundries_1.loc[boundries_1[f"is_bad_bin"] == False ]
        print(f"Total #bins to process of {comb[0]}: {boundries_1.shape[0]}")
        for i, row in boundries_1.iterrows():

            if df2.loc[i, f"is_bad_bin"] == True:
                continue

            local_df = df2.iloc[i - accepted_range:i + accepted_range+1, :]
            local_df_is_boundary = local_df[f"is_insulating"].reset_index(drop=True)
            
            if True in local_df_is_boundary.values:
                where_is_comp_boundary = np.where(local_df_is_boundary == True)[0]
                if len(where_is_comp_boundary) == 1:
                    if where_is_comp_boundary[0] == mid_point:
                        combs_data.append({
                            "sample1": comb[0],
                            "sample2": comb[1],
                            "idx1": i,
                            "idx2": i,
                            "bs1": row[f"boundary_strength_{window}"],
                            "bs2": local_df.iloc[mid_point][f"boundary_strength_{window}"],
                            "case": 'Preserved'
                        })
                    else:
                        position = i-accepted_range+where_is_comp_boundary[0]
                        combs_data.append({
                            "sample1": comb[0],
                            "sample2": comb[1],
                            "idx1": i,
                            "idx2": position,
                            "bs1": row[f"boundary_strength_{window}"],
                            "bs2": local_df.iloc[where_is_comp_boundary[0]][f"boundary_strength_{window}"],
                            "case": 'Shifted'
                        })
                else:
                    assert len(where_is_comp_boundary) > 1
                    if mid_point in where_is_comp_boundary.tolist():
                        for idx_2 in where_is_comp_boundary:
                            if idx_2 == mid_point:
                                position = i
                                combs_data.append({
                                    "sample1": comb[0],
                                    "sample2": comb[1],
                                    "idx1": i,
                                    "idx2": position,
                                    "bs1": row[f"boundary_strength_{window}"],
                                    "bs2": local_df.iloc[mid_point][f"boundary_strength_{window}"],
                                    "case": 'w-Neighbour'
                                })
                            else:
                                position = i-accepted_range+idx_2
                                if local_df.iloc[idx_2]['is_bad_bin'] == True:
                                    continue
                                combs_data.append({
                                    "sample1": comb[0],
                                    "sample2": comb[1],
                                    "idx1": i,
                                    "idx2": position,
                                    "bs1": row[f"boundary_strength_{window}"],
                                    "bs2": local_df.iloc[idx_2][f"boundary_strength_{window}"],
                                    "case": 'New-Neighbour'
                                })
                    else:
                        for idx_2 in where_is_comp_boundary:
                            position = i-accepted_range+idx_2
                            if local_df.iloc[idx_2]['is_bad_bin'] == True:
                                    continue
                            if df1.iloc[i - accepted_range:i + accepted_range+1, :]['is_insulating'].sum() > 1:
                                continue
                            combs_data.append({
                                "sample1": comb[0],
                                "sample2": comb[1],
                                "idx1": i,
                                "idx2": position,
                                "bs1": row[f"boundary_strength_{window}"],
                                "bs2": local_df.iloc[idx_2][f"boundary_strength_{window}"],
                                "case": 'Re-arranged'
                            })
            else:
                combs_data.append({
                    "sample1": comb[0],
                    "sample2": comb[1],
                    "idx1": i,
                    "idx2": -1,
                    "bs1": row[f"boundary_strength_{window}"],
                    "bs2": -1,
                    "case": 'Lost'
                })

    return pd.DataFrame(merge(combs_data))[['sample1', 'idx1', 'bs1', 'sample2', 'idx2', 'bs2', 'case']]

plt.rcParams['figure.dpi'] = 300

def make_bw_df(bw_dict, sigma, use_min_max=True):
    bw_data = {
    'n': [],
    'val': [],
    'name': []
    }

    mm = MinMaxScaler()
    labels = list(range(len(bw_dict.keys())))
    for i, k in enumerate(bw_dict.keys()):
        mean_val = np.nanmean(bw_dict[k], axis=0)
        val = mm.fit_transform(mean_val.reshape(-1, 1)).view(np.ndarray).reshape(-1).tolist() if use_min_max else mean_val
        val = gaussian_filter1d(val, sigma = sigma)
        bw_data['n'].extend(np.arange(len(val)).tolist())
        bw_data['val'].extend(val)
        #bw_data['name'].extend([labels[i]] * len(val))
        bw_data['name'].extend([k] * len(val))

    bw_df = pd.DataFrame(bw_data)
    return bw_df

def plot_bw(bw_dict, sigma=5, ax=None, use_min_max=True, ylabel=None, armlabel=500, xlabel=False):
    bw_df = make_bw_df(bw_dict, sigma, use_min_max=use_min_max)
    l = sns.lineplot(
        data=bw_df,
            x='n', 
            y='val',
            hue='name', 
            ax=ax,
            legend=True)
    end = bw_df['n'].max()
    l.set_xticks([0, end//2, end])

    l.set_xticklabels([f"-{armlabel}kb", 0, f"+{armlabel}kb"])
    l.set_xlabel('Distance from Boundary') if xlabel else l.set_xlabel('')
    l.set_ylabel(ylabel)

def plot_heatmap(
    data, 
    cmp=sns.color_palette("vlag", as_cmap=True), 
    norm=LogNorm, 
    vmin=None, 
    vmax=None, 
    ax=None,
    title=None,
    cbar_ax=None):

    if norm == LogNorm:
        scale = 'log'
        vmin, vmax = plotpup.get_min_max(
            data, 
            sym=False, 
            scale=scale,
            vmin=vmin,
            vmax=vmax,
            )
    
        norm = norm(vmin=vmin, vmax=vmax)

    elif norm == Normalize:
        scale = 'linear'
        vmin, vmax = plotpup.get_min_max(
            data, 
            sym=False, 
            scale=scale,
            vmin=vmin,
            vmax=vmax,
            )
    
        norm = norm(vmin=vmin, vmax=vmax)

    s = sns.heatmap(
        data,
        cmap=cmp,
        norm=norm,
        ax=ax,
        cbar=True if cbar_ax != None else False,
        cbar_ax=None if cbar_ax is None else cbar_ax,
    )

    s.set_title(title)
    s.set_xticks([])
    s.set_xticklabels([])
    s.set_yticks([])
    s.set_yticklabels([])
    
    return norm