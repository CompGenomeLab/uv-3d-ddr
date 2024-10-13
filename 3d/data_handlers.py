# import and center bigwigs mean to 1 genome-wide

import bbi
import bioframe
import numpy as np

chrom_sizes = bioframe.fetch_chromsizes("hg38")[:23]
def get_mean_bw(bw, chrom_sizes):
    for chrom in chrom_sizes.index:
        start, end = 0, chrom_sizes.loc[chrom] // 10000 * 10000
        yield np.nanmean(bbi.fetch(bw, chrom, start, end, bins=chrom_sizes.loc[chrom] // 1000))


BWS = {
    'xr-64-over-sim': '/home/carlos/oldies/manuscripts/notebooks/bws/XR_64_real_over_sim_res1000.bw',
    'xr-cpd-over-sim': '/home/carlos/oldies/manuscripts/notebooks/bws/XR_CPD_real_over_sim_res1000.bw',
    'ds-64-over-sim': '/home/carlos/oldies/manuscripts/notebooks/bws/DS_64_real_over_sim_res1000.bw',
    'ds-cpd-over-sim': '/home/carlos/oldies/manuscripts/notebooks/bws/DS_CPD_real_over_sim_res1000.bw',
    'rep-eff-64': '/home/carlos/oldies/manuscripts/notebooks/bws/XR_64_rep_eff_res1000.bw',
    'rep-eff-cpd': '/home/carlos/oldies/manuscripts/notebooks/bws/XR_CPD_rep_eff_res1000.bw',
    'ctcf': '/home/carlos/oldies/manuscripts/notebooks/bws/ctcf.foi.bigWig',
    'e1_0': '/home/carlos/oldies/manuscripts/review/compartments/ev_bigwigs/E1_0.res.25kb.bw',
    'e1_12': '/home/carlos/oldies/manuscripts/review/compartments/ev_bigwigs/E1_12.res.25kb.bw',
    'rna_0': '/home/carlos/oldies/manuscripts/review/compartments/merged_t00.bw',
    'rna_12': '/home/carlos/oldies/manuscripts/review/compartments/merged_t12.bw',
    'heterochromatin': '/home/carlos/Downloads/ENCFF891XLY.bigWig',
    'repressive':'/home/carlos/Downloads/ENCFF614HNF.bigWig',
    'enhancer': '/home/carlos/Downloads/ENCFF430ZMK.bigWig',
    'promoter':'/home/carlos/Downloads/ENCFF194XTD.bigWig',
    'dnase': "/home/carlos/oldies/manuscripts/notebooks/bws/dnase.rpkm.bw",
    # 'smc3': '/home/carlos/Downloads/ENCFF971PWK.bigWig',
    # 'rad21': '/home/carlos/Downloads/ENCFF357XMG.bigWig',
    'faire': '/home/carlos/Downloads/ENCFF000TKE_hg38.bigWig',
    }

MEANS_BW = {}
for k in BWS:
    curr_mean = np.mean(list(get_mean_bw(BWS[k], chrom_sizes)))
    MEANS_BW[k] = curr_mean

means_norm = True
interval_size = 1000
bws_per_key = {}
for k in BWS.keys():
    per_chr = []
    for chrName in chrom_sizes.index:
        chrLen = chrom_sizes.loc[chrName]
        start, end = 0, chrLen // interval_size * interval_size
        n_bins = end // interval_size
        arr_chr = bbi.fetch(
            BWS[k],
            chrName,
            start,
            end,
            bins=n_bins,
        )
        per_chr.append(arr_chr / MEANS_BW[k]) if means_norm else per_chr.append(arr_chr)
    bws_per_key[k] = per_chr

def _bw_data(
    curr_DF, 
    flank, 
    n_bins_start, 
    n_bins_end, 
    n_bins_main, 
    mean_norm=True):

    start_data = {k: bbi.stackup(v,
                            curr_DF.chrom,
                            curr_DF.start - flank,
                            curr_DF.start,
                            bins=n_bins_start) 
                            for k, v in BWS.items()
                            }
    main_data = {k: bbi.stackup(v,
                            curr_DF.chrom,
                            curr_DF.start,
                            curr_DF.end,
                            bins=n_bins_main) 
                            for k, v in BWS.items()
                            }
    end_data = {k: bbi.stackup(v,
                            curr_DF.chrom,
                            curr_DF.end,
                            curr_DF.end + flank,
                            bins=n_bins_end)
                            for k, v in BWS.items()
                            }

    if mean_norm:
        for k in start_data:
            start_data[k] = start_data[k] / MEANS_BW[k]
        for k in main_data:
            main_data[k] = main_data[k] / MEANS_BW[k]
        for k in end_data:
            end_data[k] = end_data[k] / MEANS_BW[k]
    return start_data, main_data, end_data

def region_per_df_wrapper(df, arr_k_dict, num_bins, interval_size=1000, flank=500_000, func=np.nanmean):
    regions = []
    for i, row in df.iterrows():
        chrom, start, end = row['chrom'], row['start'], row['end']
        chrom_idx = chrom_sizes.index.get_loc(chrom)
        start, end = start - flank, end + flank

        curr_arr = arr_k_dict[chrom_idx][start//interval_size:end//interval_size].flatten()
        if len(np.where(curr_arr == np.nan)[0]) >= len(curr_arr * 10 / 100):
            print("nan")
            continue

        bin_edges = np.linspace(0, len(curr_arr), num_bins + 1, dtype=int)
        binned_array = np.array([func(curr_arr[bin_edges[i]:bin_edges[i+1]]) for i in range(num_bins)])
        regions.append(binned_array)

    return np.array(regions)