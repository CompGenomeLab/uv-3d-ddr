# UV-DDR-3D analysis 
___
This repository contains exploratory and hierarchical genome architecture analysis code, and research code for dataset generation, training, inference for GNN algorithm, as described in the manuscript titled **"UV-induced reorganization of 3D genome mediates DNA damage response"**.

All analysis was performed on a Linux system (Ubuntu 22.04.4 LTS). To reproduce the analysis and the figures in the manuscript, code is particularly deposited in the form of Jupyter notebooks (Python), while also preserving code outputs.
___
This code repo contains two sections:
1. Python notebooks under **3D folder** contain genome architecture analysis:
   1. **boundaries/tads_rework.ipynb** contains analysis code with respect to TADs and boundaries/insulation (Fig3 b,c,d,e,f,g, Supp.Fig2a,b,c,d,e).
   2. **boundaries/bs_raincloud.ipynb** raincloud plot of boundary strengths per sample (Fig3 a).
   3. **compartments/comp_ev_check.ipynb** contains the correlation analysis of compartment profiles (EV) and phasing tracks (Supp Table).
   4. **compartments/comp_switch.ipynb** contains the analysis code for correlation of compartment profiles with respect to timecourse DEG clusters. (Fig2b,c,d) 
   5. **compartments/compartments.ipynb** contains analysis code with respect to compartments, saddle plots saddle data analysis, compartment strength. (Fig1d,e,f, Supp.Fig1a)
   6. **compartments/ev-change_deseq.ipynb & ev-tpm.ipynb** EV change profiles and DEGs between samples, and TPMs per compartment profiles (Supp.Fig1c,d).
   7. **compartments/raincloud.ipynb** Top 20% compartment profiles interactions change (Fig1g).
   8. **RNA** contains R scripts for DESeq time series/DEG analysis and **rna_3D.ipynb** for annotating GNN analysis with respect to RNA-Seq data (Fig6e,f). Also, source data files of coseq results for expression clustering related analysis. 
   9. **unibind** contains shell scripts to run TFBS differential enrichment analysis on regions of GNN comparison and chromatin loop anchors, also contains a notebook for visualization (Fig4h, Fig6c,d).
   10. **distance_decay.ipynb** notebook is for the analysis for distance decay and powerlaw estimations (Fig1 b,c).
   11. **gnn_analyse/gnn_results_analyse.ipynb** notebook is for the analysis of comparisons from GNN results (Fig5c,6a,b).
   12. **gnn_analyse/gnn_dnase_faire.ipynb** notebook is for the annotation of comparison groups with DNAse-Seq and FAIRE-Seq signals (Fig5d).
   13. **loops/loops.ipynb** notebook is for loop calling, loop and anchor analysis (Fig4a,b,c,d,e,f, Supp.Fig3).
   14. **loops/loops_rna.ipynb** notebook is for annotating specific loops with DEGs (Fig4g).
   15. **repair_ops.ipynb** notebook is for downstream processing of XR-Seq and Damage-Seq (and simulations) bed files generated with xr-ds-snakemake Snakemake pipeline, and converting to bigwigs.
2. **GNN folder** contain dataset generation, training and inference code for the GNN workflow described in the manuscript.
   1. Follow **save_obs_exp_qt.ipynb** in order to perform preprocessing steps for the matrices (Both megamap, and Hi-C matrices of our own). (You can easily modify for your own **.mcool/.cool** files.)
   2. Follow **training_dataset.ipynb** in order to generate in-memory datasets to sample target and query graphs from, for training the GNN model. (You can easily modify for your own data.)
   3. Use **train.sh** script to start training. (An example slurm batch script is provided)
   4. Use **predict.sh** for the inference of the model on the contact matrices to be compared. (GNN Model checkpoint used in the analysis can be found under chekpoints/last.ckpt, also comparison raw results between Control and 12min as reference point can be found as **results_t0_t12.tsv** ) (An example slurm batch script is provided)

# Input Data

Processed files in order to reproduce the analysis code were made available on GEO database. Multi-resolution contact matrices in **mcool** format can be found with the identifier [GSE268350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268350) for Hi-C data, and **TPM quantification** data can be found with the identifier [GSE268349](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268349) for RNA-seq data.

Also, all raw sequencing data generated in this study have been submitted to the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/bioproject/), under accession number [PRJNA1053630](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1053630). All sequencing data were first preprocessed to QC with [fastp](https://github.com/OpenGene/fastp) with flag **--detect_adapter_for_pe**. Hi-C contacts matrices were generated (**.hic** format) with Juicer v1.6 pipeline, with default parameters. **.hic** format (q30 filtered) has been converted to **.mcool** format with **hicConvertFormat** subprogram of [HiCExplorer](https://github.com/deeptools/HiCExplorer). Matrix balancing has been performed with **balance** subprogram of [cooler](https://github.com/open2c/cooler) with parameters **--cis-only --max-iters 1000 --blacklist**, where blacklist file is [here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz). Accordingly, all of the 3D genome research code presented in this code repository operates on **.cool/.mcool** data format.

## Dependencies

Below is the list of required Python packages, installing via Mamba is advised.

- python=3.11.0
- cooler=0.9.2
- cooltools=0.5.4
- powerlaw
- bioframe=0.4.1
- coolpuppy=1.1.0
- scikit-learn=1.3.0
- pytorch=2.0.1
- pytorch-geometric=2.3.1
- pytorch-lightning=2.0.6
- gseapy=1.0.6

# Contact
Please email [vogulcan@sabanciuniv.edu](vogulcan@sabanciuniv.edu) or raise an issue in the github repository with any questions about installation or usage.
