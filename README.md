# UV-DDR-3D analysis 
___
This repository contains exploratory and hierarchical genome architecture analysis code, and research code for dataset generation, training, inference for GNN algorithm, as described in the manuscript titled **"UV-induced reorganization of 3D genome mediates DNA damage response"**.

Preprint is available at [bioRxiv](https://doi.org/10.1101/2024.05.27.595922).

All analysis was performed on a Linux system (Ubuntu 22.04.4 LTS). To reproduce the analysis and the figures in the manuscript, code is particularly deposited in the form of Jupyter notebooks (Python), while also preserving code outputs.

Supplementary tables referenced in the manuscript can be found as **supp_table.xlsx**.
___
This code repo contains two sections:
1. Python notebooks under **3D folder** contain genome architecture analysis:
   1. **boundaries** contains analysis code with respect to TADs and boundaries/insulation.
   2. **compartments** contains analysis code with respect to compartments, saddle plots saddle data analysis, compartment strength, EV correlation to phasing tracks and also compartment correlation analysis code.
   3. **RNA** contains R scripts for DESeq time series/DEG analysis and for clustering, also contains a notebook for GNN, chromatin loop, TADs analysis with respect to RNA-Seq data.
   4.  **unibind** contains shell scripts to run TFBS differential enrichment analysis on regions of GNN comparison and chromatin loop anchors, also contains a notebook for visualization.
   5.  **distance_decay.ipynb** notebook is for the analysis for distance decay and powerlaw estimations.
   6.  **gnn_results_analyze.ipynb** notebook is for the analysis of comparisons from GNN results.
   7.  **loops.ipynb** notebook is for loop calling, loop and anchor analysis.
   8.  **repair_ops.ipynb** notebook is for downstream processing of XR-Seq and Damage-Seq (and simulations) bed files generated with xr-ds-snakemake Snakemake pipeline, and converting to bigwigs.
2. **GNN folder** contain dataset generation, training and inference code for the workflow described in the manuscript.
   1. Follow **save_obs_exp_qt.ipynb** in order to perform preprocessing steps for the matrices (Both megamap, and Hi-C matrices of our own). (You can easily modify for your own **.mcool/.cool** files.)
   2. Follow **training_dataset.ipynb** in order to generate in-memory datasets to sample target and query graphs from, during training. (You can easily modify for your own data.)
   3. Use **train.sh** script to start training. (An example slurm batch script is provided)
   4. Use **predict.sh** for the inference of the model on the contact matrices to be compared. (GNN Model checkpoint used in the analysis can be found under chekpoints/last.ckpt, also comparison raw results between Control and 12min as reference point can be found as **results_t0_t12.tsv** ) (An example slurm batch script is provided)

# Input Data

Processed files in order to reproduce the analysis code were made available on GEO database. Multi-resolution contact matrices in **mcool** format can be found with the identifier [GSE268350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268350) for Hi-C data, and **TPM quantification** data can be found with the identifier [GSE268349](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268349) for RNA-seq data.

Also, all raw sequencing data generated in this study have been submitted to the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/bioproject/), under accession number [PRJNA1053630](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1053630). All sequencing data were first preprocessed to QC with [fastp](https://github.com/OpenGene/fastp) with flag **--detect_adapter_for_pe**. Hi-C contacts matrices were generated (**.hic** format) with Juicer v1.6 pipeline, with default parameters. **.hic** format (q30 filtered) has been converted to **.mcool** format with **hicConvertFormat** subprogram of [HiCExplorer](https://github.com/deeptools/HiCExplorer). Matrix balancing has been performed with **balance** subprogram of [cooler](https://github.com/open2c/cooler) with parameters **--cis-only --max-iters 1000 --blacklist**, where blacklist file is [here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz). Accordingly, all of the 3D genome research code presented in this code repository operates on **.cool/.mcool** data format.

## Dependencies

Below is the list of required packages, installing via Mamba is advised.

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
