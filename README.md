## UV-DDR - 3D analysis 
___
This repository contains exploratory and hierarchical genome architecture analysis code, and research code for dataset generation, training, inference for GNN algorithm, as described in the manuscript titled **"UV-induced reorganization of 3D genome mediates DNA damage response"**.
___
This repo contains two sections:
1. Python notebooks under **3D folder** contain genome architecture analysis:
   1. **boundaries** contains analysis code with respect to TADs and boundaries/insulation.
   2. **compartments** contains analysis code with respect to compartments, saddle plots saddle data analysis, compartment strength and also compartment switch analysis code.
   3. **RNA** contains R scripts for DESeq time series/DEG analysis and for clustering, also contains a notebook for GNN, chromatin loop, TADs analysis with respect to RNA-Seq data.
   4.  **unibind** contains shell scripts to run TFBS differential enrichment analysis on regions of GNN comparison and chromatin loop anchors, also contains a notebook for visualization.
   5.  **distance_decay.ipynb** notebook is for the analysis for distance decay and powerlaw estimations.
   6.  **gnn_results_analyze.ipynb** notebook is for the analysis of comparisons from GNN results.
   7.  **loops.ipynb** notebook is for loop calling, loop and anchor analysis.
   8.  **repair_ops.ipynb** notebook is for downstream processing of XR-Seq and Damage-Seq (and simulations) bed files generated with xr-ds-snakemake Snakemake pipeline, and converting to bigwigs.
   9.  **scc.ipynb** notebook is for the analysis of the Hi-C samples with HiCRep, and MDS code.
2. **GNN folder** contain dataset generation, training and inference code for the workflow described in the manuscript.
   1. Follow **save_obs_exp_qt.ipynb** in order to perform preprocessing steps for the matrices (Both megamap, and Hi-C matrices of our own).
   2. Follow **training_dataset.ipynb** in order to generate in-memory datasets to sample target and query graphs from, during training.
   3. Use **train.sh** script to start training.
   4. Use **predict.sh** for the inference of the model on the contact matrices to be compared. 