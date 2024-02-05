## UV-DDR - 3D analysis 
___
This repository contains exploratory and hierarchical genome architecture analysis code, and research code for dataset generation, training, inference for GNN algorithm, as described in the manuscript titled **"UV-induced reorganization of 3D genome mediates DNA damage response"**.
___
This repo contains two sections:
1. Python notebooks contain genome architecture analysis:
   1. 
   2. 
2. GNN folder contain dataset generation, training and inference code for the workflow described in the manuscript.
   1. Follow **save_obs_exp_qt.ipynb** in order to perform preprocessing steps for the matrices (Both megamap, and Hi-C matrices of our own).
   2. Follow **training_dataset.ipynb** in order to generate in-memory datasets to sample target and query graphs from, during training.
   3. Use **train.sh** script to start training.
   4. Use **predict.sh** for the inference of the model on the contact matrices to be compared. 