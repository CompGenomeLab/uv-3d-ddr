#!/bin/bash
#SBATCH --account=investor

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:a100:1

#SBATCH --qos=short_investor
#SBATCH --partition=short_investor
#SBATCH --time=2:00:00
#SBATCH --mem=60G

source /cta/users/vkaya/miniconda3/bin/activate gnn

checkpoint_path="/cta/users/vkaya/gnn/diff-gnn/experiments/tag3/version_0/checkpoints/last.ckpt"

srun --unbuffered python predict.py \
    --checkpoint_path $checkpoint_path \
    --sample1 t0_q30 \
    --sample2 t12_q30 \
    --npz_path /cta/users/vkaya/gnn/matrix/matrix/HeLa_10000.obs_exp_qt.npz \
    --batch_size 512 \
    --n_workers 12 \
    --results_save /cta/users/vkaya/gnn/diff-gnn/experiments/tag3/results_t0_t12.tsv
