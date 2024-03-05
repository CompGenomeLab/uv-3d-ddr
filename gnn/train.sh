#!/bin/bash
#SBATCH --account=investor

#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:a100:1

#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --time=1-00:00:00
#SBATCH --signal=SIGHUP@90
#SBATCH --mem=60G

source /cta/users/vkaya/miniconda3/bin/activate gnn

echo "CUDA devices:" $CUDA_VISIBLE_DEVICES

dataset="/cta/users/vkaya/gnn/matrix/hic_4DNFIBM9QCFG_nMax51_nMin31_perc10"
srun --unbuffered python train.py \
    --n_workers 14 --n_devices 1 --n_nodes 1 \
    --dataset $dataset \
    --tag 3 \
    --env slurm \
    --n_epoch 20 \
    --edge_margin 0.15

