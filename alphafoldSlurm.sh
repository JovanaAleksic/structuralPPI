#!/bin/bash
#SBATCH --job-name=alphafold_run
#SBATCH --cpus-per-task=16         # max of 8 cores, AlphaFold has no benefit to use more
#SBATCH -w acnodeg01
#SBATCH --partition=gpu
#SBATCH --output=alphafold_%j.out

# Store the current directory
CURRENT_DIR=$(pwd)

# Load necessary modules
module load anaconda3
 
#activate myenv
source /acfs-home/joa4001/.bashrc
conda activate alphafold
module load alphafold/2.3.2

# Set AlphaFold src and data directory
export ALPHAFOLD_DATA_DIR=/data/apps/software/alphafold/2.3.2
export DATA_DIR=/data/apps/software/alphafold/2.3.2/db




python3 ${ALPHAFOLD_DATA_DIR}/run_alphafold.py \
  --fasta_paths=$CURRENT_DIR/mdm2_akt1.fasta \
  --output_dir=$CURRENT_DIR/output \
  --max_template_date=2022-01-01 \
  --data_dir=$DATA_DIR \
  --db_preset=full_dbs \
  --uniref30_database_path=$DATA_DIR/uniref30/UniRef30_2021_03 \
  --bfd_database_path=$DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
  --uniref90_database_path=$DATA_DIR/uniref90/uniref90.fasta \
  --mgnify_database_path=$DATA_DIR/mgnify/mgy_clusters_2022_05.fa \
  --uniprot_database_path=$DATA_DIR/uniprot/uniprot.fasta \
  --pdb_seqres_database_path=$DATA_DIR/pdb_seqres/pdb_seqres.txt \
  --use_gpu_relax=True \
  --template_mmcif_dir=$DATA_DIR/pdb_mmcif/mmcif_files \
  --obsolete_pdbs_path=$DATA_DIR/pdb_mmcif/obsolete.dat \
  --model_preset=multimer


