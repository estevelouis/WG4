#!/bin/bash

languages_str="en zh"

read -a languages <<< "$languages_str"

for lang in "${languages[@]}"; do
    echo "Processing $lang" 
    cat <<EOT >sample.sh
#!/bin/bash
#SBATCH --job-name=sp_$lang
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=30:30:00
#SBATCH --mem=512GB
#SBATCH --partition=small
#SBATCH --account=project_462000447
#SBATCH --out=./slurmlog/%x_%j.out.log
#SBATCH --err=./slurmlog/%x_%j.err.log

source /flash/project_462000506/miniconda3/etc/profile.d/conda.sh 
conda activate eda_py3.12
rm -rf /dev/shm/*
hostname

cd /scratch/project_462000447/members/shaoxion/WG4/prepare_data
time python -u sample.py --lang_code $lang --sample_rate 0.1
EOT
    sbatch sample.sh
done
