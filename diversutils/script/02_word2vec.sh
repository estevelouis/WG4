#!/bin/bash
languages_str="af az bg ca cy de es eu fi ga gu he hy is ja kk ko la lv ml mr mt nb nl pa ps ro si sl sq sw te tl tt ur vi ar be bn cs da el eo et fa fr gl hbs hi hu id it ka kn ky lt mk mn ms my ne nn pl pt ru sk so sv ta th tr uk uz"
read -a languages <<< "$languages_str"

for lang in "${languages[@]}"; do
    echo "Processing $lang" 
    cat <<EOT >train_word2vec.sh
#!/bin/bash
#SBATCH --job-name=wv_$lang
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=70:30:00
#SBATCH --mem=512GB
#SBATCH --partition=small
#SBATCH --account=project_462000447
#SBATCH --out=./slurmlog/%x_%j.out.log
#SBATCH --err=./slurmlog/%x_%j.err.log

rm -rf /dev/shm/*
hostname

hplt_data_dir=/scratch/project_462000447/members/shaoxion/WG4/hpltv1.2_samples
data_dir=/scratch/project_462000447/members/shaoxion/WG4/data

cd /scratch/project_462000447/members/shaoxion/WG4/diversutils
mkdir -p \$data_dir/word2vec/raw_txt
mkdir -p \$data_dir/word2vec/bin

raw_text=\$data_dir/word2vec/raw_txt/${lang}_raw.txt
if [ -e "\$raw_text" ]; then
    echo "Raw text exists."
else
    ./bin/jsonl_to_word2vec_format \$raw_text \$hplt_data_dir/$lang/*.jsonl
fi

dim=100
echo "The number of CPUs per task: \$SLURM_CPUS_PER_TASK"
./word2vec/bin/word2vec -train \$data_dir/word2vec/raw_txt/${lang}_raw.txt -output \$data_dir/word2vec/bin/hpltv1.2_${lang}_vec_\$dim.bin -cbow 1 -size \$dim -window 10 -negative 10 -hs 0 -threads \$SLURM_CPUS_PER_TASK -binary 1 -iter 3 -min-count 1
EOT
    sbatch train_word2vec.sh
done
