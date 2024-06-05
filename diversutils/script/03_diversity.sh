#!/bin/bash
languages_str="af az bg ca cy de en es eu fi ga gu he hy is ja kk ko la lv ml mr mt nb nl pa ps ro si sl sq sw te tl tt ur vi ar be bn cs da el eo et fa fr gl hbs hi hu id it ka kn ky lt mk mn ms my ne nn pl pt ru sk so sv ta th tr uk uz zh"
read -a languages <<< "$languages_str"
languages=("zh")
for lang in "${languages[@]}"; do
    echo "Processing $lang" 
    cat <<EOT >cal_diversity.sh
#!/bin/bash
#SBATCH --job-name=div_$lang
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --mem=224GB
#SBATCH --partition=small
#SBATCH --account=project_462000447
#SBATCH --out=./slurmlog/%x_%j.out.log
#SBATCH --err=./slurmlog/%x_%j.err.log

rm -rf /dev/shm/*
hostname

# cd /scratch/project_462000447/members/shaoxion/WG4/diversutils
# make measurement_standard
cd /scratch/project_462000447/members/shaoxion/WG4/diversutils/bin
word2vec=/scratch/project_462000447/members/shaoxion/WG4/data/word2vec/bin
data_dir=/scratch/project_462000447/members/shaoxion/WG4/hpltv1.2_samples/$lang

txt_dir=/scratch/project_462000447/members/shaoxion/WG4/data/measurement_files
output_dir=/scratch/project_462000447/members/shaoxion/WG4/data/measurement_output_disparity/$lang
mkdir -p \$output_dir

if [ ! -d "\$txt_dir" ]; then
    echo "input_path is \$txt_dir/$lang.txt"
else
    for file in "\$data_dir"/*; do
        realpath "\$file"
    done > "\$txt_dir/$lang.txt"
    echo "input_path is \$txt_dir/$lang.txt"
fi

CMD_ARGS="--enable_pairwise=0 \\
    --enable_non_disparity_functions=0 \\
    --stirling_alpha=0.25 \\
    --num_row_threads=1 \\
    --num_matrix_threads=1 \\
    --row_generation_batch_size=1 \\
    --w2v_path=\$word2vec/hpltv1.2_${lang}_vec_100.bin \\
    --input_path=\$txt_dir/$lang.txt \\
    --output_path=\$output_dir/mesurement_output.tsv \\
    --force_timing_and_memory_to_output_path=1 \\
    --enable_document_count_recompute_step=1 \\
    --document_recompute_step_use_log10=1 \\
    --document_count_recompute_step_log10=0.1
"
echo \$CMD_ARGS

time {
./main_measurement \$CMD_ARGS
}


EOT
    sbatch cal_diversity.sh
done
