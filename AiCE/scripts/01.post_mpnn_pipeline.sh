#!/bin/bash

# ------------------------------------------------------------------------------
# This version skips MPNN & DSSP. It uses existing .fa and .ss files directly.
# Usage:
#   ./01.post_mpnn_pipeline.sh <scripts_dir> <fa_input_folder> <pdb_folder> [output_folder]
# ------------------------------------------------------------------------------

if [ $# -lt 3 ]; then
    echo "Usage: $0 <scripts_dir> <fa_input_folder> <pdb_folder> [output_folder]"
    echo "Example: $0 ../scripts ./input ../input ../output"
    exit 1
fi

scripts_dir="$1"
fa_input_dir="$2"
pdb_input_dir="$3"
output_dir="${4:-../output}"
single_dir="$output_dir/single_result"

mkdir -p "$single_dir"

shopt -s nullglob
seq_files=("$fa_input_dir"/*.fa)
shopt -u nullglob

if [ ${#seq_files[@]} -eq 0 ]; then
    echo "❌ No .fa files found in $fa_input_dir"
    exit 2
fi

echo "Found ${#seq_files[@]} .fa files. Proceeding..."

for seq_file in "${seq_files[@]}"; do
    file_prefix=$(basename "$seq_file" .fa)
    echo "-----------------------------------------"
    echo "[INFO] Processing: $file_prefix.fa"

    cp "$seq_file" "$single_dir/${file_prefix}.fa"

    tmp_file="$single_dir/${file_prefix}.tmp.fa"
    awk 'BEGIN{c=0}; /^>/{c++; if(c==1){print ">ref"} else {print ">" c-1}; next}; {print}' \
        "$single_dir/${file_prefix}.fa" > "$tmp_file"
    mv "$tmp_file" "$single_dir/${file_prefix}.fa"
    echo "  [Done] Reformat FASTA headers"

    ss_file="$single_dir/${file_prefix}.ss"
    if [ ! -f "$ss_file" ]; then
        echo "❌ Missing $ss_file, please ensure .ss files exist."
        continue
    fi

    echo "[Step B] Residue frequency..."
    python "$scripts_dir/count_residue_freq.py" \
        --msa "$single_dir/${file_prefix}.fa" \
        --count "$single_dir/${file_prefix}.freq" \
        --out "$single_dir/${file_prefix}.re.freq" || {
        echo "Error: Frequency calculation failed for $file_prefix"
        exit 3
    }

    echo "[Step C] Predicting mutations..."
    python "$scripts_dir/predicted_single_HF_mutations.py" \
        -freq "$single_dir/${file_prefix}.re.freq" \
        -dssp "$ss_file" \
        -comb "$single_dir/${file_prefix}.comb" \
        -mut "$single_dir/${file_prefix}.mut" \
        -model_path "$scripts_dir" || {
        echo "Error: Mutation prediction failed for $file_prefix"
        exit 4
    }

    echo "✅ Done: $file_prefix"
done

echo "-----------------------------------------"
echo "Pipeline finished successfully!"
ls -lh "$single_dir"
