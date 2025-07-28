#!/usr/bin/env bash
set -euo pipefail

# Step 0: Convert .dssp ➝ .ss
python AiCE_scripts/convert_dssp_to_ss.py

# Step 0.1: Copy .fa, .pdb from input to output
cp input/*.fa output/
cp input/*.pdb output/

# Step 0.2: Copy converted .ss files to output
cp output/single_result/*.ss output/

# Step 1: Run single mutation prediction
bash AiCE_scripts/01.post_mpnn_pipeline.sh AiCE_scripts ../output ../output ../output

# Step 2: LD
python AiCE_scripts/02.caculated_ld.py ../output ../ld_output

# Step 3: SCA
#bash AiCE_scripts/03.caculated_sca.sh   ./pySCA   ../output   ../sca_result

# Step 4: Combined mutation prediction
#bash AiCE_scripts/04.com_mut_prediction.sh   AiCE_scripts   ../output   2   ../com_result
