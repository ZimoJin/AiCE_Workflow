#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Full AiCE pipeline
# Assumes inputs are .pdb, .fa, and .dssp files in ./input
# Output structure: ./output/{single_result, ld_output, sca_result, com_result}
# -----------------------------------------------------------------------------

set -euo pipefail

echo "===== [STEP 0] Converting .dssp ➝ .ss ====="
python scripts/convert_dssp_to_ss.py

echo "===== Copying input files to output folder ====="
mkdir -p output
cp input/*.fa output/
cp input/*.pdb output/
cp output/single_result/*.ss output/

echo "===== [STEP 1] Single mutation prediction ====="
bash scripts/01.post_mpnn_pipeline.sh scripts output output output
cp -v output/single_result/* output/

echo "===== [STEP 2] Calculating LD ====="
mkdir -p ld_output
cd scripts
python 02.caculated_ld.py ../output ../ld_output
cd ..
mv -v ld_output output/
cp -v output/ld_output/* output/


echo "===== [STEP 3] Calculating SCA ====="
mkdir -p sca_result
bash scripts/03.caculated_sca.sh scripts/pySCA output sca_result
mv -v sca_result output/
cp -v output/sca_result/* output/

echo "===== [STEP 4] Combined mutation prediction ====="
mkdir -p com_result
bash scripts/04.com_mut_prediction.sh scripts output 2 com_result
mv -v com_result output/
cp -v output/com_result/* output/

echo "🎉 All steps completed successfully!"
