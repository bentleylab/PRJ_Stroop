#!/bin/bash

# Move to scripts directory 
root_dir="/Volumes/hoycw_clust/PRJ_Stroop/"
cd $root_dir/scripts/

declare -a SBJs
SBJs=(`cat "${root_dir}scripts/SGE/SBJ_lists/MCC2019_poster_list.sbj"`)

# Run BHV scripts for these SBJs
for sbj in "${SBJs[@]}"; do
    cd ${root_dir}/data/${sbj}/05_recon/
    mkdir old_pipeline/
    mv *full* old_pipeline/
    mv *gROI* old_pipeline/
    mv *compiled* old_pipeline/
done


