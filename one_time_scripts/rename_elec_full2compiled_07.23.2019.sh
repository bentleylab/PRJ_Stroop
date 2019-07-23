#!/bin/bash

# Move to scripts directory 
root_dir="/Volumes/hoycw_clust/PRJ_Stroop/"
cd $root_dir/scripts/

declare -a SBJs
SBJs=(`cat "${root_dir}scripts/SGE/SBJ_lists/${list_name}.sbj"`)

# Run BHV scripts for these SBJs
for sbj in "${SBJs[@]}"; do
    mv ${root_dir}/data/${sbj}/05_recon/${sbj}_elec_main_ft_pat_Dx_full.mat ${root_dir}/data/${sbj}/05_recon/${sbj}_elec_main_ft_pat_Dx_compiled.mat
done


