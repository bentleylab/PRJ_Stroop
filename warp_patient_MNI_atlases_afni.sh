cd /Volumes/hoycw_clust/PRJ_Stroop/data/atlases/Yeo/Yeo_JNeurophysiol11_MNI152

# Skull stripping can be done by @auto_tlrc, so skipping this step
# 3dSkullStrip -prefix single_subj_T1_1mm_noskull.nii -input single_subj_T1_1mm.nii 

# Use @auto_tlrc with default parameteres to see how it does
#@auto_tlrc -base single_subj_T1_1mm.nii -input FSL_MNI152_FreeSurferConformed_1mm.nii.gz -suffix _colin27.nii
# this failed by not being abel to open some of the files it generated... weird.

# Try with align_epi_anat.py
#   -cost lpa: recommended in the help
# align_epi_anat.py -dset1 FSL_MNI152_FreeSurferConformed_1mm.nii.gz \
#         -dset2 single_subj_T1_1mm.nii -dset1to2 \
#         -suffix _colin27 -cost lpa
# NOTE: This used only linear transformations, so the result is basically the same but resized a bit...

# Try 3dQwarp for non-linear alignment...
#   below call took 43m on my laptop
3dQwarp -allineate -blur 3 3 \
        -base single_subj_T1_1mm_noskull.nii \
        -source FSL_MNI152_FreeSurferConformed_1mm.nii.gz \
        -prefix FSL_MNI152_FreeSurferConformed_1mm_colin27_nl3-3.nii
# Now apply the warp from the anatomical alignment to the atlases
# use NN (Nearest Neighbor) interpolation to avoid smearing the atlas labels
3dNwarpApply -nwarp FSL_MNI152_FreeSurferConformed_1mm_colin27_nl3-3_WARP.nii \
        -source Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
        -prefix Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii \
        -ainterp NN
3dNwarpApply -nwarp FSL_MNI152_FreeSurferConformed_1mm_colin27_nl3-3_WARP.nii \
        -source Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
        -prefix Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii \
        -ainterp NN

