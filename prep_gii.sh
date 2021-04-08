


wb_command -volume-to-surface-mapping ~/bin/ThalHi_MRI/images/eds_stay.nii.gz \
/data/backed_up/shared/wb_files/HCP_S1200_GroupAvg_v1/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii \
~/bin/ThalHi_MRI/images/eds_stay_l.shape.gii -trilinear

wb_command -volume-to-surface-mapping ~/bin/ThalHi_MRI/images/ids_stay.nii.gz \
/data/backed_up/shared/wb_files/HCP_S1200_GroupAvg_v1/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii \
~/bin/ThalHi_MRI/images/ids_stay_l.shape.gii -trilinear

wb_command -volume-to-surface-mapping ~/bin/ThalHi_MRI/images/eds_ids.nii.gz \
/data/backed_up/shared/wb_files/HCP_S1200_GroupAvg_v1/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii \
~/bin/ThalHi_MRI/images/eds_ids_l.shape.gii -trilinear

wb_command -volume-to-surface-mapping ~/bin/ThalHi_MRI/images/eds_stay.nii.gz \
/data/backed_up/shared/wb_files/HCP_S1200_GroupAvg_v1/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii \
~/bin/ThalHi_MRI/images/eds_stay_r.shape.gii -trilinear

wb_command -volume-to-surface-mapping ~/bin/ThalHi_MRI/images/ids_stay.nii.gz \
/data/backed_up/shared/wb_files/HCP_S1200_GroupAvg_v1/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii \
~/bin/ThalHi_MRI/images/ids_stay_r.shape.gii -trilinear

wb_command -volume-to-surface-mapping ~/bin/ThalHi_MRI/images/eds_ids.nii.gz \
/data/backed_up/shared/wb_files/HCP_S1200_GroupAvg_v1/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii \
~/bin/ThalHi_MRI/images/eds_ids_r.shape.gii -trilinear
