import nibabel as nib
from nilearn import input_data

Schaefer400 = nib.load('//mnt/nfs/lss/lss_kahwang_hpc/ROIs/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
Schaefer400_masker = input_data.NiftiLabelsMasker(Schaefer400)


PC =  Schaefer400_masker.fit_transform("/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/hubs/MGH_PC.nii")
stats = ["Context", "CxCO", "CxSH", "Task_Performed", "Resp", "CxCOxTP", "CxSHxTP"]

for stat in stats:
    beta = Schaefer400_masker.fit_transform("/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/RSA/searchlight/3dttest_20230607/GroupAnalysis_%s_stats.nii.gz" %stat)
    r = np.corrcoef(beta, PC)[0,1]
    print("%s : %s" %(stat, r))


# Context : 0.28086333553047227
# CxCO : 0.11394798320115816
# CxSH : 0.20014999658180274
# Task_Performed : -0.2251467315298735
# Resp : -0.20504152512343257
# CxCOxTP : -0.32132465675738003
# CxSHxTP : 0.07678691495025437
