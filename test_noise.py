import numpy as np
import nibabel as nib
from nilearn.image import concat_imgs, index_img, new_img_like
from nilearn import input_data
import seaborn as sns
import matplotlib.pyplot as plt
plt.ion()

subjects=['10001', '10002', '10003', '10004', '10005',
'10008', '10010', '10012', '10013', '10014',
'10017', '10018', '10019', '10020', '10022',
'10023', '10024', '10025', '10027', '10028',
'10031', '10032', '10033', '10034', '10035',
'10036', '10037', '10038', '10039', '10040',
'10041', '10042', '10043', '10044', '10054',
'10057', '10058', '10059', '10060', '10063',
'10064', '10066', '10068', '10069', '10071',
'10072', '10073', '10074', '10076', '10077',
'10080', '10162', '10169', '10170', '10173',
'10174', '10175', '10176', '10179']

# Tesing calculating noise ceiling as true variance / total varaince, as described here:
# https://github.com/maedbhk/cerebellum_connectivity/blob/master/notebooks/2.1-Evaluation_noiseceilings.ipynb

# calculate the original test retest reliability and convert to noise ceiling sqrt(r)
ses1 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/Stay_59subs_run14_predicted_ctx_evrs_max.npy')
ses2 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/Stay_59subs_run58_predicted_ctx_evrs_max.npy')

R = np.zeros(59)
for s in np.arange(59):
    R[s] = np.sqrt(np.corrcoef(ses1[s,:],ses2[s,:])[0,1])

#now calculate it as portion of total variance:  np.sqrt(vU) / np.sqrt(vU + vE), vE is the error 
R_V = np.zeros(59)
for s in np.arange(59):
    # first concat data from two sessions
    Y = np.r_[ses1[s,:],ses2[s,:]]
    X = np.r_[ses2[s,:],ses1[s,:]]

    vU = np.sum((X-np.mean(X))**2)
    vE = np.sum((Y-X)**2)
    R_V[s] = np.sqrt(vU) / np.sqrt(vU+vE)


sns.pointplot(R_V, R)


#### OLD STUFF, IGNORE
# data_path = "/home/kahwang/bsh/ThalHi_MRI_2020/3dDeconvolve/"
# Schaefer400 = nib.load('/data/backed_up/shared/ROIs/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
# Schaefer400_masker = input_data.NiftiLabelsMasker(Schaefer400)

# EDS_reliability = np.zeros(len(subjects))
# IDS_reliability = np.zeros(len(subjects))
# Stay_reliability = np.zeros(len(subjects))

# for i, s in enumerate(subjects):
#     f = nib.load(data_path+"sub-"+s+"/sub-"+s+ "_FIRmodel_MNI_stats_task_switch_r1_r4.nii.gz")
#     data = f.get_fdata()
#     data = np.squeeze(data)
#     f = new_img_like(f,data)

#     EDS_beta_ses1 = index_img(f,5) ## these are the betas from the FIR, not t stat
#     IDS_beta_ses1 = index_img(f,8)
#     Stay_beta_ses1 = index_img(f,11)

#     f = nib.load(data_path+"sub-"+s+"/sub-"+s+ "_FIRmodel_MNI_stats_task_switch_r5_r8.nii.gz")
#     data = f.get_fdata()
#     data = np.squeeze(data)
#     f = new_img_like(f,data)

#     EDS_beta_ses2 = index_img(f,5)
#     IDS_beta_ses2 = index_img(f,8)
#     Stay_beta_ses2 = index_img(f,11)

#     # split hlaf reliability.
#     EDS_reliability[i] = np.corrcoef(Schaefer400_masker.fit_transform(EDS_beta_ses1), Schaefer400_masker.fit_transform(EDS_beta_ses2))[0,1]
#     IDS_reliability[i] = np.corrcoef(Schaefer400_masker.fit_transform(IDS_beta_ses1), Schaefer400_masker.fit_transform(IDS_beta_ses2))[0,1]
#     Stay_reliability[i] = np.corrcoef(Schaefer400_masker.fit_transform(Stay_beta_ses1), Schaefer400_masker.fit_transform(Stay_beta_ses2))[0,1]

