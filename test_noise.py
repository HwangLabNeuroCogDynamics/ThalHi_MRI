import nibabel as nib
from nilearn.image import concat_imgs, index_img, new_img_like
from nilearn import input_data

subjects=['10004','10005','10007','10008','10010',
'10012','10013','10014','10017','10018',
'10019','10020','10021','10022','10023',
'10025','10026','10027','10028','10031',
'10032','10033','10035','10036','10038',
'10040','10041','10042','10043','10044',
'10054','10058','10059','10060','10063',
'10064','10066','10069','10071','10072',
'10073','10074','10076','10080','10162',
'10173','10174','10175','10176','10179']

data_path = "/home/kahwang/bsh/ThalHi_MRI_2020/3dDeconvolve/"
Schaefer400 = nib.load('/data/backed_up/shared/ROIs/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
Schaefer400_masker = input_data.NiftiLabelsMasker(Schaefer400)

EDS_reliability = np.zeros(len(subjects))
IDS_reliability = np.zeros(len(subjects))
Stay_reliability = np.zeros(len(subjects))

for i, s in enumerate(subjects):
    f = nib.load(data_path+"sub-"+s+"/sub-"+s+ "_FIRmodel_MNI_stats_task_switch_r1_r4.nii.gz")
    data = f.get_fdata()
    data = np.squeeze(data)
    f = new_img_like(f,data)

    EDS_beta_ses1 = index_img(f,5) ## these are the betas from the FIR, not t stat
    IDS_beta_ses1 = index_img(f,8)
    Stay_beta_ses1 = index_img(f,11)

    f = nib.load(data_path+"sub-"+s+"/sub-"+s+ "_FIRmodel_MNI_stats_task_switch_r5_r8.nii.gz")
    data = f.get_fdata()
    data = np.squeeze(data)
    f = new_img_like(f,data)

    EDS_beta_ses2 = index_img(f,5)
    IDS_beta_ses2 = index_img(f,8)
    Stay_beta_ses2 = index_img(f,11)

    # split hlaf reliability.
    EDS_reliability[i] = np.corrcoef(Schaefer400_masker.fit_transform(EDS_beta_ses1), Schaefer400_masker.fit_transform(EDS_beta_ses2))[0,1]
    IDS_reliability[i] = np.corrcoef(Schaefer400_masker.fit_transform(IDS_beta_ses1), Schaefer400_masker.fit_transform(IDS_beta_ses2))[0,1]
    Stay_reliability[i] = np.corrcoef(Schaefer400_masker.fit_transform(Stay_beta_ses1), Schaefer400_masker.fit_transform(Stay_beta_ses2))[0,1]

