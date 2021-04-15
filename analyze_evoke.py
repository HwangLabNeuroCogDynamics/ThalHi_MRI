# look at HRF for diff conditions
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from nibabel.processing import resample_from_to
from nilearn.image import resample_to_img
import nilearn
import scipy
import os
import glob
import statsmodels.api as sm
import statsmodels.formula.api as smf
from nilearn import masking
from nilearn import plotting
import matplotlib.pyplot as plt
from scipy.stats import zscore
sns.set_context("paper")

#setup variables
data_path = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve'
rois = '/data/backed_up/shared/ThalHi_MRI_2020/ROIs/'
# need thalamus mask.
# need morel mask

#list of subjects
subjects = [10001, 10002, 10003, 10004, 10005, 10006, 10007,
10008, 10009, 10010, 10013, 10014, 10016,
10017, 10018, 10019, 10020, 10021, 10024, 10025];
# 20 good subs

# masks
morel_mask = nib.load(rois+'Morel_2.5_mask.nii.gz')
morel_nuclei = nib.load(rois+'Morel_2.5.nii.gz')

### dump voxel-wise HRF time series
#variables
conditions = ['IDS', 'EDS', 'Stay']
ts={}

for condition in conditions:
    ts[condition] = np.zeros((len(subjects), 9, 2473)) #subject by time by voxel
    for ix, s in enumerate(subjects):
        fn = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve/sub-%s/sub-%s_%s_FIR_MNI.nii.gz' %(s, s, condition)
        fir_hrf = nib.load(fn)
        ts[condition][ix, :,:] = masking.apply_mask(fir_hrf, morel_mask)  #time by voxel


### look for main effects of condtions, map voxels
# use permutation test?
from mlxtend.evaluate import permutation_test
# or ttest
from scipy.stats import ttest_rel
from scipy.stats import ttest_1samp

# loop through each foxel
for time in np.arange(0, ts['Stay'].shape[1]-1):

    t_eds_v_ids = np.zeros(ts['Stay'].shape[2])
    t_eds_v_stay = np.zeros(ts['Stay'].shape[2])
    t_ids_v_stay = np.zeros(ts['Stay'].shape[2])
    p_eds_v_ids = np.zeros(ts['Stay'].shape[2])
    p_eds_v_stay = np.zeros(ts['Stay'].shape[2])
    p_ids_v_stay = np.zeros(ts['Stay'].shape[2])
    t_eds = np.zeros(ts['Stay'].shape[2])
    t_ids = np.zeros(ts['Stay'].shape[2])
    t_stay = np.zeros(ts['Stay'].shape[2])
    p_eds = np.zeros(ts['Stay'].shape[2])
    p_ids = np.zeros(ts['Stay'].shape[2])
    p_stay = np.zeros(ts['Stay'].shape[2])

    for vox in np.arange(0, ts['Stay'].shape[2]):
        b_stay = ts['Stay'][:,time,vox]
        b_ids = ts['IDS'][:,time,vox]
        b_eds = ts['EDS'][:,time,vox]

        t_eds_v_ids[vox], p_eds_v_ids[vox] = ttest_rel(b_eds, b_ids)
        t_eds_v_stay[vox], p_eds_v_stay[vox] = ttest_rel(b_eds, b_stay)
        t_ids_v_stay[vox], p_ids_v_stay[vox] = ttest_rel(b_ids, b_stay)

        t_eds[vox], p_eds[vox] = ttest_1samp(b_eds,0)
        t_ids[vox], p_ids[vox] = ttest_1samp(b_ids,0)
        t_stay[vox], p_stay[vox] = ttest_1samp(b_stay,0)


    ### look for switch effects, map voxels
    fn = 'images/t_eds_v_ids_%s' %time
    masking.unmask(t_eds_v_ids, morel_mask).to_filename(fn)
    fn = 'images/t_eds_v_stay_%s' %time
    masking.unmask(t_eds_v_stay, morel_mask).to_filename(fn)
    fn = 'images/t_ids_v_stay_%s' %time
    masking.unmask(t_ids_v_stay, morel_mask).to_filename(fn)

    ### look for main effects
    fn = 'images/t_eds_%s' %time
    masking.unmask(t_eds, morel_mask).to_filename(fn)
    fn = 'images/t_ids_%s' %time
    masking.unmask(t_ids, morel_mask).to_filename(fn)
    fn = 'images/t_stay_%s' %time
    masking.unmask(t_stay, morel_mask).to_filename(fn)



### Find cortical ROIs for EDS and b_ids

### Do FC to identify thalamic convergence zones with cortical EDS IDS ROIs
# if converge or if not converge?

### beta sereis
# show task depednent coupling

### plot  ts in overlapping zone















# end of line
