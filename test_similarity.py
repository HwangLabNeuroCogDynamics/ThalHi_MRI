### given two vectors of voxel-wise patterns, what would be the best way to test their similarity?
import numpy as np
import nilearn
import nibabel as nib # requires some packages (listed below import list)
from scipy import stats
import pandas as pd
import re
import glob
import os
import sys
import nilearn.masking
import rsatoolbox
import sklearn
from scipy.spatial.distance import cosine
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr, zscore
import scipy.stats as stats
from mlxtend.evaluate import permutation_test
import datetime

### notes from Xitong
# thalamus voxel-wise evoked response:
# np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run14_thal_evrs.npy'.format(con))

# thalamocortical fc:
# fc_results=read_object('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/rsa_af/fc/{}/59subs_run14_run58_tsfc_max_clust_{}.p'.format(reg,clus))

# observed cortical pattern in rsa regions:
# ctx_results=read_object('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/rsa_af/evrs/{}/59subs_run14_run58_tent_evoked_response_clust_{}.p'.format(reg,clus))

# rsa masks:
# nib.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/rsa_af/masks/{}_pos_Clust_mask_{}.nii.gz'.format(reg,clus))

# rsa af normalized results:
# pd.read_csv('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/rsa_af/af/{}/59subs_ncaf_max.csv'.format(reg))

# reg and clus
# regressors={'Context':['0001','0002','0004','0005','4clusters'],
#             'Context_by_Color':['0002','0003','2clusters'],
#             'Context_by_Shape':['0003','0004','2clusters'],
#             'Task_Performed':['0002','1clusters'],
#             'Resp':['0001','0002','2clusters']}
# reg=regressors.keys()
# clus=regressors[reg]

def cal_reliability(x,y):
    Y = np.r_[x,y]
    X = np.r_[y,x]

    vU = np.sum((Y - Y.mean())**2)
    vE = np.sum((Y-X)**2)
    R_V = np.sqrt(vU) / np.sqrt(vU+vE)
    return R_V

def cal_EV(y,x): # encountering the same amplitdue scaling issue, looks like pearsonR is still the best option
    EV = 1-(np.sum((y-x)**2) / np.sum((y)**2))
    return EV

## load thalamus evoke response.
# note, looks like these are z scored?
EDS_tha_r14 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run14_thal_evrs.npy'.format("EDS"))
IDS_tha_r14 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run14_thal_evrs.npy'.format("IDS"))
Stay_tha_r14 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run14_thal_evrs.npy'.format("Stay"))
EDS_tha_r58 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run58_thal_evrs.npy'.format("EDS"))
IDS_tha_r58 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run58_thal_evrs.npy'.format("IDS"))
Stay_tha_r58 = np.load('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/400ROIs/evrs/{}_59subs_run58_thal_evrs.npy'.format("Stay"))

for clust in ['0001','0002','0004','0005']:
    ## load tha-ctx FC
    fc = read_object("/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/rsa_af/fc/Context/59subs_run14_run58_tsfc_max_clust_{}.p".format(clust))
    # ses1 [0][3] and ses2 [0][4]

    ## AF
    EDS_ses1 = np.zeros((59,fc[0][4].shape[1]))
    EDS_ses2 = np.zeros((59,fc[0][4].shape[1]))
    IDS_ses1 = np.zeros((59,fc[0][4].shape[1]))
    IDS_ses2 = np.zeros((59,fc[0][4].shape[1]))
    Stay_ses1 = np.zeros((59,fc[0][4].shape[1]))
    Stay_ses2 = np.zeros((59,fc[0][4].shape[1]))

    for s in np.arange(59):
        EDS_ses1[s,:] = np.dot(EDS_tha_r14[s,:], fc[s][3])    
        EDS_ses2[s,:] = np.dot(EDS_tha_r58[s,:], fc[s][4])
        IDS_ses1[s,:] = np.dot(IDS_tha_r14[s,:], fc[s][3])    
        IDS_ses2[s,:] = np.dot(IDS_tha_r58[s,:], fc[s][4])
        Stay_ses1[s,:] = np.dot(Stay_tha_r14[s,:], fc[s][3])    
        Stay_ses2[s,:] = np.dot(Stay_tha_r58[s,:], fc[s][4])

    ## obseved responses
    # looks like this is not normalized should we normalize this..? in the paper it seems to suggest it is zscored
    ctx_results = read_object('/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/Activity_Flow/noise_ceiling/rsa_af/evrs/Context/59subs_run14_run58_tent_evoked_response_clust_{}.p'.format(clust))

    ### calculate similarity
    EDS_ses1_r = np.zeros(59)
    EDS_ses2_r = np.zeros(59)
    IDS_ses1_r = np.zeros(59)
    IDS_ses2_r = np.zeros(59)
    Stay_ses1_r = np.zeros(59)
    Stay_ses2_r = np.zeros(59)

    ## calculate model reliability
    EDS_afR = np.zeros(59)
    IDS_afR = np.zeros(59)
    Stay_afR = np.zeros(59)
    for s in np.arange(59):
        EDS_afR[s] =  cal_reliability(EDS_ses2[s,:], EDS_ses1[s,:])
        IDS_afR[s] =  cal_reliability(IDS_ses2[s,:], IDS_ses1[s,:])
        Stay_afR[s] =  cal_reliability(Stay_ses2[s,:], Stay_ses1[s,:])
    
    ## calculate cortex reliability
    EDS_ctxR = np.zeros(59)
    IDS_ctxR = np.zeros(59)
    Stay_ctxR = np.zeros(59)
    for s in np.arange(59):
        EDS_ctxR[s] =  cal_reliability(zscore(ctx_results[s][1]['EDS']), zscore(ctx_results[s][2]['EDS']))
        IDS_ctxR[s] =  cal_reliability(zscore(ctx_results[s][1]['IDS']), zscore(ctx_results[s][2]['IDS']))
        Stay_ctxR[s] =  cal_reliability(zscore(ctx_results[s][1]['Stay']), zscore(ctx_results[s][2]['Stay']))

    EDS_R = np.sqrt(EDS_afR*EDS_ctxR) # this is noise ceiling accounting both model reliability and observed cortical reliability, gave very similar results
    IDS_R = np.sqrt(IDS_afR*IDS_ctxR)
    Stay_R = np.sqrt(Stay_afR*Stay_ctxR)

    for s in np.arange(59):
        # cross validation, compare observed from ses2 to predicted from ses1
        EDS_ses1_r[s] =  pearsonr(zscore(ctx_results[s][2]['EDS']), EDS_ses1[s,:])[0] / EDS_afR[s]  # only normalized by model reliability, thoguh using wholeR gave similar results
        EDS_ses2_r[s] =  pearsonr(zscore(ctx_results[s][1]['EDS']), EDS_ses2[s,:])[0] / EDS_afR[s]
        IDS_ses1_r[s] =  pearsonr(zscore(ctx_results[s][2]['IDS']), IDS_ses1[s,:])[0] / IDS_afR[s]
        IDS_ses2_r[s] =  pearsonr(zscore(ctx_results[s][1]['IDS']), IDS_ses2[s,:])[0] / IDS_afR[s]
        Stay_ses1_r[s] =  pearsonr(zscore(ctx_results[s][2]['Stay']), Stay_ses1[s,:])[0] / Stay_afR[s]
        Stay_ses2_r[s] =  pearsonr(zscore(ctx_results[s][1]['Stay']), Stay_ses2[s,:])[0] / Stay_afR[s]

    print("Cluster: ", clust)
    print("EDS mean: ", np.mean((EDS_ses2_r + EDS_ses1_r)/2))
    print("IDS mean: " ,np.mean((IDS_ses2_r + IDS_ses1_r)/2))
    print("Stay mean: ", np.mean((Stay_ses2_r + Stay_ses1_r)/2))
    
    print("EDS v IDS")
    x=(EDS_ses2_r + EDS_ses1_r)/2
    y=(IDS_ses2_r + IDS_ses1_r)/2
    print(stats.ttest_rel(x,y))
    print(stats.wilcoxon(x,y))
    print(permutation_test(x,y, method='approximate', num_rounds=1000, paired = True)) ## note, the p value flcutautes ..., but clust2 seems to have a reliable effect

    print("IDS v Stay")
    x=(Stay_ses2_r + Stay_ses1_r)/2
    y=(IDS_ses2_r + IDS_ses1_r)/2
    print(stats.ttest_rel(y,x))
    print(stats.wilcoxon(y,x))
    print(permutation_test(y,x, method='approximate', num_rounds=1000, paired = True)) 

    print("EDS v Stay")
    x=(Stay_ses2_r + Stay_ses1_r)/2
    y=(EDS_ses2_r + EDS_ses1_r)/2
    print(stats.ttest_rel(y,x))
    print(stats.wilcoxon(y,x))
    print(permutation_test(y,x, method='approximate', num_rounds=1000, paired = True)) 

