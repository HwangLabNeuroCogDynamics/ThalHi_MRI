#import matplotlib as plt
import nilearn
import nibabel as nib # requires some packages (listed below import list)
import numpy as np
import pandas as pd
import glob
import nilearn.masking
from nilearn import image, maskers, masking

##data path, here I'm being lazy just testing on one subject, but you can modify the code to loop through all subjects.
data = "/home/kahwang/argon/data/ThalHi/3dDeconvolve/sub-10003/" #note this is a symlink from my home to the hpc partition.

##residuals from GLM
func_file = nib.load(data + "sub-10003_FIRmodel_errts.nii.gz" )

##from the residuals, extract time-series from 400 ROIs.
Schaefer400 = nib.load('/home/kahwang/argon/ROIs/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
Schaefer400_masker = maskers.NiftiLabelsMasker(Schaefer400)
cortical_ts = Schaefer400_masker.fit_transform(func_file)
cortical_ts = cortical_ts[np.sum(cortical_ts,axis=1)!=0] #lots of censor?

##calculate FC matrix
def generate_correlation_mat(x, y):
	"""Correlate each n with each m.
	Parameters
	----------
	x : np.array
	  Shape Num ROI N X Timepoints.
	y : np.array
	  Shape ROI num M X Timepoints.
	Returns
	-------
	np.array
	  N X M array in which each element is a correlation coefficient.
	"""
	mu_x = x.mean(1)
	mu_y = y.mean(1)
	n = x.shape[1]
	if n != y.shape[1]:
		raise ValueError('x and y must ' +
						 'have the same number of timepoints.')
	s_x = x.std(1, ddof=n - 1)
	s_y = y.std(1, ddof=n - 1)
	cov = np.dot(x,
				 y.T) - n * np.dot(mu_x[:, np.newaxis],
								  mu_y[np.newaxis, :])
	return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])

fcmat = generate_correlation_mat(cortical_ts.T, cortical_ts.T) # full correlation


def cal_fcpc(fcmat, CI):
	''' clacuate PC with thalamocortical FC data'''
	thresholds = [86,87,88,89,90,91,92,93,94,95,96,97,98]
	pc_vectors = np.zeros((fcmat.shape[0], len(thresholds)))
	for it, t in enumerate(thresholds):
		temp_mat = fcmat.copy()
		temp_mat[temp_mat<np.percentile(temp_mat, t)] = 0 #threshold
		fc_sum = np.sum(temp_mat, axis=1)
		kis = np.zeros(np.shape(fc_sum))

		for ci in np.unique(CI):
			kis = kis + np.square(np.sum(temp_mat[:,np.where(CI==ci)[0]], axis=1) / fc_sum)

		pc_vectors[:, it] = 1-kis
	return pc_vectors, np.nanmean(pc_vectors, axis = 1)


## from FC matrx, calculate PC value of each ROI.
CI = np.loadtxt('/home/kahwang/argon/ROIs/Schaeffer400_7network_CI')

pc_alldensities, ave_pc = cal_fcpc(fcmat, CI) # will have to plot this in brain space to see if it makes sense.


## now look at RSA betas, which representation most strongly correlate with network hub? probably context.
rsa_betas= "/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/RSA/searchlight/Model_with_C_RF_TP_R/"
context_map = Schaefer400_masker.fit_transform(rsa_betas + "sub-10003_Context_beta_map.nii")
feature_map = Schaefer400_masker.fit_transform(rsa_betas + "sub-10003_Relevant_Feature_beta_map.nii")
task_map = Schaefer400_masker.fit_transform(rsa_betas + "sub-10003_Task_Performed_beta_map.nii")
resp_map = Schaefer400_masker.fit_transform(rsa_betas + "sub-10003_Resp_beta_map.nii")
#then you can correlate these betas with PC values, compare the r values across subjects, or do some kind of multi level regression
#np.corrcoef is your friend here.

## now look at *overlap* of RSA betas, are there ROIs that show strong responses for multiple representations?
# here I invented a PC like measure... 
rsa_betas = np.array([context_map, feature_map, task_map, resp_map])
rsa_betas = np.squeeze(rsa_betas).T #ROI by beta
beta_sum = np.sum(np.abs(rsa_betas), axis=1) #sum total beta weight for each ROI.
kis = np.zeros(np.shape(beta_sum)[0])
for c in np.arange(4):
    kis = kis + np.square(rsa_betas[:,c]/beta_sum) 
beta_pc = 1-kis # this is very similar to how PC for network hubs are calculated, except I didn't threshold. Not sure if this makes any sense at all.  # can make a whole brain plot to test.

