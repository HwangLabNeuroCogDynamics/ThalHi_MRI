# correlate thalamocritcal FC and cortical evoke patterns
import os
#limit numpy threads
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

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
from nilearn.input_data import NiftiLabelsMasker
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
sns.set_context("paper")
import multiprocessing
import pickle

def generate_correlation_mat(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

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


def pcorr_subcortico_cortical_connectivity(subcortical_ts, cortical_ts):
	''' function to do partial correlation bewteen subcortical and cortical ROI timeseries.
	Cortical signals (not subcortical) will be removed from subcortical and cortical ROIs,
	and then pairwise correlation will be calculated bewteen subcortical and cortical ROIs
	(but not between subcortico-subcortical or cortico-cortical ROIs).
	This partial correlation/regression approach is for cleaning subcortico-cortical
	conectivity, which seems to be heavily influenced by a global noise.
	usage: pcorr_mat = pcorr_subcortico-cortical(subcortical_ts, cortical_ts)
	----
	Parameters
	----
	subcortical_ts: txt file of timeseries data from subcortical ROIs/voxels, each row is an ROI
	cortical_ts: txt file of timeseries data from cortical ROIs, each row is an ROI
	pcorr_mat: output partial correlation matrix
	'''
	from scipy import stats, linalg
	from sklearn.decomposition import PCA

	# # transpose so that column is ROI, this is because output from 3dNetcorr is row-based.
	# subcortical_ts = subcortical_ts.T
	# cortical_ts = cortical_ts.T
	cortical_ts[np.isnan(cortical_ts)]=0
	subcortical_ts[np.isnan(subcortical_ts)]=0

	# check length of data
	assert cortical_ts.shape[0] == subcortical_ts.shape[0]
	num_vol = cortical_ts.shape[0]

	#first check that the dimension is appropriate
	num_cort = cortical_ts.shape[1]
	num_subcor = subcortical_ts.shape[1]
	num_total = num_cort + num_subcor

	#maximum number of regressors that we can use
	max_num_components = int(num_vol/20)
	if max_num_components > num_cort:
		max_num_components = num_cort-1

	pcorr_mat = np.zeros((num_total, num_total), dtype=np.float)

	for j in range(num_cort):
		k = np.ones(num_cort, dtype=np.bool)
		k[j] = False

		#use PCA to reduce cortical data dimensionality
		pca = PCA(n_components=max_num_components)
		pca.fit(cortical_ts[:,k])
		reduced_cortical_ts = pca.fit_transform(cortical_ts[:,k])

		#print("Amount of varaince explanined after PCA: %s" %np.sum(pca.explained_variance_ratio_))

		# fit cortical signal to cortical ROI TS, get betas
		beta_cortical = linalg.lstsq(reduced_cortical_ts, cortical_ts[:,j])[0]

		#get residuals
		res_cortical = cortical_ts[:, j] - reduced_cortical_ts.dot(beta_cortical)

		for i in range(num_subcor):
			# fit cortical signal to subcortical ROI TS, get betas
			beta_subcortical = linalg.lstsq(reduced_cortical_ts, subcortical_ts[:,i])[0]

			#get residuals
			res_subcortical = subcortical_ts[:, i] - reduced_cortical_ts.dot(beta_subcortical)

			#partial regression
			pcorr_mat[i+num_cort, j] = stats.pearsonr(res_cortical, res_subcortical)[0] #linalg.lstsq(res_subcortical.reshape(-1,1), res_cortical.reshape(-1,1))[0]
			pcorr_mat[j,i+num_cort ] = pcorr_mat[i+num_cort, j]

	return pcorr_mat

def save_object(obj, filename):
	''' Simple function to write out objects into a pickle file
	usage: save_object(obj, filename)
	'''
	with open(filename, 'wb') as output:
		pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
	#M = pickle.load(open(f, "rb"))


def read_object(filename):
	''' short hand for reading object because I can never remember pickle syntax'''
	o = pickle.load(open(filename, "rb"))
	return o

conditions = ['EDS', 'IDS','Stay']
subjects = ['sub-10023',
 'sub-10020',
 'sub-10014',
 'sub-10019',
 'sub-10022',
 'sub-10008',
 'sub-10026',
 'sub-10024',
 'sub-10025',
 'sub-10002',
 'sub-10005',
 'sub-10017',
 'sub-10012',
 'sub-10018',
 'sub-10007',
 'sub-10013',
 'sub-10010',
 'sub-10004',
 'sub-10003',
 'sub-10021']





####################################################
## function to do activity flow, thalamux evoke multiply by thalamocortical FC
####################################################
rois = '/data/backed_up/shared/ThalHi_MRI_2020/ROIs/'
thalamus_mask = nib.load(rois+'Morel_2.5_mask.nii.gz')
cortex_mask = nib.load(rois+'Schaefer400_2.5.nii.gz')
cortex_masker = NiftiLabelsMasker(labels_img=cortex_mask, standardize=False)
striatum_mask = nib.load(rois+'Striatum_2.5_mask.nii.gz')

#global rois, thalamus_mask, cortex_mask, cortex_masker, striatum_mask

def run_fc_evoke_corr(inputs):
    data_path = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve'
    rois = '/data/backed_up/shared/ThalHi_MRI_2020/ROIs/'

    # three elements expected in input
    s = inputs[0] # subject name
    subcortical_mask = inputs[1]  # subocortical mask
    cortex_masker = inputs[2] # cortex masker

    #thalamus_mask = nib.load(rois+'Morel_2.5_mask.nii.gz')
    #cortex_mask = nib.load(rois+'Schaefer400_2.5.nii.gz')
    #cortex_masker = NiftiLabelsMasker(labels_img=cortex_mask, standardize=False)

    subcortical_mask_size = np.sum(subcortical_mask.get_fdata()>0)
    roi_size = len(np.unique(cortex_masker.labels_img.get_fdata()))-1

    fcmat = np.zeros((subcortical_mask_size,roi_size))
    conditions = ['IDS', 'EDS', 'Stay']
    subcortical_evoke = {}
    ctx_evoke = {}
    for condition in conditions:
        subcortical_evoke[condition] = np.zeros((9, subcortical_mask_size)) #subject by time by voxel
        ctx_evoke[condition] = np.zeros((9, roi_size)) #subject by time by cortical ROI

    # FC
    fn = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve/%s/%s_FIRmodel_errts.nii.gz' %(s,s)
    functional_data = nib.load(fn)
    cortex_ts = cortex_masker.fit_transform(functional_data)
    subcortical_ts = masking.apply_mask(functional_data, subcortical_mask)

    # remove censored timepoints
    mts = np.mean(cortex_ts, axis = 1)
    if any(mts==0):
        del_i = np.where(mts==0)
        cortex_ts = np.delete(cortex_ts, del_i, axis = 0)
        subcortical_ts = np.delete(subcortical_ts, del_i, axis = 0)

    # now try principal compoment regression
    ts_len = cortex_ts.shape[0]
    n_comps = np.amin([ts_len, roi_size, subcortical_mask_size]) // 20
    pca = PCA(n_comps)
    reduced_mat = pca.fit_transform(subcortical_ts) # Time X components
    components = pca.components_
    regrmodel = LinearRegression()
    reg = regrmodel.fit(reduced_mat, cortex_ts) #cortex ts also time by ROI
    #project regression betas from component
    fcmat[:, :] = pca.inverse_transform(reg.coef_).T #reshape to cortex

    # correlation and partial correlation
    #fcmat[:, :] = generate_correlation_mat(thalamus_ts.T, cortex_ts.T) #th by ctx
    #fcmat[:, :] = pcorr_subcortico_cortical_connectivity(thalamus_ts, cortex_ts)[400:, 0:400]

    #Extract tha and cortical evoke
    for condition in conditions:
        fn = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve/%s/%s_%s_FIR_MNI.nii.gz' %(s, s, condition)
        fir_hrf = nib.load(fn)
        subcortical_evoke[condition][:,:] = masking.apply_mask(fir_hrf, subcortical_mask)  #time by voxel
        ctx_evoke[condition][:,:] = cortex_masker.fit_transform(fir_hrf) #time by cortical ROI

    return s, fcmat, subcortical_evoke, ctx_evoke


####################################################
## thalamux evoke multiply by thalamocortical FC
####################################################

## parallel
pool = multiprocessing.Pool(10)
results = pool.map(run_fc_evoke_corr, zip(subjects, [thalamus_mask]*len(subjects), [cortex_masker]*len(subjects)))
pool.close()
pool.join()

##### unpack results
print("correlation between observed and predicted cortical evoked pattern")
print(" ")
conditions = ['EDS', 'IDS','Stay']
predicted_ctx = np.zeros((3, results[0][1].shape[1]))
observed_ctx = np.zeros((3, results[0][1].shape[1]))
observed_tha = np.zeros((3, results[0][1].shape[0]))
corr = np.zeros((3, len(subjects)))
for ic, cond in enumerate(conditions):
    for ix , res in enumerate(results):
        fc = res[1]
        tha_b = np.mean(results[ix][2][cond][:,:],axis=0)
        ctx_b = np.mean(results[ix][3][cond][:,:],axis=0)
        # correlate predicted cortical evoke vs observed
        predicted_ctx[ic, :] = zscore(np.dot(tha_b, fc)) #
        observed_ctx[ic, :] = zscore(ctx_b)
        observed_tha[ic, :] = tha_b
        corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

    print(cond)
    print(corr[ic, :])
    print(np.mean(corr[ic, :]))

### some plotting
#sns.set_context('talk')
#sns.regplot(x=predicted_ctx[1,:], y=observed_ctx[1,:])
#plt.show()

# kde plot
df = pd.DataFrame()
for ic, cond in enumerate(conditions):
    tdf = pd.DataFrame()
    tdf['Correlation'] = corr[ic, :]
    tdf['Condition'] = cond
    df= pd.concat([df, tdf])

# condtions
sns.set_context('talk')
fig1=sns.pointplot(x="Condition", y="Correlation", hue="Condition",
			  data=df, dodge=False, join=False)

fig1=sns.stripplot(x="Condition", y="Correlation", hue="Condition",
			  data=df, dodge=False, alpha=.5)
fig1.legend_.remove()
fig1.set_ylim([0.3, 0.9])
plt.tight_layout()
plt.show()

####################################################
# set thalamus evoke input to uniform
####################################################
print(" ")
print(" ")
print("correlation between observed and predicted cortical evoked pattern after setting thalamus to unifomred activity")
conditions = ['EDS', 'IDS','Stay']
predicted_ctx = np.zeros((3, results[0][1].shape[1]))
observed_ctx = np.zeros((3, results[0][1].shape[1]))
observed_tha = np.zeros((3, results[0][1].shape[0]))
corr = np.zeros((3, len(subjects)))
for ic, cond in enumerate(conditions):
    for ix , res in enumerate(results):
        fc = res[1]
        tha_b = np.mean(results[ix][2][cond][:,:],axis=0)
        tha_b = np.ones(tha_b.shape) #uniform
        ctx_b = np.mean(results[ix][3][cond][:,:],axis=0)
        # correlate predicted cortical evoke vs observed
        predicted_ctx[ic, :] = zscore(np.dot(tha_b, fc))
        observed_ctx[ic, :] = zscore(ctx_b)
        observed_tha[ic, :] = tha_b
        corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

    print(cond)
    print(corr[ic, :])
    print(np.mean(corr[ic, :]))

### some plotting
# sns.set_context('talk')
# sns.regplot(x=predicted_ctx[1,:], y=observed_ctx[1,:])
# plt.show()

# kde plot of nul
df = pd.DataFrame()
for ic, cond in enumerate(conditions):
    tdf = pd.DataFrame()
    tdf['Correlation'] = corr[ic, :]
    tdf['Condition'] = cond
    df= pd.concat([df, tdf])

sns.kdeplot(x = 'Correlation', data=df)
plt.tight_layout()
plt.show()


####################################################
## permutation thalamus indicex to get null distribution
####################################################
print(" ")
print(" ")
print("creating an empircal null distribution by randomly swapping the thalamus evoke vector")
print(" ")

conditions = ['EDS', 'IDS','Stay']
predicted_ctx = np.zeros((3, results[0][1].shape[1]))
observed_ctx = np.zeros((3, results[0][1].shape[1]))
observed_tha = np.zeros((3, results[0][1].shape[0]))
corr = np.zeros((3, len(subjects)))
for ic, cond in enumerate(conditions):
    for ix , res in enumerate(results):
        fc = res[1]
        tha_b = np.random.permutation(np.mean(results[ix][2][cond][:,:],axis=0))
        ctx_b = np.mean(results[ix][3][cond][:,:],axis=0)
        # correlate predicted cortical evoke vs observed
        predicted_ctx[ic, :] = zscore(np.dot(tha_b, fc))
        observed_ctx[ic, :] = zscore(ctx_b)
        observed_tha[ic, :] = tha_b
        corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

    print(cond)
    print(corr[ic, :])
    print(np.mean(corr[ic, :]))


df = pd.DataFrame()
for ic, cond in enumerate(conditions):
    tdf = pd.DataFrame()
    tdf['Correlation'] = corr[ic, :]
    tdf['Condition'] = cond
    df= pd.concat([df, tdf])

sns.kdeplot(x = 'Correlation', data=df)
plt.tight_layout()
plt.show()


print(" ")
print(" ")
print("creating an empircal null distribution by randomly swapping the FC matrix")
print(" ")

conditions = ['EDS', 'IDS','Stay']
predicted_ctx = np.zeros((3, results[0][1].shape[1]))
observed_ctx = np.zeros((3, results[0][1].shape[1]))
observed_tha = np.zeros((3, results[0][1].shape[0]))
corr = np.zeros((3, len(subjects)))
for ic, cond in enumerate(conditions):
    for ix , res in enumerate(results):
        fc = np.random.permutation(res[1]) ### will only randomly shuffle the first axis, which is the thalamus vector
        tha_b = np.mean(results[ix][2][cond][:,:],axis=0)
        # correlate predicted cortical evoke vs observed
        predicted_ctx[ic, :] = zscore(np.dot(tha_b, fc))
        observed_ctx[ic, :] = zscore(ctx_b)
        observed_tha[ic, :] = tha_b
        corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

    print(cond)
    print(corr[ic, :])
    print(np.mean(corr[ic, :]))


df = pd.DataFrame()
for ic, cond in enumerate(conditions):
    tdf = pd.DataFrame()
    tdf['Correlation'] = corr[ic, :]
    tdf['Condition'] = cond
    df= pd.concat([df, tdf])

sns.kdeplot(x = 'Correlation', data=df)
plt.tight_layout()
plt.show()

####################################################
## random lesion the thalamus. 480 (20%) voxels at a time. Will need to keep track of which voxels.
####################################################
# lesion_idx = np.zeros((10000,240))
# corr = np.zeros((3, len(subjects), 10000))
# for n in np.arange(0,10000):
#     lesion_idx[n,:] = np.random.randint(2472, size=240)
#     for ic, cond in enumerate(conditions):
#         for ix , res in enumerate(results):
#             fc = res[1] ### will only randomly shuffle the first axis, which is the thalamus vector
#             tha_b = np.mean(results[ix][2][cond][:,:],axis=0)
#             tha_b[lesion_idx[0,:].astype('int')]=0
#             ctx_b = np.mean(results[ix][3][cond][:,:],axis=0)
#             # correlate predicted cortical evoke vs observed
#             corr[ic, ix, n] = np.corrcoef(zscore(np.dot(zscore(tha_b), fc)), zscore(ctx_b))[0,1]
#
global conditions
global results
lesion_size = 480
global lesion_size

def sim_lesion(idx):
    corr = np.zeros((3, len(subjects)))
    lesion_idx = np.random.permutation(np.arange(0,2473))[0:lesion_size]             #np.arange(idx, idx+1000)
    for ic, cond in enumerate(conditions):
        for ix , res in enumerate(results):
            fc = res[1]
            tha_b = np.mean(results[ix][2][cond][:,:],axis=0)
            tha_b[lesion_idx.astype('int')]=0
            ctx_b = np.mean(results[ix][3][cond][:,:],axis=0)
            # correlate predicted cortical evoke vs observed
            corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]
    return corr, lesion_idx

# run simulation
pool = multiprocessing.Pool(os.cpu_count() // 3 )
idxs = np.arange(2473*3)
lesion_results = pool.map(sim_lesion, idxs)
pool.close()
pool.join()


## unpack simulation
corr = np.zeros((3, len(subjects), len(idxs)))
lesion_idx = np.zeros((len(idxs),lesion_size))
for ix , res in enumerate(lesion_results):
    corr[:,:,ix] = res[0]
    lesion_idx[ix,:] = res[1]

for ic, cond in enumerate(conditions):
    print(cond)
    print(np.percentile(np.mean(corr[ic, :,:], axis=0), 99))
    print(np.percentile(np.mean(corr[ic, :,:], axis=0), 1))


####################################################
#### Do the same thing for BG

## parallel

pool = multiprocessing.Pool(10)
striatum_results = pool.map(run_fc_evoke_corr, zip(subjects, [striatum_mask] * len(subjects), [cortex_masker]* len(subjects)))
pool.close()
pool.join()


##### unpack results
print("correlation between observed and predicted cortical evoked pattern, for striatum")
print(" ")
conditions = ['EDS', 'IDS','Stay']
predicted_ctx = np.zeros((3, striatum_results[0][1].shape[1]))
observed_ctx = np.zeros((3, striatum_results[0][1].shape[1]))
observed_tha = np.zeros((3, striatum_results[0][1].shape[0]))
corr = np.zeros((3, len(subjects)))
for ic, cond in enumerate(conditions):
    for ix , res in enumerate(striatum_results):
        fc = res[1]
        tha_b = np.mean(striatum_results[ix][2][cond][:,:],axis=0)
        ctx_b = np.mean(striatum_results[ix][3][cond][:,:],axis=0)
        # correlate predicted cortical evoke vs observed
        predicted_ctx[ic, :] = zscore(np.dot(tha_b, fc))
        observed_ctx[ic, :] = zscore(ctx_b)
        observed_tha[ic, :] = tha_b
        corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

    print(cond)
    print(corr[ic, :])
    print(np.mean(corr[ic, :]))

df = pd.DataFrame()
for ic, cond in enumerate(conditions):
    tdf = pd.DataFrame()
    tdf['Correlation'] = corr[ic, :]
    tdf['Condition'] = cond
    df= pd.concat([df, tdf])
sns.set_context('talk')
fig1 = sns.pointplot(x="Condition", y="Correlation", hue="Condition",
			  data=df, dodge=False, join=False)

fig1 = sns.stripplot(x="Condition", y="Correlation", hue="Condition",
			  data=df, dodge=False, alpha=.5)
fig1.legend_.remove()
fig1.set_ylim([0.3, 0.9])
plt.tight_layout()
plt.show()

########################################################################################################
#### Use resting-state FC as the FC matrix for activity flow
########################################################################################################
MGH_fn = '/home/kahwang/bsh/MGH/MGH/*/MNINonLinear/rfMRI_REST_ncsreg.nii.gz'
MGH_files = glob.glob(MGH_fn)

def rsfc(input):
    file = input[0]
    subcortical_mask = input[1]
    cortex_masker = input[2]
    functional_data = nib.load(file)

    functional_data = nilearn.image.resample_to_img(functional_data, subcortical_mask)
    cortex_ts = cortex_masker.fit_transform(functional_data)
    subcortical_ts = masking.apply_mask(functional_data, subcortical_mask)

    # remove censored timepoints
    mts = np.mean(cortex_ts, axis = 1)
    if any(mts==0):
        del_i = np.where(mts==0)
        cortex_ts = np.delete(cortex_ts, del_i, axis = 0)
        subcortical_ts = np.delete(subcortical_ts, del_i, axis = 0)

    # now try principal compoment regression
    pca = PCA(239)
    reduced_mat = pca.fit_transform(subcortical_ts) # Time X components
    components = pca.components_
    regrmodel = LinearRegression()
    reg = regrmodel.fit(reduced_mat, cortex_ts) #cortex ts also time by ROI
    #project regression betas from component
    fcmat = pca.inverse_transform(reg.coef_).T

    return fcmat

pool = multiprocessing.Pool(10)
rsfc_results = pool.map(rsfc, zip(MGH_files, [thalamus_mask] * len(MGH_files), [cortex_masker] * len(MGH_files)))
pool.close()
pool.join()

rsfc_corrmat = np.zeros((len(MGH_files),2473,400))
for ir, res in enumerate(rsfc_results):
    rsfc_corrmat[ir, :,:] = res
rsfc_avemat = np.mean(rsfc_corrmat, axis=0)

# resting-state thalamocortical FC activityflow
print("correlation between observed and predicted cortical evoked pattern, using rsfc")
print(" ")
conditions = ['IDS', 'EDS', 'Stay']
for cond in conditions:
    corr = np.zeros(len(subjects))
    for ix , res in enumerate(results):
        fc = rsfc_avemat
        tha_b = np.mean(results[ix][2][cond][:,:],axis=0)
        ctx_b = np.mean(results[ix][3][cond][:,:],axis=0)
        # correlate predicted cortical evoke vs observed
        corr[ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

    print(cond)
    print(corr)
    print(np.mean(corr))

############################################################################################################
### restrict analysis in ROIs
############################################################################################################

rois = '/data/backed_up/shared/ThalHi_MRI_2020/ROIs/'
thalamus_mask = nib.load(rois+'Morel_2.5_mask.nii.gz')
striatum_mask = nib.load(rois+'Striatum_2.5_mask.nii.gz')
lpfc_mask = nib.load(rois+'lpfc.nii.gz')
rlpfc_mask = nib.load(rois+'rlplfc.nii.gz')
preCS_mask = nib.load(rois+'preCS.nii.gz')


def run_roi_fc_evoke_corr(inputs):
    data_path = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve'
    rois = '/data/backed_up/shared/ThalHi_MRI_2020/ROIs/'

    # three elements expected in input
    s = inputs[0] # subject name
    subcortical_mask = inputs[1]  # subocortical mask
    cortex_mask = inputs[2] # cortex masker

    #thalamus_mask = nib.load(rois+'Morel_2.5_mask.nii.gz')
    #cortex_mask = nib.load(rois+'Schaefer400_2.5.nii.gz')
    #cortex_masker = NiftiLabelsMasker(labels_img=cortex_mask, standardize=False)

    subcortical_mask_size = np.sum(subcortical_mask.get_fdata()>0)
    roi_size = np.sum(cortex_mask.get_fdata()>0)

    fcmat = np.zeros((subcortical_mask_size,roi_size))
    conditions = ['IDS', 'EDS', 'Stay']

    subcortical_evoke = {}
    ctx_evoke = {}
    for condition in conditions:
        subcortical_evoke[condition] = np.zeros((9, subcortical_mask_size)) #subject by time by voxel
        ctx_evoke[condition] = np.zeros((9, roi_size)) #subject by time by cortical ROI

    # FC
    fn = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve/%s/%s_FIRmodel_errts.nii.gz' %(s,s)
    functional_data = nib.load(fn)
    cortex_ts = masking.apply_mask(functional_data, cortex_mask)
    subcortical_ts = masking.apply_mask(functional_data, subcortical_mask)

    # remove censored timepoints
    mts = np.mean(cortex_ts, axis = 1)
    if any(mts==0):
        del_i = np.where(mts==0)
        cortex_ts = np.delete(cortex_ts, del_i, axis = 0)
        subcortical_ts = np.delete(subcortical_ts, del_i, axis = 0)

    # now try principal compoment regression
    ts_len = cortex_ts.shape[0]
    n_comps = np.amin([ts_len, roi_size, subcortical_mask_size]) // 20
    pca = PCA(n_comps)
    reduced_mat = pca.fit_transform(subcortical_ts) # Time X components
    components = pca.components_
    regrmodel = LinearRegression()
    reg = regrmodel.fit(reduced_mat, cortex_ts) #cortex ts also time by ROI
    #project regression betas from component
    fcmat[:, :] = pca.inverse_transform(reg.coef_).T #reshape to cortex

    # correlation and partial correlation
    #fcmat[:, :] = generate_correlation_mat(thalamus_ts.T, cortex_ts.T) #th by ctx
    #fcmat[:, :] = pcorr_subcortico_cortical_connectivity(thalamus_ts, cortex_ts)[400:, 0:400]

    #Extract tha and cortical evoke
    for condition in conditions:
        fn = '/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve/%s/%s_%s_FIR_MNI.nii.gz' %(s, s, condition)
        fir_hrf = nib.load(fn)
        subcortical_evoke[condition][:,:] = masking.apply_mask(fir_hrf, subcortical_mask)  #time by voxel
        ctx_evoke[condition][:,:] = masking.apply_mask(fir_hrf, cortex_mask)  #time by cortical ROI

    return s, fcmat, subcortical_evoke, ctx_evoke

## parallel
pool = multiprocessing.Pool(10)
lpfc_results = pool.map(run_roi_fc_evoke_corr, zip(subjects, [thalamus_mask]*len(subjects), [lpfc_mask]*len(subjects)))
pool.close()
pool.join()
save_object(lpfc_results, 'data/lpfc_results')

pool = multiprocessing.Pool(10)
rlpfc_results = pool.map(run_roi_fc_evoke_corr, zip(subjects, [thalamus_mask]*len(subjects), [rlpfc_mask]*len(subjects)))
pool.close()
pool.join()
save_object(rlpfc_results, 'data/rlpfc_results')

pool = multiprocessing.Pool(10)
preCS_results = pool.map(run_roi_fc_evoke_corr, zip(subjects, [thalamus_mask]*len(subjects), [preCS_mask]*len(subjects)))
pool.close()
pool.join()
save_object(preCS_results, 'data/preCS_results')

##### unpack results
print("correlation between observed and predicted cortical evoked pattern in LPFC")
print(" ")

def run_acfl(input_results):
    conditions = ['EDS', 'IDS','Stay']
    predicted_ctx = np.zeros((3, input_results[0][1].shape[1]))
    observed_ctx = np.zeros((3, input_results[0][1].shape[1]))
    observed_tha = np.zeros((3, input_results[0][1].shape[0]))
    corr = np.zeros((3, len(subjects)))
    fc_sum = np.zeros((20, input_results[0][1].shape[0], input_results[0][1].shape[1] ))
    for ic, cond in enumerate(conditions):
        for ix , res in enumerate(input_results):
            fc = res[1]
            fc_sum[ix, :,:] = fc = res[1]
            tha_b = np.mean(input_results[ix][2][cond][:,:],axis=0)
            ctx_b = np.mean(input_results[ix][3][cond][:,:],axis=0)
            # correlate predicted cortical evoke vs observed
            predicted_ctx[ic, :] = zscore(np.dot(tha_b, fc))
            observed_ctx[ic, :] = zscore(ctx_b)
            observed_tha[ic, :] = tha_b
            corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]

        print(cond)
        print(corr[ic, :])
        print(np.mean(corr[ic, :]))

    return corr, fc_sum

#run_acfl(rlpfc_results)
#run_acfl(preCS_results)
rlpfc_results = read_object('data/lpfc_results')
corr, fc_sum = run_acfl(rlpfc_results)
# condtions
df = pd.DataFrame()
for ic, cond in enumerate(conditions):
    tdf = pd.DataFrame()
    tdf['Correlation'] = abs(corr[ic, :])
    tdf['Condition'] = cond
    df= pd.concat([df, tdf])

sns.set_context('talk')
fig1=sns.pointplot(x="Condition", y="Correlation", hue="Condition",
			  data=df, dodge=False, join=False)

fig1=sns.stripplot(x="Condition", y="Correlation", hue="Condition",
			  data=df, dodge=False, alpha=.2)
fig1.legend_.remove()
fig1.set_ylim([0, 0.9])
plt.tight_layout()
plt.show()

from scipy.stats import ttest_rel
ttest_rel(corr[0,:], corr[1,:])
ttest_rel(corr[0,:], corr[2,:])
ttest_rel(corr[0,:], corr[2,:])






global conditions
global rlpfc_results
lesion_size = 480
global lesion_size

def sim_lesion(idx):
    corr = np.zeros((3, len(subjects)))
    lesion_idx = np.random.permutation(np.arange(0,2473))[0:lesion_size]             #np.arange(idx, idx+1000)
    for ic, cond in enumerate(conditions):
        for ix , res in enumerate(rlpfc_results):
            fc = res[1]
            tha_b = np.mean(rlpfc_results[ix][2][cond][:,:],axis=0)
            tha_b[lesion_idx.astype('int')]=0
            ctx_b = np.mean(rlpfc_results[ix][3][cond][:,:],axis=0)
            # correlate predicted cortical evoke vs observed
            corr[ic, ix] = np.corrcoef(zscore(np.dot(tha_b, fc)), zscore(ctx_b))[0,1]
    return corr, lesion_idx

# run simulation
pool = multiprocessing.Pool(24)
idxs = np.arange(2473*50)
lesion_results = pool.map(sim_lesion, idxs)
pool.close()
pool.join()


## unpack simulation
corr = np.zeros((3, len(subjects), len(idxs)))
lesion_idx = np.zeros((len(idxs),lesion_size))
for ix , res in enumerate(lesion_results):
    corr[:,:,ix] = res[0]
    lesion_idx[ix,:] = res[1]

for ic, cond in enumerate(conditions):
    print(cond)
    print(np.percentile(np.mean(corr[ic, :,:], axis=0), 99))
    print(np.percentile(np.mean(corr[ic, :,:], axis=0), 1))


# get pseudo lesion mask
a = masking.apply_mask(thalamus_mask, thalamus_mask)
li = lesion_idx[idx, :]
counts = np.zeros(a.shape[0])
for i in np.arange(a.shape[0]):
    counts[i] = np.sum(li==i)
n = masking.unmask(counts, thalamus_mask)





###### Test multiprocessing
# def func(a):
#     b = {'r': np.zeros(a)}
#     c = {'r': np.zeros(a)}
#     d = {'r': np.zeros(a)}
#
#     return b, c, d
# import os
# import multiprocessing
# import numpy as np
# pool = multiprocessing.Pool(os.cpu_count() // 2 )
#
# data = np.random.randint(1, 1000,100)
# results = pool.imap_unordered(func, data)
#



# end of line
