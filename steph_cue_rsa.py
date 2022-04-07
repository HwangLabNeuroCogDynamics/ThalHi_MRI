## testing Steph's RSA regression

from nilearn.image import resample_to_img
from thalpy import base
from sklearn.linear_model import RidgeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.svm import SVC
import sys
from thalhi.decoding import SubjectLssTentData
from thalhi.rsa.repsimanalyze import RSA_cues,apply_mask,get_voxels_to_exclude_SL_version
from sklearn.multiclass import OneVsRestClassifier
import os
import argparse
from thalpy.decoding import searchlight
from thalpy import masks
import nibabel as nib
import pickle
from nilearn.datasets import load_mni152_template
from nilearn.input_data import NiftiMasker
import nilearn.masking
import rsatoolbox
import numpy as np
#import pandas as pd
#import re
#import glob2 as glob

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Decode subject",
        usage="[Subject] [Classifier] [OPTIONS] ... ",
    )
    #parser.add_argument("subject", help="id of subject ie 10001")
    #parser.add_argument(
    #    "--statmethod", help="stats method to use: Options Pearson, Spearman, Mahalanobis, Crossnobis")
    parser.add_argument(
        "--njobs", help="number of cores to use: int between 1 and 80")
    parser.add_argument("--searchlight",
                        help="run searchlight, default is false",
                        default=False, action="store_true")
    return parser


def main():
    CUES = ["dcb", "fcb", "dpb", "fpb", "dcr", "fcr", "dpr", "fpr"]

    if os.path.exists("/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/"):
        THAL_HI_DIR = "/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/"
    elif os.path.exists("/Shared/lss_kahwang_hpc/data/ThalHi/"):
        THAL_HI_DIR = "/Shared/lss_kahwang_hpc/data/ThalHi/"

    # Parse command line arguments and separate subcommands
    parser = init_argparse()
    args = parser.parse_args(sys.argv[1:])
    #stats_method = args.statmethod
    num_cores = int(args.njobs)

    dir_tree = base.DirectoryTree(THAL_HI_DIR)
    unusable_subs = ['10006','10009',"10011","10015","10016","10029","10030","10034","10038","10039","10055","10061","10062","10065", "JH"]
    subjects = base.get_subjects(dir_tree.deconvolve_dir, dir_tree, excluded=unusable_subs) # name is weird, but "completed_subs" will exclude list passed to it

    #subject = next(sub for sub in subjects if sub.name == args.subject)
    subject_list = [sub.name for sub in subjects]
    sub_list = [] #[sub.name for sub in subjects]
    #sub_list=sorted(set(sub_list).difference(unusable_subs))

    if args.searchlight:
        #### if searchlight argument entered, run code below

        #### this will get two masks... one is just the brain vs not the brain and the other is the roi regions (I think)
        coritcal_mask = nib.load(masks.CORITCAL_BINARY_PATH) # get mask that will be used as the current searchlight mask

        # # # Load template for getting size and converting to niftii format
        img = nib.load(
            THAL_HI_DIR + "fmriprep/sub-10001/func/sub-10001_task-ThalHi_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
        dims=img.get_fdata().shape
        # # # Load 8 cue files and create voxel by cue matrix (stacked by subject)
        # # # Load errors from regressor for each subject too while we are at it
        os.chdir(os.path.join("/Shared","lss_kahwang_hpc","ThalHi_MRI_2020","RSA","searchlight")) # change directory to the deconvolve directory
        if os.path.exists("voxels3d_by_cuesXsubjs.npy"):
            #### if the LSS file has already been created, just load that
            voxels3d_by_cuesXsubjs = np.load("voxels3d_by_cuesXsubjs.npy")
            resids_data_list = pickle.load(open("searchlight_resids_list.p", "rb"))
            sub_list = subject_list
            print("loaded voxels3d_by_cuesXsubjs (numpy) and resids_data_list files")
            print("shape of voxels3d_by_cuesXsubjs", voxels3d_by_cuesXsubjs.shape)
            print("size of resids_data_list", len(resids_data_list))
        else:
            voxels3d_by_cuesXsubjs=np.zeros( (dims[0],dims[1],dims[2],(8*len(subject_list))) )
            resids_data_list = [] # make empty list... will be size==number of participants
            for ind_s, subject in enumerate(subjects):
                # change directory to the deconvolve directory
                subname = subject.name
                sub_list.append(subname)
                os.chdir(subject.deconvolve_dir) 
                print("\nLoading cue files for subject ",subname)
                # loop through cues and load and save them in larger 4D matrix 
                for ind_c, cue in enumerate(CUES):
                    print("Loading cue file ", cue)
                    iresp_file = os.path.join(subject.deconvolve_dir,("sub-"+subname+"_"+cue+"_FIR_MNI.nii.gz"))
                    data_file = nib.load(iresp_file) # load data file
                    fdata=data_file.get_fdata()
                    print("cue file is shape: ",fdata.shape)
                    fdata_avg = np.mean(fdata, axis=3) # take average over masked data (assumes 4th dimension is tents)
                    print("avg cue file is shape: ",fdata_avg.shape)
                    #data=nilearn.image.new_img_like(data_file, fdata_avg)
                    index = ((ind_s*len(CUES))+ind_c)
                    print(index)
                    voxels3d_by_cuesXsubjs[:,:,:,index] = fdata_avg
                # also load residuals file for this subject and mask so it's just the voxels
                print("\nloading residuals file for subject ",subname)
                resids_file = os.path.join(subject.deconvolve_dir, ("sub-"+subname+"_FIRmodel_errts_cues.nii.gz"))
                resids_data_nii = nib.load(resids_file)
                r_data = resids_data_nii.get_fdata()
                cortical_masker = NiftiMasker(coritcal_mask)
                resids_data = cortical_masker.fit_transform(resids_data_nii)
                #resids_data = nilearn.masking.apply_mask(resids_data_nii,coritcal_masker) 
                print("resids size: ",resids_data.shape,"\n")
                resids_data_list.append(resids_data)
                voxels_to_exclude = get_voxels_to_exclude_SL_version(r_data)
                # if ind_s==1:
                #     break # for testing purposes break early
            # save voxels3d_by_cuesXsubjs (numpy array) and resids_data_list (list type)
            np.save(os.path.join("/Shared","lss_kahwang_hpc","ThalHi_MRI_2020","RSA","searchlight", "voxels3d_by_cuesXsubjs.npy"), voxels3d_by_cuesXsubjs)
            print(voxels3d_by_cuesXsubjs.shape)
            pickle.dump(resids_data_list, open(os.path.join("/Shared","lss_kahwang_hpc","ThalHi_MRI_2020","RSA","searchlight", "searchlight_resids_list.p"), "wb"), protocol=4)
        
        #### this chunk just downsamples to 3m x 3m x 3m 
        #img = nib.load(
        #    THAL_HI_DIR + "fmriprep/sub-10001/func/sub-10001_task-ThalHi_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
        # # # Convert above 4D array to a niftii image
        imgs = nib.Nifti1Image(voxels3d_by_cuesXsubjs,
                               affine=img.affine, header=img.header)
        # resample imgs to 3x3x3
        template = load_mni152_template(resolution=3)

        resampled_imgs = resample_to_img(imgs, template)

        print(resampled_imgs.get_fdata().shape)

            X, A = _apply_mask_and_get_affinity(
        process_mask_coords, imgs, self.radius, True,
        mask_img=self.mask_img)
        
        #### Searchlight notes from meeting with Evan
        #   1st arg == sphere == brain_mask
        #   2nd arg == function that gets applied to each sphere == decode_cues
        # I am modifying the 2nd arg to calculate RSA instead
        # The 2nd argument will be run on the searchlight data passed to the SearchLight class (aka the data the func will have is a single search light)
        # my function: RSA_cues(searchlight_data, central_voxel_ind, sphere_voxel_inds, stats_method, cue_list, subject, noise)
        #sl_obj = searchlight.SearchLight(
        #    coritcal_mask, RSA_cues,
        #    [CUES, subject_list, resids_data_list], verbose=1, radius=10., n_jobs=num_cores)
        #    #subject_lss_data.trial_df.Cue.to_numpy(), CUES, subject_lss_data.trial_df.Run.to_numpy()], verbose=1, radius=10., n_jobs=args.cores, process_mask_img=coritcal_mask)
        #sl_obj.run(resampled_imgs)
        # Evan's output was a list
        # output is a list that stores the output from my RSA_cues function 
        print(sl_obj.output)
        #### this uses pickle to save my output
        pickle.dump(sl_obj, open(os.path.join("/Shared","lss_kahwang_hpc","ThalHi_MRI_2020","RSA","searchlight" "searchlight_rsa.p"), "wb"), protocol=4)


### here are the searchlight functions from Evan
"""
The searchlight is a widely used approach for the study of the
fine-grained patterns of information in fMRI analysis, in which
multivariate statistical relationships are iteratively tested in the
neighborhood of each location of a domain.
"""
# Authors : Vincent Michel (vm.michel@gmail.com)
#           Alexandre Gramfort (alexandre.gramfort@inria.fr)
#           Philippe Gervais (philippe.gervais@inria.fr)
#
# License: simplified BSD

import time
import sys
import warnings

import numpy as np

from joblib import Parallel, delayed, cpu_count
from sklearn.exceptions import ConvergenceWarning

from nilearn import masking
from nilearn.image.resampling import coord_transform
from nilearn.input_data.nifti_spheres_masker import _apply_mask_and_get_affinity
from nilearn._utils import check_niimg_4d
import sys


def search_light(X, func, args, A, n_jobs=-1, verbose=0):
    """Function for computing a search_light
    Parameters
    ----------
    X : array-like of shape at least 2D
        data to fit.
    func : function
        to apply to X and y
    args : list
        additional arguments for func
    A : scipy sparse matrix.
        adjacency matrix. Defines for each feature the neigbhoring features
        following a given structure of the data.
    %(n_jobs_all)s
    %(verbose0)s
    Returns
    -------
    output : list of length (number of rows in A)
    """

    group_iter = GroupIterator(A.shape[0], n_jobs)
    with warnings.catch_warnings():  # might not converge
        warnings.simplefilter('ignore', ConvergenceWarning)
        output = Parallel(n_jobs=n_jobs, verbose=verbose, backend="multiprocessing")(
            delayed(_group_iter_search_light)(
                A.rows[list_i],
                X, func, args, thread_id + 1, A.shape[0], verbose)
            for thread_id, list_i in enumerate(group_iter))
    return output


class GroupIterator(object):
    """Group iterator
    Provides group of features for search_light loop
    that may be used with Parallel.
    Parameters
    ----------
    n_features : int
        Total number of features
    %(n_jobs)s
    """

    def __init__(self, n_features, n_jobs=1):
        self.n_features = n_features
        if n_jobs == -1:
            n_jobs = cpu_count()
        self.n_jobs = n_jobs

    def __iter__(self):
        split = np.array_split(np.arange(self.n_features), self.n_jobs)
        for list_i in split:
            yield list_i


def _group_iter_search_light(list_rows, X, func, args, thread_id, total, verbose=0):
    """Function for grouped iterations of search_light
    Parameters
    -----------
    list_rows : array of arrays of int
        adjacency rows. For a voxel with index i in X, list_rows[i] is the list
        of neighboring voxels indices (in X).
    estimator : estimator object implementing 'fit'
        object to use to fit the data
    X : array-like of shape at least 2D
        data to fit.
    func : function
        to apply to X
    args : list
        additional arguments for func
    thread_id : int
        process id, used for display.
    total : int
        Total number of voxels, used for display
    verbose : int, optional
        The verbosity level. Default is 0
    Returns
    -------
    par_scores : numpy.ndarray
        score for each voxel. dtype: float64.
    """
    output = []
    t0 = time.time()
    for i, row in enumerate(list_rows):
        output.append(func(X[:, row], row, *args))

        if verbose > 0:
            # One can't print less than each 10 iterations
            step = 11 - min(verbose, 10)
            if (i % step == 0):
                # If there is only one job, progress information is fixed
                if total == len(list_rows):
                    crlf = "\r"
                else:
                    crlf = "\n"
                percent = float(i) / len(list_rows)
                percent = round(percent * 100, 2)
                dt = time.time() - t0
                # We use a max to avoid a division by zero
                remaining = (100. - percent) / max(0.01, percent) * dt
                sys.stderr.write(
                    "Job #%d, processed %d/%d voxels "
                    "(%0.2f%%, %i seconds remaining)%s"
                    % (thread_id, i, len(list_rows), percent, remaining, crlf))
                sys.stdout.flush()
                sys.stderr.flush()
    return output

##############################################################################
# Class for search_light #####################################################
##############################################################################


class SearchLight():
    """Implement search_light analysis using an arbitrary type of classifier.
    Parameters
    -----------
    mask_img : Niimg-like object
        See http://nilearn.github.io/manipulating_images/input_output.html
        boolean image giving location of voxels containing usable signals.
    func : function 
        to apply to X and y
    args : list
        additional arguments for func
    process_mask_img : Niimg-like object, optional
        See http://nilearn.github.io/manipulating_images/input_output.html
        boolean image giving voxels on which searchlight should be
        computed.
    radius : float, optional
        radius of the searchlight ball, in millimeters. Defaults to 2.
    estimator : 'svr', 'svc', or an estimator object implementing 'fit'
        The object to use to fit the data
    %(n_jobs)s
    scoring : string or callable, optional
        The scoring strategy to use. See the scikit-learn documentation
        If callable, takes as arguments the fitted estimator, the
        test data (X_test) and the test target (y_test) if y is
        not None.
    cv : cross-validation generator, optional
        A cross-validation generator. If None, a 3-fold cross
        validation is used or 3-fold stratified cross-validation
        when y is supplied.
    %(verbose0)s
    Notes
    ------
    The searchlight [Kriegeskorte 06] is a widely used approach for the
    study of the fine-grained patterns of information in fMRI analysis.
    Its principle is relatively simple: a small group of neighboring
    features is extracted from the data, and the prediction function is
    instantiated on these features only. The resulting prediction
    accuracy is thus associated with all the features within the group,
    or only with the feature on the center. This yields a map of local
    fine-grained information, that can be used for assessing hypothesis
    on the local spatial layout of the neural code under investigation.
    Nikolaus Kriegeskorte, Rainer Goebel & Peter Bandettini.
    Information-based functional brain mapping.
    Proceedings of the National Academy of Sciences
    of the United States of America,
    vol. 103, no. 10, pages 3863-3868, March 2006
    """

    def __init__(self, mask_img, func, args, process_mask_img=None, radius=2.,
                 n_jobs=1,
                 verbose=0):
        self.mask_img = mask_img
        self.func = func
        self.args = args
        self.process_mask_img = process_mask_img
        self.radius = radius
        self.n_jobs = n_jobs
        self.verbose = verbose

    def run(self, imgs):
        """Fit the searchlight
        Parameters
        ----------
        imgs : Niimg-like object
            See http://nilearn.github.io/manipulating_images/input_output.html
            4D image.

        _apply_mask_and_get_affinity(seeds, niimg, radius, allow_overlap,
                                 mask_img=None):
        example code for function at following link
        https://github.com/nilearn/nilearn/blob/757257e46d381b0933209fccf59f336ff306ee8f/nilearn/maskers/nifti_spheres_masker.py#L24
        lines 24 - 136 cover this function
        """

        # check if image is 4D
        imgs = check_niimg_4d(imgs)

        # Get the seeds
        process_mask_img = self.process_mask_img
        if self.process_mask_img is None:
            process_mask_img = self.mask_img

        # Compute world coordinates of the seeds
        process_mask, process_mask_affine = masking._load_mask_img(
            process_mask_img)
        process_mask_coords = np.where(process_mask != 0)
        process_mask_coords = coord_transform(
            process_mask_coords[0], process_mask_coords[1],
            process_mask_coords[2], process_mask_affine)
        process_mask_coords = np.asarray(process_mask_coords).T

        
        X, A = _apply_mask_and_get_affinity(
            process_mask_coords, imgs, self.radius, True,
            mask_img=self.mask_img)
        
        self.output = search_light(X, self.func, self.args, A,
                                   self.n_jobs, self.verbose)
        return self.output

#### more regression functions from Steph
from thalpy import masks
from thalpy.constants.paths import SCRIPTS_DIR
import nibabel as nib
import nilearn
import nilearn.masking
import numpy as np
import os
import sys
import pandas as pd
import pickle
import rsatoolbox
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import statistics
import statsmodels.api as sm
import time

"""
This script creates an RSA object that contains the following information
   * self - the object itself
   * sub_deconvolve_dir - path to where the 3dDeconvolve data is stored
   * sub_rsa_dir - path to where the rsa results are stored
   * cues - a list of the cues
   * rois - a list of the rois used
   * noise_list - a list of the noise objects (rsatoolbox format) used to compute RSA
   * path - the path to where this object is stored and the object file name 
This script contains necessary functions to calculate RSA (on fMRI data at least)
"""

# # # # CREATE RSA CLASS GENERATED FROM PREDEFINED ROIS
class SubjectRsaData:
    def __init__(self, sub_deconvolve_dir, sub_rsa_dir, cues, rois, path="RSA_obj.p"):
        # directories
        self.sub_deconvolve_dir = sub_deconvolve_dir
        self.sub_rsa_dir = sub_rsa_dir
        # lists of info I need for analysis
        self.cues = cues # list of cues
        self.rois = rois # list of rois used
        # noise I need for RSA toolbox
        self.noise_list_rois = None # list of noise objects for RSA
        # RSA outputs
        self.RSA_array = None # numpy array of RSA matrices for each roi

        self.__Noise_cues() # create noise w/ searchlight
        self.__RSA_cues() # create RSA matrices w/ searchlight
        self.save()

    @staticmethod
    def load(filepath):
        sys.path.append(os.path.join(SCRIPTS_DIR, "thalhi/"))
        return pickle.load(open(filepath, "rb"))

    def save(self, path=None):
        if path:
            self.path = path
        pickle.dump(self, open(self.path, "wb"), protocol=4)

# # # # CREATE RSA CLASS GENERATED FROM SEARCHLIGHT SPHERES
class SubjectRsaData_SL:
    def __init__(self, sub_deconvolve_dir, sub_rsa_dir, cues, path="RSA_obj.p"):
        # directories
        self.sub_deconvolve_dir = sub_deconvolve_dir
        self.sub_rsa_dir = sub_rsa_dir
        # lists of info I need for analysis
        self.cues = cues # list of cues
        # noise I need for RSA toolbox
        self.noise_list_sl = None # list of noise objects for earch searchlight
        # RSA outputs
        self.RSA_list = None # list of RSA matrices for each searchlight

        self.__Noise_cues() # create noise w/ searchlight
        self.__RSA_cues() # create RSA matrices w/ searchlight
        self.save()

    @staticmethod
    def load(filepath):
        sys.path.append(os.path.join(SCRIPTS_DIR, "thalhi/"))
        return pickle.load(open(filepath, "rb"))

    def save(self, path=None):
        if path:
            self.path = path
        pickle.dump(self, open(self.path, "wb"), protocol=4)

# ----- LOADING FILES/DATA ------ #
def load_3dDeconvolve(THAL_HI_DIR,subject,cue):
    iresp_file = os.path.join(THAL_HI_DIR,"3dDeconvolve",("sub-"+subject),("sub-"+subject+"_"+cue+"_FIR_MNI.nii.gz"))
    print("loading ... ", iresp_file)
    data_file = nib.load(iresp_file)
    return data_file

def load_mask(mask_file):
    print("\n\nPulling mask file from ... ", mask_file)
    mask = nib.load(mask_file)
    mask_data = mask.get_fdata()
    print(mask_data.shape)
    return mask, mask_data

def load_trial_level_data(subject,cue_list):
    os.chdir(subject.deconvolve_dir) # change directory to the deconvolve directory
    print("changed directory to ", subject.deconvolve_dir)
    if os.path.exists("LSS_TENT.p"):
        #### if the LSS file has already been created, just load that
        subject_lss_data = SubjectLssTentData.load("LSS_TENT.p")
    else:
        #### if the LSS file has not been created, grab lss files, create it, and save it
        print("Converting LSS files")
        subject_lss_data = SubjectLssTentData(subject.deconvolve_dir, cue_list)
        subject_lss_data.save()
    return subject_lss_data


# ----- MASKING FUNCTIONS ------ #
def get_binary_mask(mask_data,voxels_to_exclude,roi):
    print("number of voxels in current ROI mask:",np.where(mask_data==roi,1,0).sum())
    voxels_orginal = np.where(mask_data==roi,1,0).sum()
    mask_binary=np.where(np.logical_and(mask_data==roi,voxels_to_exclude==0),1,0) # make sure to also exclude voxels with all zeros
    print("number of usable voxels in current ROI mask:",mask_binary.sum(),"\n")
    voxels_usable = mask_binary.sum()
    return mask_binary, voxels_orginal, voxels_usable

def apply_mask(data_file,mask_binary):
    print("data size: ", data_file.get_fdata().shape, "\napplying mask ...")
    mask_binary_nif=nilearn.image.new_img_like(data_file, mask_binary)
    masked_data=nilearn.masking.apply_mask(data_file,mask_binary_nif) 
    #masked_data=input_data.NiftiMasker(mask_binary_nif)
    print("masked data is a numpy array of size: ", masked_data.shape)
    return masked_data


# ----- USABLE/EXLUDE FUNCTIONS ------ #
def get_voxels_to_exclude(resids_file):
    print("loading ... ", resids_file)
    r=nib.load(resids_file)
    print(r.shape)
    # check for zeros
    r_data = r.get_fdata()
    voxels_to_exclude = np.zeros((r_data.shape[0], r_data.shape[1], r_data.shape[2])) # initialize 3D matrix of what voxels to exclude
    for x in range(r_data.shape[0]):
        for y in range(r_data.shape[1]):
            for z in range(r_data.shape[2]):
                if r_data[x,y,z,:].sum()==0:
                    # voxel had 0 for all time points ... exclude from further analysis
                    voxels_to_exclude[x,y,z]=1
    print("A total of",voxels_to_exclude.sum(),"voxels will be EXCLUDED due to 0s for all time points")
    print(voxels_to_exclude.sum(),"voxels exluded out of",(r_data.shape[0]*r_data.shape[1]*r_data.shape[2]),"total voxels")
    #proportion_excluded=voxels_to_exclude.sum()/(r_data.shape[0]*r_data.shape[1]*r_data.shape[2])
    print((voxels_to_exclude.sum()/(r_data.shape[0]*r_data.shape[1]*r_data.shape[2])),"proportion of voxels excluded\n")
    return voxels_to_exclude

def remove_censored_data(resids_file,mask_binary):
    #    CHECK FOR CENSORED PERIODS IN VOXELS and remove censored portions
    #    note, VOXELSxCUES will be 2D [voxels in current ROI x 8 cues] AND this will identify periods of time where voxels are censored (should have some usable times)
    print("loading ... ", resids_file)
    r=nib.load(resids_file)
    print(r.shape)
    #    FIRST - apply mask to resids file
    mask_binary_nif=nilearn.image.new_img_like(r,mask_binary)
    resids_data=nilearn.masking.apply_mask(r,mask_binary_nif) 
    print("masked residual data is a numpy array of size: ", resids_data.shape) # gives you [time by voxels]
    #    SECOND - look through good voxels and find censored times and remove them
    reduced_resids_data=[]
    for tpt in range(resids_data.shape[0]):
        if resids_data[tpt,1]!=0:
            reduced_resids_data.append(resids_data[tpt,:])
    reduced_resids_data=np.array(reduced_resids_data)
    print("masked residual data is NOW a numpy array of size: ", reduced_resids_data.shape)
    return reduced_resids_data
    
def get_voxels_to_exclude_SL_version(r_data):
    voxels_to_exclude = np.zeros((r_data.shape[0], r_data.shape[1], r_data.shape[2])) # initialize 3D matrix of what voxels to exclude
    for x in range(r_data.shape[0]):
        for y in range(r_data.shape[1]):
            for z in range(r_data.shape[2]):
                if r_data[x,y,z,:].sum()==0:
                    # voxel had 0 for all time points ... exclude from further analysis
                    voxels_to_exclude[x,y,z]=1
    print("A total of",voxels_to_exclude.sum(),"voxels will be EXCLUDED due to 0s for all time points")
    print(voxels_to_exclude.sum(),"voxels exluded out of",(r_data.shape[0]*r_data.shape[1]*r_data.shape[2]),"total voxels")
    #proportion_excluded=voxels_to_exclude.sum()/(r_data.shape[0]*r_data.shape[1]*r_data.shape[2])
    print((voxels_to_exclude.sum()/(r_data.shape[0]*r_data.shape[1]*r_data.shape[2])),"proportion of voxels excluded\n")
    return voxels_to_exclude

def remove_censored_data_noise_version(resids_data):
    reduced_resids_data=[]
    for tpt in range(resids_data.shape[0]):
        if resids_data[tpt,1]!=0:
            reduced_resids_data.append(resids_data[tpt,:])
    reduced_resids_data=np.array(reduced_resids_data)
    #print("masked residual data is NOW a numpy array of size: ", reduced_resids_data.shape)
    return reduced_resids_data

# ----- FORMATTING/CONVERSION FUNCTIONS ------ #
def make_rsatoolbox_format(CONDxVOXELS,sub,roi):
    #    ---------------------------------------------
    #    descriptors: subject = #
    #                 ROI = #
    #    channel_descriptors: voxels = ['voxel_0' 'voxel_1' 'voxel_2' ... 'voxel_N']
    #    obs_descriptors: measure = ['measure_0' 'measure_1' ... 'measure_N']
    #    noise: 
    #       * will need for mahalanobis or crossnobis options
    #       * can be the residual file from the regression
    #    cv_descriptor: trials = [0 1 2 3 4 0 1 2 3 4 ... 0 1 2 3 4]
    #       * will need for crossnobis option (tells us what trial a value came from)
    #       * will be same length as obs_descriptors and channel_descriptors
    #    ---------------------------------------------
    des={'subject': sub, 'ROI': roi}
    chn_des={'voxels': np.array(['voxel_' + str(x) for x in np.arange(CONDxVOXELS.shape[1])])}
    obs_des={'measure': np.array(['cue_' + str(x) for x in np.arange(CONDxVOXELS.shape[0])])} # ['trial_activity','cue_template']} #
    #print("input CONDxVOXELS shape: ",CONDxVOXELS.shape)
    rsatoolbox_formatted_data=rsatoolbox.data.Dataset(CONDxVOXELS, descriptors=des, obs_descriptors=obs_des, channel_descriptors=chn_des)
    return rsatoolbox_formatted_data

def convert_numpy_to_niftii(numpy_formatted_data, niftii_to_use_as_template):
    nif_format = nilearn.image.new_img_like(niftii_to_use_as_template,numpy_formatted_data)
    return nif_format

def get_noise_metric(resids_file,mask_binary):
    # # # # GET NOISE
    reduced_resids_data = remove_censored_data(resids_file,mask_binary)
    #print("reduced_resids_data\n",reduced_resids_data)
    print("Reduced residual data is size: ", reduced_resids_data.shape)
    try:
        #noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, dof=72, method='full')
        noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, method='diag')
        print("noise_pres_res shape: ",noise_pres_res.shape)
    except:
        noise_pres_res = np.nan
    return noise_pres_res

# def get_noise_metric_noise_version(central_voxel_ind, sphere_voxel_inds, resids_file):
#     # # # # GET NOISE
#     reduced_resids_data = remove_censored_data_noise_version(resids_file)
#     #print("reduced_resids_data\n",reduced_resids_data)
#     print("Reduced residual data is size: ", reduced_resids_data.shape)
#     try:
#         #noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, dof=72, method='full')
#         noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, method='diag')
#         print("noise_pres_res shape: ",noise_pres_res.shape)
#     except:
#         noise_pres_res = np.nan
#     return noise_pres_res


# ----- ANALYSIS SETUP FUNCTIONS ------ #
def get_cond_inds(cue_list):
    filled_inds = [] # 1=hollow, 2=filled
    hollow_inds = []
    for ind, cue in enumerate(cue_list):
        if cue[0] == "d":
            hollow_inds.append(ind)
        elif cue[0] == "f":
            filled_inds.append(ind)
    print("Context condition order in cue list: ",filled_inds,hollow_inds)
    return filled_inds, hollow_inds

def calculate_RSA(data_rsatool_format,noise_pres_res,roi,trl):
    #print(len(noise_pres_res.shape))
    if len(noise_pres_res.shape) > 1:
        try: 
            # # # # Calculate RSA
            # run rsatoolbox to get coeff matrix
            rdm = rsatoolbox.rdm.calc_rdm(data_rsatool_format, method='mahalanobis', descriptor='measure', noise=noise_pres_res)
            #print(rdm)
            Coeff_mat = rdm.get_matrices()[0] # the distance matrix
            print("\nCoeff_mat as 1-rdm\n\n",(Coeff_mat))
            print("\nCoefficient matrix calculated for ROI", roi, "and TRIAL", trl)
        except:
            CUESxVOXELS = VOXELSxCUES.T
            Coeff_mat = np.empty((CUESxVOXELS.shape[0],CUESxVOXELS.shape[0]))
            Coeff_mat[:] = np.nan
            #ROIs_with_inversion_error.append(roi)
            print("\n Coefficient matrix NOT calculated for ROI", roi, "and TRIAL", trl)
            print("UNKNOWN ERROR ... coefficient matrix saved as NaN matrix and ROI recorded") 
    else:
        CUESxVOXELS=VOXELSxCUES.T
        Coeff_mat = np.empty((CUESxVOXELS.shape[0],CUESxVOXELS.shape[0]))
        Coeff_mat[:] = np.nan
        #ROIs_with_inversion_error.append(roi)
        print("\n Coefficient matrix NOT calculated for ROI", roi, "and TRIAL", trl)
        print("INVERSION ERROR ... coefficient matrix saved as NaN matrix and ROI recorded") 
    return Coeff_mat
    
def create_matching_and_mismatching(masked_trl_lvl_data_avg,filled_template,hollow_template):
    # set up an array of zeros to fill in below
    matching_mat = np.zeros((filled_template.shape[0],2))
    mismatching_mat = np.zeros((filled_template.shape[0],2))
    # fill in matching and mismatching arrays
    matching_mat[:,0] = masked_trl_lvl_data_avg
    mismatching_mat[:,0] = masked_trl_lvl_data_avg
    if trl_cue[0] == "f":
        print("current trial is a filled cue type") # if filled trl
        matching_mat[:,1] = filled_template
        mismatching_mat[:,1] = hollow_template
    elif trl_cue[0] == "d":
        print("current trial is a hollow cue type") # if hollow trl
        matching_mat[:,1] = hollow_template
        mismatching_mat[:,1] = filled_template
    #print("Matching matrix shape: ",matching_mat.shape)
    #print("Mismatching matrix shape: ",mismatching_mat.shape)
    matching_mat=matching_mat.T
    mismatching_mat=mismatching_mat.T
    return matching_mat, mismatching_mat




# --------------------------------------------------------------------- #
# ------------------------- Set up Models ----------------------------- #
# --------------------------------------------------------------------- #
# Note, RSA models have 1=more similar and 0=less similar
def extract_char(list_of_strings,index_to_pull_from):
    element_list = [e[index_to_pull_from] for e in list_of_strings]
    return element_list

def gen_version_models(version_completed,element_order):
    # -- Get relevant feature model
    RelevantFeature_model = np.zeros((len(element_order),len(element_order)))
    rel_feat_vec=[]
    for indx, cur_cue in enumerate(element_order):
        # check if current cue was d | f ...
        if cur_cue[0]==version_completed[0]:
            # check what the relevant feature was
            if version_completed[1]=="c":
                rel_feat_vec.append("C") # this cue focused on color
            else:
                rel_feat_vec.append("S") # this cue focused on shape
        elif cur_cue[0]==version_completed[2]:
            # check what the relevant feature was
            if version_completed[3]=="c":
                rel_feat_vec.append("C") # this cue focused on color
            else:
                rel_feat_vec.append("S") # this cue focused on shape
    for ccol, col_cue in enumerate(element_order):
        for crow, row_cue in enumerate(element_order):
            # get the currently relevant feature
            rel_feat = rel_feat_vec[ccol]
            if rel_feat == "C":
                if col_cue[2]==row_cue[2]:
                    RelevantFeature_model[crow][ccol] = 1          
            elif rel_feat == "S":
                if col_cue[1]==row_cue[1]:
                    RelevantFeature_model[crow][ccol] = 1
    #print("\n",version_completed,"Relevant Feature Model\n",RelevantFeature_model)

    # -- Get task performed model
    TaskPerformed_model = np.zeros((len(element_order),len(element_order)))
    if version_completed=="dcfs":
        # hard coding order because I didn't save that info anywhere
        task_performed_order = ["F","S","F","S","S","F","S","F"]
    elif version_completed=="dsfc":
        task_performed_order = ["S","S","F","F","S","S","F","F"]
    for ccol, col_task in enumerate(task_performed_order):
        for crow, row_task in enumerate(task_performed_order):
            if col_task == row_task:
                TaskPerformed_model[crow][ccol] = 1          
            elif col_task == row_task:
                TaskPerformed_model[crow][ccol] = 1
    #print("\n",version_completed,"Task Performed Model\n",TaskPerformed_model)

    return RelevantFeature_model, TaskPerformed_model

def gen_RSA_models(element_order):
    # This function takes a list that tells us the order of the elements in the RSA coeff matrix
    # Note, this is pretty specific to ThalHi and it's cue list...
    # For ThalHi, the cues have 3 pieces of info
    #   1. d | f ... donut  | filled
    #   2. c | p ... circle | polygon
    #   3. r | b ...  red   |  blue
    #print(extract_char(element_order,0),"\n",extract_char(element_order,1),"\n",extract_char(element_order,2)) 

    # Set up basic models that don't depend on the version (Context & Identity)
    data_Ctxt = np.zeros((len(element_order),len(element_order)))
    for crow, celmnt1 in enumerate(extract_char(element_order,0)):
        for ccol, celmnt2 in enumerate(extract_char(element_order,0)):
            if extract_char(element_order,0)[crow]==extract_char(element_order,0)[ccol]:
                data_Ctxt[crow][ccol] = 1
    #print("\nContext Model \n",data_Ctxt)

    data_Iden = np.identity(len(element_order))
    #print("\nIdentity Model \n",data_Iden)

    # Set up models that depend on version (Relevant_Feature & Task_Performed)
    dcfs_relF, dcfs_tskR = gen_version_models("dcfs",element_order)
    dsfc_relF, dsfc_tskR = gen_version_models("dsfc",element_order)

    # flip 0s and 1s so that it can work with the distance measure RSA
    data_Ctxt = (data_Ctxt*-1)+1
    data_Iden = (data_Iden*-1)+1
    dcfs_relF = (dcfs_relF*-1)+1
    dcfs_tskR = (dcfs_tskR*-1)+1
    dsfc_relF = (dsfc_relF*-1)+1
    dsfc_tskR = (dsfc_tskR*-1)+1

    return data_Ctxt, data_Iden, dcfs_relF, dcfs_tskR, dsfc_relF, dsfc_tskR

def set_up_model_vars(mod2skip):
    if mod2skip == "TP":
        regressor_list=["Intercept","Context","Relevant_Feature"]
    elif mod2skip == "RF":
        regressor_list=["Intercept","Context","Task_Performed"]
    else: 
        regressor_list=["Intercept","Context","Relevant_Feature","Task_Performed"]
    stat_list=["_Beta","_T-stat","_P-value"]
    #print("\nGenerating column headers based on given regressors and stats to save out ... \n")
    column_headers=[]
    for rr in regressor_list:
        for ss in stat_list:
            column_headers.append((rr+ss))
    #print(column_headers,"\n")

    temp_vec=np.tril(np.random.rand(8,8), k=-1).flatten() # k=-1 should reduce to only entries below the diagonal
    lower_triangle_inds=np.where(temp_vec!=0)[0] # will use this to pull out the lower triangle from the coeff mats generated below
    #print("length of lower triangle inds vector",len(lower_triangle_inds))

    return regressor_list, stat_list, column_headers, lower_triangle_inds

def get_start_and_end(ind, len_of_data_to_add):
    start_ind = ind*len(len_of_data_to_add)
    stop_ind = start_ind+len(len_of_data_to_add)
    return start_ind, stop_ind

def create_model_vecs(sub_list,lower_triangle_inds):
    y_vec = np.zeros((len(sub_list)*len(lower_triangle_inds))) # will be size num_subs*correlation_pairs (not included correlations with itself)
    context_vec=np.zeros((len(sub_list)*len(lower_triangle_inds))) 
    relFeat_vec=np.zeros((len(sub_list)*len(lower_triangle_inds))) 
    taskPerform_vec=np.zeros((len(sub_list)*len(lower_triangle_inds)))
    return y_vec, context_vec, relFeat_vec, taskPerform_vec

def create_regressor_dataframe(mod2skip, regressor_list, y_vec, context_vec, relFeat_vec, taskPerform_vec):
    # ---- create data frame for regressors
    if mod2skip == "TP":
        data=np.array([y_vec,np.ones((len(y_vec))),context_vec,relFeat_vec]).T
        df=pd.DataFrame(data=data, columns=["Y","Intercept","Context","Relevant_Feature"])
    elif mod2skip == "RF":
        data=np.array([y_vec,np.ones((len(y_vec))),context_vec,taskPerform_vec]).T
        df=pd.DataFrame(data=data, columns=["Y","Intercept","Context","Task_Performed"])
    else: 
        data=np.array([y_vec,np.ones((len(y_vec))),context_vec,relFeat_vec,taskPerform_vec]).T
        df=pd.DataFrame(data=data, columns=["Y","Intercept","Context","Relevant_Feature","Task_Performed"])
    # ---- save out data frame for current ROI (usefull for double checking data input to model if needed)
    #ROI_dataset_filename = "ROI_" + str(roi) + "_Dataset_" + stats_method + "_argon.csv"
    #df.to_csv(os.path.join(data_dir,"RSA","ROI_Level_Datasets",ROI_dataset_filename))
    # ---- run stats for current roi
    model = sm.OLS(df["Y"], df[regressor_list])
    res = model.fit()
    # print("summary",res.summary())
    # print("params",res.params)
    # print("tvalues",res.tvalues)
    # print("pvalues",res.pvalues)
    # ---- Pull out Betas, T-stats, and P-values
    #      Note, order is: [0]=Intercept [1]=Context, [2]=Relevant_Feature, [3]=Task_Performed
    cur_row_list=[]
    for cur_regressor in range(0,len(regressor_list)):
        cur_row_list.append(res.params[cur_regressor])
        cur_row_list.append(res.tvalues[cur_regressor])
        cur_row_list.append(res.pvalues[cur_regressor])
    #print("cur_row_list",cur_row_list)
    cur_dict = {'context_beta': cur_row_list[0],'context_tval': cur_row_list[1],'context_pval': cur_row_list[2],
                'relfeat_beta': cur_row_list[3],'relfeat_tval': cur_row_list[4],'relfeat_pval': cur_row_list[5],
                'taskper_beta': cur_row_list[6],'taskper_tval': cur_row_list[7],'taskper_pval': cur_row_list[8]}
    return cur_dict




# ------------------------------------------------- #
# ----- INTEGRATE FUNCTIONS TO CALCULATE RSA ------ #
def RSA_cues(searchlight_data, sphere_voxel_inds, cue_list, sub_list, resids_data_list):
    # initialize some directories
    if os.path.exists("/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/"):
        mnt_dir = "/mnt/cifs/rdss/rdss_kahwang/ThalHi_data/MRI_data/"
    elif os.path.exists("/Shared/lss_kahwang_hpc/data/ThalHi/"):
        mnt_dir = "/Shared/lss_kahwang_hpc/ThalHi_MRI_2020/"
    # initialize lists, arrays, and variables
    mod2skip=""
    regressor_list, stat_list, column_headers, lower_triangle_inds = set_up_model_vars(mod2skip)
    y_vec, context_vec, relFeat_vec, taskPerform_vec = create_model_vecs(sub_list,lower_triangle_inds)
    # initialize model info
    data_Ctxt, data_Iden, dcfs_relF, dcfs_tskR, dsfc_relF, dsfc_tskR = gen_RSA_models(cue_list)
    version_info = pd.read_csv(os.path.join(mnt_dir, "Version_Info.csv"))
    print(sub_list)
    print(version_info)
    # loop through subjects, get rsa on sphere, and run model on all subjects
    for ind, subject in enumerate(sub_list):
        # # # # Loop where we reduce to just current subejct to get RSA for each subj
        #print("\nsearchlight data is size",searchlight_data.shape,"\nsphere voxel arg is size",len(sphere_voxel_inds)) # searchlight data = [subject_x_cues , voxels]
        start_ind, stop_ind = get_start_and_end(ind, cue_list)
        #print("start:",start_ind,"\tstop:",stop_ind)
        cur_sub_data = searchlight_data[start_ind:stop_ind,:]
        #print("searchlight data for current subject is size",cur_sub_data.shape) # cur_sub_data = [cues , voxels]
        resids_data = resids_data_list[ind]
        #                     remove_censored_data_noise_version(resids_data)
        reduced_resids_data = remove_censored_data_noise_version(resids_data[:,sphere_voxel_inds])
        noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, method='diag')
        #                     make_rsatoolbox_format(CONDxVOXELS, sub,  roi)
        data_rsatool_format = make_rsatoolbox_format(cur_sub_data,subject,1) # entering 1 as roi (b/c we don't need it)
        # # # # Calculate RSA
        # RSA will be calculated for each subject... the coeffs will be pulled out and added to a vector
        # the stats model will then be applied to the vector of all subject coeffs for this searchlight
        try: 
            rdm = rsatoolbox.rdm.calc_rdm(data_rsatool_format, method='mahalanobis', descriptor='measure', noise=noise_pres_res)
            Coeff_mat = rdm.get_matrices()[0] # the distance matrix
        except:
            Coeff_mat = np.empty((cur_sub_data.shape[0],cur_sub_data.shape[0]))
            Coeff_mat[:] = np.nan
            #ROIs_with_inversion_error.append(roi)
        #print("Coeff mat is size:",Coeff_mat.shape,"\nCoeff mat looks like:",Coeff_mat)
        # # # # Pull out lower triangle from coeff mat for this subject
        coeff_vec=np.tril(Coeff_mat, k=-1).flatten()
        start_pt, end_pt = get_start_and_end(ind, lower_triangle_inds)
        #print("start point =",start_pt,"\tend point =",end_pt)
        y_vec[start_pt:end_pt]=coeff_vec[lower_triangle_inds] # add to y vector for model
        #    set up other model vectors based on what version the current sub did
        context_vec[start_pt:end_pt]=np.tril(data_Ctxt).flatten()[lower_triangle_inds]
        if version_info["version"][ind]=="DCFS":
            #    not swapped version
            relFeat_vec[start_pt:end_pt]=np.tril(dcfs_relF).flatten()[lower_triangle_inds]
            taskPerform_vec[start_pt:end_pt]=np.tril(dcfs_tskR).flatten()[lower_triangle_inds]
        else:
            relFeat_vec[start_pt:end_pt]=np.tril(dsfc_relF).flatten()[lower_triangle_inds]
            taskPerform_vec[start_pt:end_pt]=np.tril(dsfc_tskR).flatten()[lower_triangle_inds]
    #                         create_regressor_dataframe(mod2skip, regressor_list, y_vec, context_vec, relFeat_vec, taskPerform_vec)
    searchlight_stat_output = create_regressor_dataframe(mod2skip, regressor_list, y_vec, context_vec, relFeat_vec, taskPerform_vec)
    print(searchlight_stat_output) # dictionary type

    # output variable is list format
    return searchlight_stat_output

def searchlight_tester(searchlight_data, sphere_voxel_inds, start_time):
    cur_time = time.time() - start_time
    # just return the current voxel
    return cur_time





if __name__ == "__main__":
    main()