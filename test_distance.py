# import packages we'll need
#import matplotlib as plt
import nilearn
import nibabel as nib # requires some packages (listed below import list)
import numpy as np
from scipy import stats
import pandas as pd
import re
import glob
import os
import sys
import nilearn.masking
import rsatoolbox
from thalpy import base
import argparse

"""
This script is meant to load iresp output from 3dDeconvolve and calcualte RSA

to call and run from a bash terminal, run as follows
 python create_RSA_Arrays.py [subject] --statsmethod [stats_method] --datadir [data_dir]

[subject] = the current subject to run 
            formatted as a five digit ID (i.e., 10001)
[stats_method] = the current stats method to use
                 options: Pearson  Spearman  Mahalanobis  *Crossnobis*
                 *currently in progress and not fully supported yet
[data_dir] = the directory for the data that will be used
             formatted as a typical data path (i.e., /full/path/to/data/folder/)

Last edited: Jan 7th 2022
         by: Stephanie C Leach
"""

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Decode subject",
        usage="[Subject] [--STATSMETHOD]",
    )
    parser.add_argument("subject", help="id of subject ie 10001")
    parser.add_argument("--statsmethod", help="stats_method to use: Options Pearson Spearman Mahalanobis")
    parser.add_argument("--datadir", help="data directory: e.g., /data/backed_up/shared/ThalHi_MRI_2020/")
    return parser
#     unusable_subs = ['10006','10009',"10011","10015","10016","10029","10030","10034","10038","10039","10055","10061","10062","10065"]

def load_and_apply_mask(data_filename,mask_binary):
    #    load file
    print("loading ... ", data_filename)
    data_file=nib.load(data_filename)
    print("data size: ", data_file.get_fdata().shape, "\n")
    print("applying mask ...")
    mask_binary_nif=nilearn.image.new_img_like(data_file, mask_binary)
    masked_data=nilearn.masking.apply_mask(data_file,mask_binary_nif) 
    print("masked data is a numpy array of size: ", masked_data.shape)
    return masked_data

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

def make_rsatoolbox_format(VOXELSxCUES,sub,roi,cue_list):
    #    ---------------------------------------------
    #    descriptors: subject = #
    #                 ROI = #
    #    channel_descriptors: voxels = ['voxel_0' 'voxel_1' 'voxel_2' ... 'voxel_N']
    #    obs_descriptors: conds = ['cond_0' 'cond_1' 'cond_2' ... 'cond_N']
    #    noise: 
    #       * will need for mahalanobis or crossnobis options
    #       * can be the residual file from the regression
    #    cv_descriptor: trials = [0 1 2 3 4 0 1 2 3 4 ... 0 1 2 3 4]
    #       * will need for crossnobis option (tells us what trial a value came from)
    #       * will be same length as obs_descriptors and channel_descriptors
    #    ---------------------------------------------
    des={'subject': sub, 'ROI': roi}
    chn_des={'voxels': np.array(['voxel_' + str(x) for x in np.arange(VOXELSxCUES.shape[0])])}
    obs_des={'cues': cue_list} #np.array(['cue_' + str(x) for x in np.arange(VOXELSxCUES.shape[1])])}
    CUESxVOXELS=VOXELSxCUES.T
    print("CUESxVOXELS shape: ",CUESxVOXELS.shape)
    rsatoolbox_formatted_data=rsatoolbox.data.Dataset(CUESxVOXELS,descriptors=des,obs_descriptors=obs_des,channel_descriptors=chn_des)
    return rsatoolbox_formatted_data

def apply_mask_get_coeff(data_dir,sub, cue_list,mask_data,corr_method):
    # # set up an empty array/list to keep track of ROIs with the inversion error
    ROIs_with_inversion_error=[]
    # # set up empty matrix for filling in ROI by cue by cue (400 X 8 X 8)
    ROIxCOEFF_mat=np.zeros((400,len(cue_list),len(cue_list)))
    
    # # load residuals file and run some checks (for zeros)
    #   note, this should identify voxels that were outside of brain and always 0
    resids_file = os.path.join(data_dir,"3dDeconvolve/",("sub-"+sub),("sub-"+sub+"_FIRmodel_errts_cues.nii.gz"))
    voxels_to_exclude = get_voxels_to_exclude(resids_file,)
    
    # make a list of voxels originally in roi and voxel in roi after excluding zeros
    voxels_orginal=[]
    voxels_usable=[]
    for roi in range(319,401):
        print("\ncurrently on ROI:", roi)
        
        #    pull out current roi from mask AND get usable voxels
        print("number of voxels in current ROI mask:",np.where(mask_data==roi,1,0).sum())
        voxels_orginal.append(np.where(mask_data==roi,1,0).sum())
        mask_binary=np.where(np.logical_and(mask_data==roi,voxels_to_exclude==0),1,0) # make sure to also exclude voxels with all zeros
        print("number of usable voxels in current ROI mask:",mask_binary.sum(),"\n")
        voxels_usable.append(mask_binary.sum())
        
        #    create an empty matrix of size VOXELSx8 (VOXELSxCUES)
        num_voxels=(mask_binary==1).sum()
        VOXELSxCUES=np.zeros((len(cue_list),num_voxels)).T
        print("voxels in current ROI:", num_voxels, "  VOXELSxCUES matrix dimensions:", VOXELSxCUES.shape)

        #    now loop through the 8 cues in cue_list
        for cue in cue_list:
            #    load iresp file for current cue AND run function that outputs masked data
            iresp_file = os.path.join(data_dir,"3dDeconvolve/",("sub-"+sub),("sub-"+sub+"_"+cue+"_FIR_MNI.nii.gz"))
            masked_data = load_and_apply_mask(iresp_file,mask_binary) 
            masked_data_avg=np.mean(masked_data, axis=0) # take average over masked data
            #print(masked_data_avg.shape)
            #    add 1D masked data vector (avg) to my numpy matrix of size VOXELSxCUES
            VOXELSxCUES[:,(cue_list.index(cue))]=masked_data_avg
            print(VOXELSxCUES.shape)
        
        print(VOXELSxCUES[:,[0, 1, 2, 3]])
        print(VOXELSxCUES[:,[4, 5, 6, 7]])

        #    Now that the entire Cue by Voxel matrix has been filled in, 
        #    calculate correlation coefficients
        if (corr_method=="Pearson"):
            #    use pearson correlation
            Coeff_mat=np.corrcoef(VOXELSxCUES, rowvar=False)
        elif (corr_method=="Spearman"):
            #    use spearman correlation
            Coeff_mat=stats.spearmanr(VOXELSxCUES,axis=0)[0] # get just coefficient matrix
        elif (corr_method=="Mahalanobis"):
            data_rsatool_format = make_rsatoolbox_format(VOXELSxCUES,sub,roi,cue_list)
            reduced_resids_data = remove_censored_data(resids_file,mask_binary)
            print("reduced_resids_data\n",reduced_resids_data)
            print(reduced_resids_data.shape)
            noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, dof=72, method='full')
            print(noise_pres_res)
            try: 
                #    get noise information for mahalanobis option
                #noise_pres_res = rsatoolbox.data.noise.prec_from_residuals(reduced_resids_data, dof=72, method='full')
                print("noise_pres_res shape: ",noise_pres_res.shape)
                #    run rsatoolbox to get coeff matrix
                rdm=rsatoolbox.rdm.calc_rdm(data_rsatool_format, descriptor='cues', method='mahalanobis', noise=noise_pres_res)
                print(rdm)
                #Coeff_mat=1-(rdm.get_matrices()[0]) # 1 - the dissimilarity matrix
                Coeff_mat=rdm.get_matrices()[0]
                print("\n\nCoeff_mat as 1-rdm\n\n",(Coeff_mat))
                print(Coeff_mat.shape)
                print("\n Coefficient matrix calculated for ROI", roi)
            except:
                CUESxVOXELS=VOXELSxCUES.T
                Coeff_mat = np.empty((CUESxVOXELS.shape[0],CUESxVOXELS.shape[0]))
                Coeff_mat[:] = np.nan
                ROIs_with_inversion_error.append(roi)
                print("\n Coefficient matrix NOT calculated for ROI", roi)
                print("INVERSION ERROR ... coefficient matrix saved as NaN matrix and ROI recorded") 
        elif (corr_method=="Crossnobis"):
            # figure out later... need trial level first
            #data_rsatool_format = make_rsatoolbox_format(VOXELSxCUES,sub,roi,cue_list)
            #reduced_resids_data = remove_censored_data(r,mask_binary)
            #Coeff_mat=rsatoolbox.rdm.calc.calc_rdm_crossnobis(data_rsatool_format,descriptor='cues')
            #Coeff_mat=rsatoolbox.rdm.calc_rdm(data_rsatool_format,method='crossnobis',descriptor='cues',noise=None,cv_descriptor=None)
            print("Crossnobis currently not supported...")

        print(Coeff_mat)
        "crash"+2
        #    add coefficient matrix to ROI by coefficient matrix array
        ROIxCOEFF_mat[(roi-1),:,:]=Coeff_mat
    # return calculated variables of interest
    return ROIxCOEFF_mat, ROIs_with_inversion_error, voxels_orginal, voxels_usable

def main():
    # # # # run parser to get input arguments
    parser = init_argparse()
    args = parser.parse_args(sys.argv[1:])
    sub = args.subject
    stats_method = args.statsmethod
    print("working on ...",sub,"\nusing stats method:",stats_method)

    # # # # set up data directory
    data_dir = args.datadir

    # # # # manually set cue list ... cue_list=["fpr","fpb","fcr","fcb","dpr","dpb","dcr","dcb"]
    cue_list=["dcr","dcb","dpr","dpb","fcr","fpr","fcb","fpb"]
    
    # # # # load mask file
    mask_file = os.path.join(data_dir,"ROIs","Schaefer400_2.5.nii.gz")
    print("\n\nPulling mask file from ... ", mask_file)
    mask=nib.load(mask_file)
    mask_data=mask.get_fdata()
    print(mask_data.shape)
    
    # # # # apply mask to data and get coeeficient matrix for each ROI
    ROIxCOEFF_mat, ROIs_with_inversion_error, voxels_orginal, voxels_usable = apply_mask_get_coeff(data_dir, sub, cue_list, mask_data, stats_method)
    
    # # # # SAVE OUT FINAL MATRIX
    output_filename="sub-"+sub+"_ROI_by_"+stats_method+"_Coeff_Matrix.npy"
    np.save(os.path.join(data_dir,"RSA/Subject_Level_Arrays/",output_filename),ROIxCOEFF_mat)
    
    # # # # save out list of ROIs with inversion issues as a numpy array
    ROIs_with_inversion_error = np.asarray(ROIs_with_inversion_error)
    if ROIs_with_inversion_error.shape[0]==0:
        print("cue order: ",cue_list,"\nJust finished running subject",sub)
    else:
        array_filename="sub-"+sub+"_ROIs_with_Inversion_Error.npy"
        np.save(os.path.join(data_dir,"RSA/Subject_Level_Arrays/",array_filename),ROIs_with_inversion_error)
        print("cue order: ",cue_list,"\nJust finished running subject",sub)

    # # # # save out information on roi voxels as .csv files
    voxel_info=pd.DataFrame({'ROI':range(1,401),'Voxel_Num_Originally': voxels_orginal,'Voxel_Num_Usable':voxels_usable})
    voxel_info.to_csv(os.path.join(data_dir,"RSA/Subject_Level_Arrays/",("sub-"+sub+"voxel_info.csv")))
    print(voxel_info)

main() # run the main function to execute code