# import packages we'll need
#import matplotlib as plt
import nilearn
import nibabel as nib # requires some packages (listed below import list)
import numpy as np
import pandas as pd
import statsmodels.api as sm
import re
import glob
import os
import nilearn.masking

"""
This file is meant to load the 4D numpy array (generated by agg_resp.py), 
convert it to a pandas data frame, run stats on that data frame, and then 
save out the betas, t-stats, and p-values for each ROI in a data frame and
in a nifti image (for easier visualization of sig. ROIs)

psedo code:
>> set up models for regression
>> set up directories, lists, etc.
>> load np array ...   np.load("filename")  
>> for roi in roi_list
>>    set up lists for current roi
>>    for cue in cuelist
>>       create vectors for regression
>>    run statistics for current roi
>>    save out currrent roi stats (will save in data frame)
>> save out data frame with stats for each roi
>> save out a nifti image for each stat (t-stat, p-value, and betas) for each condition
"""
# ---------------------- Set up Models ------------------------- #
# 0 = assumed to be more similar
# 1 = assumed to be less similar
data_Ctxt=np.array([[0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 1, 1, 1, 1], 
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [0, 0, 0, 0, 1, 1, 1, 1],
                    [1, 1, 1, 1, 0, 0, 0, 0],
                    [1, 1, 1, 1, 0, 0, 0, 0],
                    [1, 1, 1, 1, 0, 0, 0, 0],
                    [1, 1, 1, 1, 0, 0, 0, 0]])

data_Iden=np.array([[-1, 0, 0, 0, 0, 0, 0, 0],
                    [0, -1, 0, 0, 0, 0, 0, 0], 
                    [0, 0, -1, 0, 0, 0, 0, 0],
                    [0, 0, 0, -1, 0, 0, 0, 0],
                    [0, 0, 0, 0, -1, 0, 0, 0],
                    [0, 0, 0, 0, 0, -1, 0, 0],
                    [0, 0, 0, 0, 0, 0, -1, 0],
                    [0, 0, 0, 0, 0, 0, 0, -1]])+1


# -- # DCFS
# 0 means looking at same feature (color or shape) of cue
# 1 means looking at diff feature (color or shape) of cue
dcfs_relF= np.array([[0, 1, 0, 1, 0, 1, 0, 1],
                     [1, 0, 1, 0, 0, 1, 0, 1], 
                     [0, 1, 0, 1, 1, 0, 1, 0],
                     [1, 0, 1, 0, 1, 0, 1, 0],
                     [0, 1, 0, 1, 0, 1, 0, 1],
                     [0, 1, 0, 1, 1, 0, 1, 0],
                     [1, 0, 1, 0, 0, 1, 0, 1],
                     [1, 0, 1, 0, 1, 0, 1, 0]])

# 0 means performing same task (face or scene) with this cue
# 1 means performing diff task (face or scene) with this cue
dcfs_tskR= np.array([[0, 1, 0, 1, 1, 0, 1, 0],
                     [1, 0, 1, 0, 0, 1, 0, 1], 
                     [0, 1, 0, 1, 1, 0, 1, 0],
                     [1, 0, 1, 0, 0, 1, 0, 1],
                     [1, 0, 1, 0, 0, 1, 0, 1],
                     [0, 1, 0, 1, 1, 0, 1, 0],
                     [1, 0, 1, 0, 0, 1, 0, 1],
                     [0, 1, 0, 1, 1, 0, 1, 0]])


# -- # DSFC
# 0 means looking at same feature (color or shape) of cue
# 1 means looking at diff feature (color or shape) of cue
dsfc_relF= np.array([[0, 0, 1, 1, 0, 0, 1, 1],
                    [0, 0, 1, 1, 1, 1, 0, 0], 
                    [1, 1, 0, 0, 0, 0, 1, 1],
                    [1, 1, 0, 0, 1, 1, 0, 0],
                    [0, 0, 1, 1, 0, 0, 1, 1],
                    [1, 1, 0, 0, 0, 0, 1, 1],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                    [1, 1, 0, 0, 1, 1, 0, 0]])

# 0 means performing same task (face or scene) with this cue
# 1 means performing diff task (face or scene) with this cue
dsfc_tskR= np.array([[0, 0, 1, 1, 0, 0, 1, 1],
                     [0, 0, 1, 1, 0, 0, 1, 1], 
                     [1, 1, 0, 0, 1, 1, 0, 0],
                     [1, 1, 0, 0, 1, 1, 0, 0],
                     [0, 0, 1, 1, 0, 0, 1, 1],
                     [0, 0, 1, 1, 0, 0, 1, 1],
                     [1, 1, 0, 0, 1, 1, 0, 0],
                     [1, 1, 0, 0, 1, 1, 0, 0]])


# ---------------------- Run Statistics ------------------------- #
# set up directories
data_dir="/data/backed_up/shared/ThalHi_MRI_2020/"
mnt_dir="/mnt/cifs/rdss/rdss_kahwang/ThalHi_data/MRI_data/"

# load numpy array from agg_resp.py
array_dir=data_dir+"RSA/ThalHi_SortedSubject_ROI_Coeff_Matrix_TEST.npy"
SUBxROIxCOEFF_mat=np.load(array_dir)

# load subject order ... should contain version info
sub_list_dir=data_dir+"RSA/ThalHi_SortedSubject_ROI_Coeff_Matrix_SubjectList.csv"
version_info=pd.read_csv(sub_list_dir)
#print(version_info.shape)

# set up regression variables
regressor_list=["Intercept","Context","Identity","Relevant_Feature","Task_Performed"]
stat_list=["_Beta","_T-stat","_P-value"]
print("\nGenerating column headers based on given regressors and stats to save out ... \n")
column_headers=[]
for rr in regressor_list:
    for ss in stat_list:
        column_headers.append((rr+ss))
print(column_headers,"\n")

# initialize lists for final data frame
Context_col=[]
Identity_col=[]
Relevant_Feature_col=[]
Task_Performed_col=[]
stats_data=np.zeros((400,len(regressor_list)*3))

# loop through ROIs ... reducing to subXcueXcue
for roi in range(SUBxROIxCOEFF_mat.shape[1]):
    print("\nCurrently on ROI",roi)
    SUBxCOEFF_mat=SUBxROIxCOEFF_mat[:,(roi-1),:,:] # pull out sub and coeff mat for current roi
    #print(SUBxCOEFF_mat.shape) # now 38x8x8
    #print(SUBxCOEFF_mat[1,:,:],"\n")
    #print(SUBxCOEFF_mat[2,:,:],"\n")
    # -- Set up empty vectors to fill as I loop through subjects
    Y_vec=np.zeros((SUBxCOEFF_mat.shape[0]*SUBxCOEFF_mat.shape[1]*SUBxCOEFF_mat.shape[2])) # set up empty vector to fill in as I loop through subjects and flatten the coeff matrix
    Intercept_vec=np.ones((SUBxCOEFF_mat.shape[0]*SUBxCOEFF_mat.shape[1]*SUBxCOEFF_mat.shape[2])) # set up vector of ones for the model intercept
    Context_vec=np.zeros((SUBxCOEFF_mat.shape[0]*SUBxCOEFF_mat.shape[1]*SUBxCOEFF_mat.shape[2])) # set up empty vector to fill in as I loop through subjects and flatten the coeff matrix
    Identity_vec=np.zeros((SUBxCOEFF_mat.shape[0]*SUBxCOEFF_mat.shape[1]*SUBxCOEFF_mat.shape[2])) # set up empty vector to fill in as I loop through subjects and flatten the coeff matrix
    RelFeat_vec=np.zeros((SUBxCOEFF_mat.shape[0]*SUBxCOEFF_mat.shape[1]*SUBxCOEFF_mat.shape[2])) # set up empty vector to fill in as I loop through subjects and flatten the coeff matrix
    TaskPerform_vec=np.zeros((SUBxCOEFF_mat.shape[0]*SUBxCOEFF_mat.shape[1]*SUBxCOEFF_mat.shape[2])) # set up empty vector to fill in as I loop through subjects and flatten the coeff matrix
    print("Y is a column vector of length: ", Y_vec.shape)
    for sub in range(SUBxCOEFF_mat.shape[0]):
        #print(SUBxCOEFF_mat[sub,:,:].flatten()[:]) # flatten the matrix
        #print(SUBxCOEFF_mat[sub,:,:].flatten().shape)
        #    flatten the coeff matrix
        coeff_vec=SUBxCOEFF_mat[sub,:,:].flatten()
        #print(np.where(Y_vec==0)[0][0])
        start_pt=np.where(Y_vec==0)[0][0]
        end_pt=np.where(Y_vec==0)[0][0]+len(coeff_vec)
        Y_vec[np.where(Y_vec==0)[0][0]:(np.where(Y_vec==0)[0][0]+len(coeff_vec))]=coeff_vec

        #    set up models based on what version the current sub did
        Context_vec[start_pt:end_pt]=data_Ctxt.flatten()
        Identity_vec[start_pt:end_pt]=data_Iden.flatten()
        if version_info["version"][sub]=="DCFS":
            #    not swapped version
            #print("not swapped:",version_info["version"][sub])
            RelFeat_vec[start_pt:end_pt]=dcfs_relF.flatten()
            TaskPerform_vec[start_pt:end_pt]=dcfs_tskR.flatten()
        else:
            #print("swapped:",version_info["version"][sub])
            RelFeat_vec[start_pt:end_pt]=dsfc_relF.flatten()
            TaskPerform_vec[start_pt:end_pt]=dsfc_tskR.flatten()
    # print(Y_vec)
    # print(Context_vec)
    # print(Identity_vec)
    # print(RelFeat_vec)
    # print(TaskPerform_vec)
    
    # -- create data frame for regressors
    data=np.array([Y_vec,Intercept_vec,Context_vec,Identity_vec,RelFeat_vec,TaskPerform_vec]).T
    df=pd.DataFrame(data=data, columns=["Y","Intercept","Context","Identity","Relevant_Feature","Task_Performed"])
    #print(df)

    # -- save out data frame for current ROI (usefull for double checking data input to model if needed)
    ROI_dataset_filename = "ROI_" + str(roi) + "_Dataset.csv"
    df.to_csv(os.path.join(data_dir,"RSA","ROI_Level_Datasets",ROI_dataset_filename))

    # -- run stats for current roi
    model = sm.OLS(df["Y"], df[regressor_list])
    #model = sm.OLS(df["Y"], df[["Context","Relevant_Feature","Task_Performed"]])
    res = model.fit()
    # print(res.summary())
    # print(res.params)
    # print(res.tvalues)
    # print(res.pvalues)

    # -- Pull out Betas, T-stats, and P-values
    #    Note, order is: [0]=Intercept [1]=Context, [2]=Identity, [3]=Relevant_Feature, [4]=Task_Performed
    cur_row_list=[]
    for cur_regressor in range(0,5):
        cur_row_list.append(res.params[cur_regressor])
        cur_row_list.append(res.tvalues[cur_regressor])
        cur_row_list.append(res.pvalues[cur_regressor])
    #print(cur_row_list)
    
    #    add stats for current roi to final data frame
    stats_data[(roi-1),:]=cur_row_list
    #print("\n\n",stats_data[:2,:])
    
# ---------------------- Save Statistics ------------------------- #
print("\nPreview of stats data frame...\n")
test_df=pd.DataFrame(data=stats_data, columns=column_headers)
test_df["ROI"]=range(1,401) # add ROI info
print(test_df)
print("Saving stats data frame to",os.path.join(data_dir,"RSA"))
test_df.to_csv(os.path.join(data_dir,"RSA","ThalHi_MRI_2020_RSA_Stats.csv"))


# -- NOW turn each column of the statistics data frame into a nifty image
#    load roi mask file
mask_file = data_dir+"ROIs/Schaefer400_2.5.nii.gz"
print("\n\nPulling mask file from ... ", mask_file)
mask=nib.load(mask_file)
mask_data=mask.get_fdata()
#print(mask_data.shape)

#    loop through ROIs and stats and create a nifti image with stats
for cur_stat in column_headers[3:]:
    mask_file = data_dir+"ROIs/Schaefer400_2.5.nii.gz"
    print("\n\nPulling mask file from ... ", mask_file)
    mask=nib.load(mask_file)
    mask_data=mask.get_fdata()
    print(mask_data.shape)
    for roi in test_df["ROI"]:
        mask_data[mask_data==float(roi)] = test_df.loc[test_df.ROI==roi, cur_stat]
    print("Min and Max for",cur_stat,"was",mask_data.min(),mask_data.max())
    mask_w_stats=nilearn.image.new_img_like(mask, mask_data)
    print("Mask dimensions still ", mask_data.shape)
    mask_filename=data_dir+"RSA/Statistics/"+cur_stat+".nii.gz"
    print("Saving",mask_filename)
    mask_w_stats.to_filename(mask_filename)
