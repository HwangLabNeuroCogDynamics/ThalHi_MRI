# import packages we'll need
#import matplotlib as plt
import nilearn
import nibabel as nib # requires some packages (listed below import list)
import numpy as np
import pandas as pd
import re
import glob
import os
import nilearn.masking

"""
This file is meant to load iresp output from 3dDeconvolve and ...
psedo code:
>> set up directories, lists, etc.
>> load mask file ...   nib.load("/data/backed_up/shared/ThalHi_MRI_2020/ROIs/Schaefer400_2.5.nii.gz)
>> for sub in sublist
>>    for roi in roi_list
>>       for cue in cuelist
>>          load iresp file ...  nib.load("sub-${sub}_${cue}_FIR_MNI.nii.gz")  
>>          apply mask to file ... nilearn.masking.apply_mask(data,mask)
>>       calculate coeff matrix & save in 2D array
>>       add coeff matrix to 3D array (roi by coeff matrix)
>>    add roi_by_coeff_matrix to 4D array (sub by roi by coeff matrix)   
>> save out 4D array and subject order info
"""
#    set up data directory
data_dir="/data/backed_up/shared/ThalHi_MRI_2020/"
mnt_dir="/mnt/cifs/rdss/rdss_kahwang/ThalHi_data/MRI_data/"

#    manually set cue list
#cue_list=["fpr","fpb","fcr","fcb","dpr","dpb","dcr","dcb"]
cue_list=["dcr","dcb","dpr","dpb","fcr","fpr","fcb","fpb"]

#    create subject list
sub_list=glob.glob("/data/backed_up/shared/ThalHi_MRI_2020/3dDeconvolve/*",recursive=False)
subid_list=[]
for sub in sub_list:
    sid=re.search("[0-9]{5}",sub)
    if sid:
        subid_list.append(sid.group(0))
#print(subid_list.sort())
unusable_subs = ['10006','10009',"10011","10015","10016","10029","10030","10034","10039","10055","10062","10065"]
subid_list2=set(subid_list).difference(unusable_subs)
print(subid_list2)

#    load mask file
mask_file = data_dir+"ROIs/Schaefer400_2.5.nii.gz"
print("\n\nPulling mask file from ... ", mask_file)
mask=nib.load(mask_file)
mask_data=mask.get_fdata()
print(mask_data.shape)

#    get set up to loop through subjects, ROIs, and cues
total_subs = len(subid_list2)
# # set up empty matrix for filling in ROI by cue by cue (400 X 8 X 8)
ROIxCOEFF_mat=np.zeros((400,len(cue_list),len(cue_list)))
# # set up empty matrix for filling in subject by ROI by cue by cue (num_subs X num_ROI X 8 X 8)
SUBxROIxCOEFF_mat=np.zeros((total_subs,400,len(cue_list),len(cue_list)))

print("Looping through", total_subs, "subjects")
counter = 0
#subid_list2=sorted(subid_list2)
for sub in sorted(subid_list2):
    counter+=1
    print("\n\nworking on subject ", sub, " ... ", counter, "/", total_subs)
    for roi in range(1,401):
        print("\ncurrently on ROI:", roi)
        #    pull out current roi
        mask_binary=np.where(mask_data==roi,1,0)
        #    create an empty matrix of size 8xROI (cuexROI)
        num_voxels=(mask_binary==1).sum()
        CUExVOXELS=np.zeros((len(cue_list),num_voxels)).T
        print("voxels in current ROI:", num_voxels, "  CUExVOXELS matrix dimensions:", CUExVOXELS.shape)
        for cue in cue_list:
            #    load iresp file
            iresp_file = data_dir+"3dDeconvolve/sub-"+sub+"/"+"sub-"+sub+"_"+cue+"_FIR_MNI.nii.gz"
            print("loading ... ", iresp_file)
            d=nib.load(iresp_file)
            #print("data size: ", d.get_fdata().shape, "\n")
            #print("applying mask ...")
            mask_binary_nif=nilearn.image.new_img_like(d, mask_binary)
            #print(type(mask_binary_nif))
            masked_data=nilearn.masking.apply_mask(d,mask_binary_nif) 
            print("masked data is a numpy array of size: ", masked_data.shape)
            #    masked_data now a numpy array (2D) ... just average over 2nd dimension for now
            masked_data_avg=np.mean(masked_data, axis=0)
            #print(masked_data_avg.shape)
            #    once i get my vector add this vector to my numpy matrix of size cuexROI
            CUExVOXELS[:,(cue_list.index(cue))]=masked_data_avg
            #print(CUExVOXELS)
        
        #    Now that the entire Cue by Voxel matrix has been filled in, 
        #    calculate correlation coefficients
        Coeff_mat=np.corrcoef(CUExVOXELS, rowvar=False)
        print("\n Coefficient matrix calculated for ROI", roi)
        #print(Coeff_mat)
        print(Coeff_mat.shape)
        #    put my 8x8 corr coef matrix inside a 400x8x8 matrix (all_ROIxcuexcue)
        # if (roi==1):
        #     #    if the first sub, initiate the matrix
        #     ROIxCOEFF_mat=Coeff_mat[None]
        # else:
        #     #    if not first sub, add to existing matrix
        #     ROIxCOEFF_mat=np.vstack((ROIxCOEFF_mat,Coeff_mat[None]))
        ROIxCOEFF_mat[(roi-1),:,:]=Coeff_mat
        #print("\n ROI by Coefficient matrix dimensions:",ROIxCOEFF_mat.shape,"\n")
        #print(ROIxCOEFF_mat)

    #    repeat ROIxCOEFF_mat merge but for subj
    # if (counter==1):
    #     #    if the first sub, initiate the matrix
    #     SUBxROIxCOEFF_mat=ROIxCOEFF_mat[None]
    # else:
    #     #    if not first sub, add to existing matrix
    #     SUBxROIxCOEFF_mat=np.vstack((SUBxROIxCOEFF_mat,ROIxCOEFF_mat[None]))
    SUBxROIxCOEFF_mat[(counter-1),:,:,:]=ROIxCOEFF_mat
    #print("\n\n Subject by ROI by Coefficient matrix dimensions:",SUBxROIxCOEFF_mat.shape)
    #print(SUBxROIxCOEFF_mat)

# # # # SAVE OUT FINAL MATRIX
data_dir+="RSA/"
os.chdir(data_dir) 
#np.save("ThalHi_SortedSubject_ROI_Coeff_Matrix_TEST.npy",SUBxROIxCOEFF_mat)

# # # # SAVE OUT SUBJECT LIST
sub_df=pd.DataFrame(sorted(subid_list2), columns=["sub"]) # make sub list into data frame
sub_df["sub"]=sub_df["sub"].astype("float")
version_dir=mnt_dir+"ThalHi_MRI_2020_RTs.csv"
version_info=pd.read_csv(version_dir)
version_info["version"][version_info["version"].isna()]="DCFS"
version_info["version"][version_info["version"]=="Donut=Shape,Filled=Color"]="DSFC"
version_info=version_info[["sub","version"]].drop_duplicates()
#    merge the data frames
version_info=sub_df.merge(version_info)
#    save the subject list
print("subject list saved to csv file\n",sorted(subid_list2))
version_info.to_csv("ThalHi_SortedSubject_ROI_Coeff_Matrix_SubjectList.csv")
#    also print the order of cues
print("cue order: ",cue_list)