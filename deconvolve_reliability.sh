#subject=10012
subjects=(10001 10002 10003 10004 10005 \
10008 10010 10012 10013 10014 \
10017 10018 10019 10020 10022 \
10023 10024 10025 10027 10028 \
10031 10032 10033 10034 10035 \
10036 10037 10038 10039 10040 \
10041 10042 10043 10044 10054 \
10057 10058 10059 10060 10063 \
10064 10066 10068 10069 10071 \
10072 10073 10074 10076 10077 \
10080 10162 10169 10170 10173 \
10174 10175 10176 10179)

for subject in ${subjects[@]}
do 
outputpath="/home/kahwang/argon/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}"
testpath="/home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/"
mkdir ${testpath}
data="/home/kahwang/argon/data/ThalHi/fmriprep/sub-${subject}/func"

# do whole brain signal regression
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz > /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-2*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz >> /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-3*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz >> /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-4*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz >> /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
1dcat ${testpath}/wb.1D ${outputpath}/nuisance_r1_r4.1D > ${testpath}/nuisance_r1_r4.1D

3dDeconvolve -input ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
${data}/sub-${subject}_task-ThalHi_run-2*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
${data}/sub-${subject}_task-ThalHi_run-3*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
${data}/sub-${subject}_task-ThalHi_run-4*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
-mask ${outputpath}/combined_mask+tlrc.BRIK \
-polort A \
-censor ${outputpath}/censor_r1_r4.1D \
-ortvec ${testpath}/nuisance_r1_r4.1D \
-local_times \
-num_stimts 3 \
-stim_times 1 ${outputpath}/EDS_r1_r4.1D.txt 'SPMG2' -stim_label 1 EDS \
-stim_times 2 ${outputpath}/IDS_r1_r4.1D.txt 'SPMG2' -stim_label 2 IDS \
-stim_times 3 ${outputpath}/Stay_r1_r4.1D.txt 'SPMG2' -stim_label 3 Stay \
-num_glt 8 \
-gltsym 'SYM: +1*EDS - 1*Stay' -glt_label 1 EDS-Stay \
-gltsym 'SYM: +1*EDS' -glt_label 2 EDS \
-gltsym 'SYM: +1*IDS' -glt_label 3 IDS \
-gltsym 'SYM: +1*Stay' -glt_label 4 Stay \
-gltsym 'SYM: +1*EDS - 1*IDS' -glt_label 5 EDS-IDS \
-gltsym 'SYM: +1*IDS - 1*Stay' -glt_label 6 IDS-Stay \
-gltsym 'SYM: +1*EDS + 1*IDS + 1*Stay' -glt_label 7 All \
-gltsym 'SYM: +1*EDS + 1*IDS - 2*Stay' -glt_label 8 Switch \
-rout \
-tout \
-bucket ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r1_r4.nii.gz \
-noFDR \
-jobs 16 \
-ok_1D_text

3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-5*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz > /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-6*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz >> /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-7*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz >> /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
3dmaskave -mask ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz[0] -quiet ${data}/sub-${subject}_task-ThalHi_run-8*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz >> /home/kahwang/argon/data/ThalHi/reliability_test/sub-${subject}/wb.1D
1dcat ${testpath}/wb.1D ${outputpath}//nuisance_r5_r8.1D > ${testpath}/nuisance_r5_r8.1D

3dDeconvolve -input ${data}/sub-${subject}_task-ThalHi_run-5*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
${data}/sub-${subject}_task-ThalHi_run-6*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
${data}/sub-${subject}_task-ThalHi_run-7*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
${data}/sub-${subject}_task-ThalHi_run-8*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz \
-mask ${outputpath}/combined_mask+tlrc.BRIK \
-polort A \
-censor ${outputpath}/censor_r5_r8.1D \
-ortvec ${testpath}/nuisance_r5_r8.1D \
-local_times \
-num_stimts 3 \
-stim_times 1 ${outputpath}/EDS_r5_r8.1D.txt 'SPMG2' -stim_label 1 EDS \
-stim_times 2 ${outputpath}/IDS_r5_r8.1D.txt 'SPMG2' -stim_label 2 IDS \
-stim_times 3 ${outputpath}/Stay_r5_r8.1D.txt 'SPMG2' -stim_label 3 Stay \
-num_glt 8 \
-gltsym 'SYM: +1*EDS - 1*Stay' -glt_label 1 EDS-Stay \
-gltsym 'SYM: +1*EDS' -glt_label 2 EDS \
-gltsym 'SYM: +1*IDS' -glt_label 3 IDS \
-gltsym 'SYM: +1*Stay' -glt_label 4 Stay \
-gltsym 'SYM: +1*EDS - 1*IDS' -glt_label 5 EDS-IDS \
-gltsym 'SYM: +1*IDS - 1*Stay' -glt_label 6 IDS-Stay \
-gltsym 'SYM: +1*EDS + 1*IDS + 1*Stay' -glt_label 7 All \
-gltsym 'SYM: +1*EDS + 1*IDS - 2*Stay' -glt_label 8 Switch \
-rout \
-tout \
-bucket ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r5_r8.nii.gz \
-noFDR \
-jobs 16 \
-ok_1D_text

# 3dresample -master ${data}/sub-${subject}_task-ThalHi_run-1*space-MNI152NLin2009cAsym_desc-preproc_bold*.nii.gz -inset /home/kahwang/argon/ROIs/mni_atlas/MNI_thalamus_2mm.nii.gz \
# -prefix /home/kahwang/argon/ROIs/thmask.nii.gz

3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r1_r4.nii.gz[2] > ${testpath}/EDS_1-4_beta.1D
3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r5_r8.nii.gz[2] > ${testpath}/EDS_5-8_beta.1D

3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r1_r4.nii.gz[7] > ${testpath}/IDS_1-4_beta.1D
3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r5_r8.nii.gz[7] > ${testpath}/IDS_5-8_beta.1D

3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r1_r4.nii.gz[12] > ${testpath}/Stay_1-4_beta.1D
3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r5_r8.nii.gz[12] > ${testpath}/Stay_5-8_beta.1D

3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r1_r4.nii.gz[35] > ${testpath}/All_1-4_beta.1D
3dmaskdump -mask /home/kahwang/argon/ROIs/thmask.nii.gz -noijk ${testpath}/sub-${subject}_SPMGmodel_MNI_stats_task_switch_r5_r8.nii.gz[35] > ${testpath}/All_5-8_beta.1D

done