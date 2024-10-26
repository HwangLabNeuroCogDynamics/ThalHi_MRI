# include RT regressor into the 3dDeconvolve model.
#$ -N FuRT
#$ -q SEASHORE
#$ -pe smp 4
#$ -t 1-73
#$ -tc 19
#$ -ckpt user
#$ -o /Users/kahwang/sge_logs/
#$ -e /Users/kahwang/sge_logs/
/bin/echo Running on compute node: `hostname`.
/bin/echo Job: $JOB_ID
/bin/echo Task: $SGE_TASK_ID
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`


subjects=(10001  10006  10011  10017  10022  10027  10032  10037  10042  10057  10062  10068  10074  10169  10175 10002  10007  10012  10018  10023  10028  10033  10038  10043  10058  10063  10069  10076  10170  10176 10003  10008  10013  10019  10024  10029  10034  10039  10044  10059  10064  10071  10077  10172  10179 10004  10009  10014  10020  10025  10030  10035  10040  10054  10060  10065  10072  10080  10173 10005  10010  10016  10021  10026  10031  10036  10041  10055  10061  10066  10073  10162  10174)
echo subjects: ${subjects[@]}
echo total_subjects=${#subjects[@]}
subject="${subjects[$SGE_TASK_ID-1]}"
echo "Starting 3dDeconvolve/LSS on $subject"
deconvolve_path="/Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/"

# get rid of nans
sed -i "s/nan//g" ${deconvolve_path}sub-${subject}/RT.1D
sed -i "s/nan//g" ${deconvolve_path}sub-${subject}/All_Trials.1D
sed -i "s/nan//g" ${deconvolve_path}sub-${subject}/EDS.acc.1D
sed -i "s/nan//g" ${deconvolve_path}sub-${subject}/IDS.acc.1D
sed -i "s/nan//g" ${deconvolve_path}sub-${subject}/Stay.acc.1D

# create amplitude modulated regressor
singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
1dMarry ${deconvolve_path}sub-${subject}/All_Trials.1D ${deconvolve_path}sub-${subject}/RT.1D > ${deconvolve_path}sub-${subject}/tmp.1D
cat ${deconvolve_path}sub-${subject}/tmp.1D | tail -n 8 > ${deconvolve_path}sub-${subject}/trial_RT.1D

# 3dDeconvolve
singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-2_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-3_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-4_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-5_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-6_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-7_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-8_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
-mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
-polort A \
-censor ${deconvolve_path}sub-${subject}/censor.1D \
-ortvec ${deconvolve_path}sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 4 \
-stim_times 1 ${deconvolve_path}sub-${subject}/EDS.acc.1D 'TENT(6, 20.4, 9)' -stim_label 1 EDS \
-stim_times 2 ${deconvolve_path}sub-${subject}/IDS.acc.1D 'TENT(6, 20.4, 9)' -stim_label 2 IDS \
-stim_times 3 ${deconvolve_path}sub-${subject}/Stay.acc.1D 'TENT(6, 20.4, 9)' -stim_label 3 Stay \
-stim_times_AM1 4 ${deconvolve_path}sub-${subject}/trial_RT.1D 'TENT(6, 20.4, 9)' -stim_label 4 RT \
-iresp 1 ${deconvolve_path}sub-${subject}/sub-${subject}_EDS_acc_rt_FIR_MNI.nii.gz \
-iresp 2 ${deconvolve_path}sub-${subject}/sub-${subject}_IDS_acc_rt_FIR_MNI.nii.gz \
-iresp 3 ${deconvolve_path}sub-${subject}/sub-${subject}_Stay_acc_rt_FIR_MNI.nii.gz \
-iresp 4 ${deconvolve_path}sub-${subject}/sub-${subject}_RT_acc_rt_FIR_MNI.nii.gz \
-num_glt 9 \
-gltsym 'SYM: +1*EDS - 1*Stay' -glt_label 1 EDS-Stay \
-gltsym 'SYM: +1*EDS' -glt_label 2 EDS \
-gltsym 'SYM: +1*IDS' -glt_label 3 IDS \
-gltsym 'SYM: +1*Stay' -glt_label 4 Stay \
-gltsym 'SYM: +1*EDS - 1*IDS' -glt_label 5 EDS-IDS \
-gltsym 'SYM: +1*IDS - 1*Stay' -glt_label 6 IDS-Stay \
-gltsym 'SYM: +1*EDS + 1*IDS + 1*Stay' -glt_label 7 All \
-gltsym 'SYM: +1*EDS + 1*IDS - 2*Stay' -glt_label 8 Switch \
-gltsym 'SYM: +1*RT' -glt_label 9 RT \
-rout \
-tout \
-bucket ${deconvolve_path}sub-${subject}/sub-${subject}_FIRmodel_MNI_stats_RT_tent_model.nii.gz \
-errts ${deconvolve_path}sub-${subject}/sub-${subject}_FIRmodel_RT_tent_model_errts.nii.gz \
-noFDR \
-nocout \
-jobs 8 \
-ok_1D_text \
-GOFORIT 99 

chmod 775 ${deconvolve_path}sub-${subject}/sub-${subject}_FIRmodel_MNI_stats_RT_tent_model.REML_cmd
singularity exec --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
${deconvolve_path}sub-${subject}/sub-${subject}_FIRmodel_MNI_stats_RT_tent_model.REML_cmd