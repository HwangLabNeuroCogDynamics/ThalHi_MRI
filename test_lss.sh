# run LSS on thalahi dataset.
#$ -N thalhi_lss
#$ -q SEASHORE
#$ -pe smp 4
#$ -t 1-73
#$ -tc 20
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

#lss
rm ${deconvolve_path}sub-${subject}/lss_errts.nii.gz

#if [ ! -f ${deconvolve_path}sub-${subject}/lss_errts.nii.gz ]; then
singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dTproject -input \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-2_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-3_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-4_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-5_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-6_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-7_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
/Shared/lss_kahwang_hpc/data/ThalHi/fmriprep/sub-${subject}/func/sub-${subject}_task-ThalHi_run-8_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz \
-mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
-polort 4 \
-ort ${deconvolve_path}sub-${subject}/nuisance.1D \
-prefix ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
-censor ${deconvolve_path}sub-${subject}/censor.1D \
-cenmode ZERO \
-overwrite
#fi


singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve -input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/combined_mask+tlrc.BRIK \
-polort 0 \
-concat '1D: 0 216 432 648 864 1080 1296 1512' \
-num_stimts 2 -local_times \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/censor.1D \
-stim_times_IM 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/cue_stimtimes.1D SPMG -stim_label 1 cue \
-stim_times_AM1 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/trial_RT.1D SPMG -stim_label 2 RT \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/cue_rtcov_IM.xmat.1D -x1D_stop -allzero_OK -jobs 4

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dLSS -input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/combined_mask+tlrc.BRIK \
-matrix /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/cue_rtcov_IM.xmat.1D \
-prefix /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${subject}/cue_rtcov_IM.LSS.nii.gz



### graveyard backup of old approach
# #3dDeconvolve stop at xmat, have to loop through 8 cue types....
# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times_IM 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/dcb_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times_IM 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/dcr_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times_IM 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/dpb_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times_IM 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/dpr_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times_IM 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/fcb_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times_IM 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/fcr_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times_IM 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/fpb_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1

# singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
# 3dDeconvolve \
# -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
# -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
# -polort 0 \
# -concat '1D: 0 216 432 648 864 1080 1296 1512' \
# -local_times \
# -num_stimts 8 \
# -stim_times 1 ${deconvolve_path}sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
# -stim_times 2 ${deconvolve_path}sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
# -stim_times 3 ${deconvolve_path}sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
# -stim_times 4 ${deconvolve_path}sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
# -stim_times 5 ${deconvolve_path}sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
# -stim_times 6 ${deconvolve_path}sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
# -stim_times 7 ${deconvolve_path}sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
# -stim_times_IM 8 ${deconvolve_path}sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
# -x1D ${deconvolve_path}sub-${subject}/fpr_IM.xmat.1D \
# -x1D_stop \
# -allzero_OK \
# -jobs 1


# # then run lss. See:
# # see https://www.sciencedirect.com/science/article/pii/S1053811911010081
# # the outputs, with spmg option, is ordered as amp then deriv for each trial. 
# if [ ! -f ${deconvolve_path}sub-${subject}/fpr.LSS.nii.gz ]; then
#     for cue in dcb dcr dpb dpr fcb fcr fpb fpr; do
#     singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
#     3dLSS -input ${deconvolve_path}sub-${subject}/lss_errts.nii.gz \
#     -mask ${deconvolve_path}sub-${subject}/combined_mask+tlrc \
#     -matrix ${deconvolve_path}sub-${subject}/${cue}_IM.xmat.1D \
#     -prefix ${deconvolve_path}sub-${subject}/${cue}.LSS.nii.gz \
#     -overwrite
#     done
# fi
