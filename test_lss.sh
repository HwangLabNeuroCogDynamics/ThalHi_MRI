# run LSS on thalahi dataset.
#$ -N thalhi_lss
#$ -q SEASHORE
#$ -pe smp 3
#$ -t 1-51
#$ -ckpt user
#$ -o /Users/kahwang/sge_logs/
#$ -e /Users/kahwang/sge_logs/
/bin/echo Running on compute node: `hostname`.
/bin/echo Job: $JOB_ID
/bin/echo Task: $SGE_TASK_ID
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

subjects=(10001  10014  10027  10038  10068 10002  10017  10028  10040  10071 10003  10018  10029  10041  10073 10004  10019  10030  10042  10074 10005  10020  10031  10043  10076 10007  10021  10032  10054  10077 10008  10022  10033  10059  10080 10010  10023  10034  10060 10011  10024  10035  10063 10012  10025  10036  10064 10013  10026  10037  10066)
echo subjects: ${subjects[@]}
echo total_subjects=${#subjects[@]}
subject="${subjects[$SGE_TASK_ID-1]}"
echo "Starting 3dDeconvolve/LSS on $subject"

#lss
if [ ! -f /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz ]; then
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
    -mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
    -polort 3 \
    -ort /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
    -prefix /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
    -overwrite
fi

#3dDeconvolve stop at xmat, have to loop through 8 cue types....
singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times_IM 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times_IM 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times_IM 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times_IM 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times_IM 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times_IM 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times_IM 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1

singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
3dDeconvolve \
-input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
-mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
-polort A \
-concat '1D: 0 215 430 645 860 1075 1290 1505' \
-censor /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/censor.1D \
-ortvec /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/nuisance.1D \
-local_times \
-num_stimts 8 \
-stim_times 1 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcb.1D.txt 'SPMG' -stim_label 1 dcb \
-stim_times 2 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dcr.1D.txt 'SPMG' -stim_label 2 dcr \
-stim_times 3 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpb.1D.txt 'SPMG' -stim_label 3 dpb \
-stim_times 4 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/dpr.1D.txt 'SPMG' -stim_label 4 dpr \
-stim_times 5 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcb.1D.txt 'SPMG' -stim_label 5 fcb \
-stim_times 6 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fcr.1D.txt 'SPMG' -stim_label 6 fcr \
-stim_times 7 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpb.1D.txt 'SPMG' -stim_label 7 fpb \
-stim_times_IM 8 /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.1D.txt 'SPMG' -stim_label 8 fpr \
-x1D /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr_IM.xmat.1D \
-x1D_stop \
-allzero_OK \
-jobs 1


# then run lss. See:
# see https://www.sciencedirect.com/science/article/pii/S1053811911010081
# the outputs, with spmg option, is ordered as amp then deriv for each trial. 
if [ ! -f /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/fpr.LSS.nii.gz ]; then
    for cue in dcb dcr dpb dpr fcb fcr fpb fpr; do
    singularity run --cleanenv /Shared/lss_kahwang_hpc/opt/afni/afni.sif \
    3dLSS -input /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/lss_errts.nii.gz \
    -mask /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/combined_mask+tlrc \
    -matrix /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/${cue}_IM.xmat.1D \
    -prefix /Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve/sub-${subject}/${cue}.LSS.nii.gz \
    -overwrite
    done
fi
