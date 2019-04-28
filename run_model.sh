# model for 3dDeconvolve


cd /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/

lnum=$(cut -f13-18 sub-20190409_task-MB2_run-001_desc-confounds_regressors.tsv | wc -l)
cut -f13-18 sub-20190409_task-MB2_run-001_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb2compcor_run1.1d
lnum=$(cut -f13-18 sub-20190409_task-MB2_run-002_desc-confounds_regressors.tsv | wc -l)
cut -f13-18 sub-20190409_task-MB2_run-002_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb2compcor_run2.1d
cat mb2compcor_run1.1d mb2compcor_run2.1d > mb2compcor

lnum=$(cut -f13-18 sub-20190409_task-MB4_run-001_desc-confounds_regressors.tsv | wc -l)
cut -f13-18 sub-20190409_task-MB4_run-001_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb4compcor_run1.1d
lnum=$(cut -f13-18 sub-20190409_task-MB4_run-002_desc-confounds_regressors.tsv | wc -l)
cut -f13-18 sub-20190409_task-MB4_run-002_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb4compcor_run2.1d
cat mb4compcor_run1.1d mb4compcor_run2.1d > mb4compcor

lnum=$(cut -f13-18 sub-20190409_task-MB2_run-001_desc-confounds_regressors.tsv | wc -l)
cut -f26-31 sub-20190409_task-MB2_run-001_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb2moco_run1.1d
lnum=$(cut -f13-18 sub-20190409_task-MB2_run-002_desc-confounds_regressors.tsv | wc -l)
cut -f26-31 sub-20190409_task-MB2_run-002_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb2moco_run2.1d
cat mb2moco_run1.1d mb2moco_run2.1d > mb2moco

lnum=$(cut -f13-18 sub-20190409_task-MB4_run-001_desc-confounds_regressors.tsv | wc -l)
cut -f31-36 sub-20190409_task-MB4_run-001_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb4moco_run1.1d
lnum=$(cut -f13-18 sub-20190409_task-MB4_run-002_desc-confounds_regressors.tsv | wc -l)
cut -f31-36 sub-20190409_task-MB4_run-002_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb4moco_run2.1d
cat mb4moco_run1.1d mb4moco_run2.1d > mb4moco

lnum=$(cut -f13-18 sub-20190409_task-MB4_run-001_desc-confounds_regressors.tsv | wc -l)
cut -f3 sub-20190409_task-MB4_run-001_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb4gsr_run1.1d
lnum=$(cut -f13-18 sub-20190409_task-MB4_run-002_desc-confounds_regressors.tsv | wc -l)
cut -f3 sub-20190409_task-MB4_run-002_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb4gsr_run2.1d
cat mb4gsr_run1.1d mb4gsr_run2.1d > mb4gsr

lnum=$(cut -f13-18 sub-20190409_task-MB2_run-001_desc-confounds_regressors.tsv | wc -l)
cut -f3 sub-20190409_task-MB2_run-001_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb2gsr_run1.1d
lnum=$(cut -f13-18 sub-20190409_task-MB2_run-002_desc-confounds_regressors.tsv | wc -l)
cut -f3 sub-20190409_task-MB2_run-002_desc-confounds_regressors.tsv | tail -n $(($lnum-1)) > mb2gsr_run2.1d
cat mb2gsr_run1.1d mb2gsr_run2.1d > mb2gsr

rm -rf /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409
mkdir /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409

3dDeconvolve -input $(/bin/ls /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/sub-20190409_task-MB2_run-0*_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz | sort -V) \
-automask \
-polort A \
-num_stimts 2 \
-CENSORTR 1:0..4 \
-CENSORTR 2:0..4 \
-ortvec /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/mb2moco moco \
-ortvec /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/mb2compcor confounds \
-ortvec /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/mb2gsr gsr \
-stim_times 1 /data/backed_up/shared/ThalHi_MRI/ScanLogs/D001_MB2_Stay_stimtime.1D 'TENT(0, 15.75, 9)' -stim_label 1 Stay \
-stim_times 2 /data/backed_up/shared/ThalHi_MRI/ScanLogs/D001_MB2_Switch_stimtime.1D 'TENT(0, 15.75, 9)' -stim_label 2 Switch \
-iresp 1 /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB2_Stay_FIR.nii.gz \
-iresp 2 /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB2_Switch_FIR.nii.gz \
-gltsym 'SYM: +1*Switch[ -1*Stay[ ' -glt_label 1 MB2_Switch-Stay \
-gltsym 'SYM: +1*Stay[ ' -glt_label 2 MB2_Stay \
-gltsym 'SYM: +1*Switch[ ' -glt_label 3 MB2_Switch \
-rout \
-tout \
-bucket /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB2_FIRmodel_stats.nii.gz \
-GOFORIT 100 \
-noFDR \
-nocout \
-errts /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB2_FIR_errts.nii.gz \
-allzero_OK	-jobs 16

. /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB2_FIRmodel_stats.REML_cmd


3dDeconvolve -input $(/bin/ls /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/sub-20190409_task-MB4_run-0*_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz | sort -V) \
-automask \
-polort A \
-num_stimts 2 \
-CENSORTR 1:0..5 \
-CENSORTR 2:0..5 \
-ortvec /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/mb4moco moco \
-ortvec /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/mb4compcor confounds \
-ortvec /data/backed_up/shared/ThalHi_MRI/fmriprep/fmriprep/sub-20190409/func/mb4gsr gsr \
-stim_times 1 /data/backed_up/shared/ThalHi_MRI/ScanLogs/D001_MB4_Stay_stimtime.1D 'TENT(0, 15, 15)' -stim_label 1 Stay \
-stim_times 2 /data/backed_up/shared/ThalHi_MRI/ScanLogs/D001_MB4_Switch_stimtime.1D 'TENT(0, 15, 15)' -stim_label 2 Switch \
-iresp 1 /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB4_Stay_FIR.nii.gz \
-iresp 2 /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB4_Switch_FIR.nii.gz \
-gltsym 'SYM: +1*Switch -1*Stay ' -glt_label 1 MB4_Switch-Stay \
-gltsym 'SYM: +1*Stay ' -glt_label 2 MB4_Stay \
-gltsym 'SYM: +1*Switch ' -glt_label 3 MB4_Switch \
-rout \
-tout \
-bucket /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB4_FIRmodel_stats.nii.gz \
-GOFORIT 100 \
-noFDR \
-nocout \
-errts /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB4_FIR_errts.nii.gz \
-allzero_OK	-jobs 16

. /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409/MB4_FIRmodel_stats.REML_cmd

cd /data/backed_up/shared/ThalHi_MRI/Results/sub-20190409

ln -s /data/backed_up/shared/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii .





