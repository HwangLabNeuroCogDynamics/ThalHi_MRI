from nilearn.image import resample_to_img
from thalpy import base
from sklearn.linear_model import RidgeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.svm import SVC
import sys
from thalhi.decoding import SubjectLssTentData, decode_cues
from sklearn.multiclass import OneVsRestClassifier
import os
import argparse
from thalpy.decoding import searchlight
from thalpy import masks
import nibabel as nib
import pickle
from nilearn.datasets import load_mni152_template
import numpy as np


def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Decode subject",
        usage="[Subject] [Classifier] [OPTIONS] ... ",
    )
    parser.add_argument("subject", help="id of subject ie 10001")
    parser.add_argument(
        "--classifier", help="Classifier to use: Options LDA, SVM, Ridge")
    parser.add_argument("--searchlight",
                        help="run searchlight, default is false",
                        default=False, action="store_true")
    parser.add_argument("--cores",
                        help="number of cores",
                        default=4)
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

    dir_tree = base.DirectoryTree(THAL_HI_DIR)
    subjects = base.get_subjects(dir_tree.deconvolve_dir, dir_tree)

    subject = next(sub for sub in subjects if sub.name == args.subject)

    if args.classifier == "LDA":
        clf = OneVsRestClassifier(LinearDiscriminantAnalysis())
    elif args.classifier == "SVM":
        clf = OneVsRestClassifier(SVC())
    elif args.classifier == "Ridge":
        clf = OneVsRestClassifier(RidgeClassifier())

    os.chdir(subject.deconvolve_dir)

    if os.path.exists("LSS_TENT.p"):
        subject_lss_data = SubjectLssTentData.load("LSS_TENT.p")
        subject_lss_data.remove_nan_trials()
    else:
        print("Converting LSS files")
        subject_lss_data = SubjectLssTentData(subject.deconvolve_dir, CUES)
        if args.classifier:
            subject_lss_data.mask_rois(masks.SCHAEFER_400_7N_PATH)
        subject_lss_data.save()
        exit()

    if args.searchlight:
        # coritcal_mask = nib.load(masks.CORITCAL_BINARY_PATH)
        frontal_mask = nib.load(
            "/Shared/lss_kahwang_hpc/ROIs/MNI_frontal.nii.gz")

        img = nib.load(
            THAL_HI_DIR + "fmriprep/sub-10001/func/sub-10001_task-ThalHi_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")
        imgs = nib.Nifti1Image(subject_lss_data.data,
                               affine=img.affine, header=img.header)

        # resample imgs and mask to 3x3x3
        template = load_mni152_template(resolution=3)
        frontal_mask = resample_to_img(
            frontal_mask, template, interpolation='nearest')
        resampled_imgs = resample_to_img(imgs, template)

        sl_obj = searchlight.SearchLight(
            frontal_mask, decode_cues,
            [subject_lss_data.trial_df.Cue.to_numpy(), clf, CUES, subject_lss_data.trial_df.Run.to_numpy()], verbose=1, radius=10., n_jobs=int(args.cores))
        sl_obj.run(resampled_imgs)
        print(sl_obj.output)
        pickle.dump(sl_obj, open(os.path.join(
            dir_tree.deconvolve_dir, "searchlight_lda.p"), "wb"), protocol=4)
    elif args.classifier:
        subject_lss_data.fit_model(clf, f"LSS_RESULTS_{args.classifier}.p")


if __name__ == "__main__":
    main()
