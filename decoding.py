from thalpy import masks
from thalpy.constants.paths import SCRIPTS_DIR
import nibabel as nib
import numpy as np
import os
import sys
import pandas as pd
import pickle
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import statistics
from sklearn.decomposition import PCA


class SubjectLssTentData:
    def __init__(self, sub_deconvolve_dir, cues, tent_length=9, path="LSS_TENT.p"):
        self.cues = cues
        self.tent_length = tent_length
        self.sub_deconvolve_dir = sub_deconvolve_dir

        lss_img = nib.load(
            os.path.join(sub_deconvolve_dir, self.cues[0] + ".LSS+tlrc.BRIK")
        )
        lss_sample_data = lss_img.get_fdata()
        self.affine = lss_img.affine

        self.trial_df = pd.read_csv(
            os.path.join(self.sub_deconvolve_dir, "conditions.csv")
        )
        self.num_trials = len(self.trial_df.index)
        self.path = os.path.join(self.sub_deconvolve_dir, path)
        self.voxel_shape = lss_sample_data.shape[:3]

        self.__convert_lss_files()
        self.__avg_tent_matrix()
        self.save()

    @staticmethod
    def load(filepath):
        sys.path.append(os.path.join(SCRIPTS_DIR, "thalhi/"))
        return pickle.load(open(filepath, "rb"))

    def save(self, path=None):
        if path:
            self.path = path
        pickle.dump(self, open(self.path, "wb"), protocol=4)

    def __convert_lss_files(self):
        self.data = np.ones(
            [
                self.voxel_shape[0],
                self.voxel_shape[1],
                self.voxel_shape[2],
                self.tent_length * self.num_trials,
            ]
        )

        for cue, group in self.trial_df.groupby(["Cue"]):
            lss_data = nib.load(
                self.sub_deconvolve_dir + f"{cue}.LSS+tlrc.BRIK"
            ).get_fdata()

            if lss_data.shape[:3] != self.data.shape[:3]:
                raise Exception(
                    "Voxel shape between lss_data and trial matrix does not match."
                )

            for brik_idx, (trial_idx, trial) in enumerate(group.iterrows()):
                trial_tent_idx = trial_idx * self.tent_length
                brik_tent_idx = brik_idx * self.tent_length
                self.data[
                    :, :, :, trial_tent_idx: trial_tent_idx + self.tent_length
                ] = lss_data[:, :, :, brik_tent_idx: brik_tent_idx + self.tent_length]

    def __avg_tent_matrix(self):
        self.avg_data = np.ones(
            [
                self.voxel_shape[0],
                self.voxel_shape[1],
                self.voxel_shape[2],
                self.num_trials,
            ]
        )
        for trial in range(self.num_trials):
            self.avg_data[:, :, :, trial] = np.nanmean(
                self.data[:, :, :, trial * 9: (trial + 1) * 9], axis=-1
            )
            
    def remove_nan_trials(self):
        tents_to_remove = []
        trials_to_remove = []
        for trial in range(self.num_trials):
            if np.any(np.isnan(self.data[:, :, :, trial * 9: (trial + 1) * 9])):
                tents_to_remove = tents_to_remove + [ trial * 9 + tent_idx for tent_idx in range(9)]
                print(tents_to_remove)
                trials_to_remove.append(trial)
                
        self.data = np.delete(self.data, tents_to_remove, axis=3)
        self.trial_df = self.trial_df.drop(trials_to_remove)
        
        
    def mask_rois(self, mask_path):
        nii_img = nib.Nifti1Image(self.data, self.affine)
        self.rois = []
        num_rois = masks.masker_count(masks.binary_masker(mask_path))
        for i in range(num_rois):
            masker = masks.binary_masker(
                mask_path, img_math=f"img=={i + 1}"
            )
            self.rois.append(masker.fit_transform(nii_img))

    def eight_fold_runs(self, matrix):
        self.cue_array = label_binarize(
            self.trial_df.Cue.to_numpy(), self.cues)
        group_fold = GroupKFold(n_splits=self.trial_df.Run.nunique())
        splits = group_fold.split(
            matrix,
            self.cue_array,
            groups=self.trial_df.Run,
        )
        return splits

    def fit_model(self, classifier, output_path):
        num_rois = len(self.rois)
        accuracy = np.empty([num_rois, self.tent_length])
        probabilites = np.empty([num_rois, self.tent_length, len(self.cues)])
        lss_results = DecodingResult(accuracy, probabilites)

        for roi_idx in range(len(self.rois)):
            roi = self.rois[roi_idx]

            for tent_idx in range(self.tent_length):
                pca = PCA(0.95)
                PCA_components = pca.fit_transform(
                    roi[tent_idx:: self.tent_length, :])

                splits = self.eight_fold_runs(PCA_components)
                tent_accuracy = []
                tent_y_tests = []
                tent_y_scores = []

                for train_index, test_index in splits:
                    train_tent_index = [
                        roi_idx * self.tent_length + tent_idx for roi_idx in train_index
                    ]
                    test_tent_index = [
                        roi_idx * self.tent_length + tent_idx for roi_idx in test_index
                    ]

                    X_train, X_test = (
                        roi[train_tent_index, :],
                        roi[test_tent_index, :],
                    )
                    Y_train, Y_test = (
                        self.cue_array[train_index, :],
                        self.cue_array[test_index, :],
                    )

                    clf = classifier.fit(X_train, Y_train)

                    tent_y_scores.append(clf.decision_function(X_test))
                    tent_y_tests.append(Y_test)
                    # np.append(tent_probabilities, clf.predict_proba(X_test), axis=0)
                    tent_accuracy.append(clf.score(X_test, Y_test))

                lss_results.accuracy[roi_idx,
                                     tent_idx] = statistics.mean(tent_accuracy)
                # lss_results.probability[roi_idx, tent_idx] = np.mean(np.stack(tent_probabilities))
                lss_results.y_scores.append(tent_y_scores)
                lss_results.y_tests.append(tent_y_tests)

        lss_results.save(output_path)


class DecodingResult:
    def __init__(self, accuracy, probability):
        self.accuracy = accuracy
        self.probability = probability
        self.y_scores = []
        self.y_tests = []

    def roc_all(self):
        for roi_idx in range(self.y_scores.shape[0]):
            self.roc_roi(roi_idx)

    def roc_roi(self, roi_idx):
        roc(self.y_tests[roi_idx, :, :],
            self.y_scores[roi_idx, :, :], self.cues)


def roc(Y_test, y_score, classes):
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    colors = [
        "aqua",
        "darkorange",
        "cornflowerblue",
        "red",
        "yellow",
        "purple",
        "brown",
        "green",
    ]
    plt.figure()
    for i in range(len(classes)):
        fpr[i], tpr[i], _ = roc_curve(Y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
        # plot roc curves
        plt.plot(
            fpr[i],
            tpr[i],
            color=colors[i],
            lw=2,
            label=f"{classes[i]}" "".format(i, roc_auc[i]),
        )

    plt.plot([0, 1], [0, 1], "k--", lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(
        "Some extension of Receiver operating characteristic to multi-class")
    plt.legend(loc="lower right")
    plt.show()


def decode_cues(X, sphere_voxel_idxs, y, classifier, cues, runs, tent_length=9):
    accuracy = np.empty([tent_length])
    probabilites = np.empty([tent_length, len(cues)])
    lss_result = DecodingResult(accuracy, probabilites)

    for tent_idx in range(tent_length):
        splits = fold_runs(X[tent_idx:: tent_length, :], y, cues, runs)
        tent_accuracy = []
        tent_y_tests = []
        tent_y_scores = []

        for train_index, test_index in splits:
            X_train, X_test = (
                X[train_index, :],
                X[test_index, :],
            )
            Y_train, Y_test = (
                y[train_index],
                y[test_index],
            )

            clf = classifier.fit(X_train, Y_train)

            tent_y_scores.append(clf.decision_function(X_test))
            tent_y_tests.append(Y_test)
            # np.append(tent_probabilities, clf.predict_proba(X_test), axis=0)
            tent_accuracy.append(clf.score(X_test, Y_test))

        lss_result.accuracy[tent_idx] = statistics.mean(tent_accuracy)
        # lss_results.probability[roi_idx, tent_idx] = np.mean(np.stack(tent_probabilities))
        lss_result.y_scores.append(tent_y_scores)
        lss_result.y_tests.append(tent_y_tests)

    return lss_result


def fold_runs(X, y, cues, runs):
    y_binary = label_binarize(y, classes=cues)
    group_fold = GroupKFold(n_splits=len(np.unique(runs)))
    splits = group_fold.split(
        X,
        y_binary,
        groups=runs,
    )
    return splits
