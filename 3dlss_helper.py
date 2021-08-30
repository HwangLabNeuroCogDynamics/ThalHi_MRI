from thalpy import base
from thalpy.constants import paths, wildcards
from thalpy.analysis import motion, regressors

import pandas as pd
import os
import argparse
import warnings
import numpy as np
import glob as glob
import nibabel as nib


# Classes
class Stimfile:
    def __init__(self, task_name, sub_deconvolve_dir):
        self.name = task_name
        self.file = f"{task_name}.1D.txt"
        self.filepath = f"{sub_deconvolve_dir}{self.file}"
        self.runs = list()

    def write_file(self):
        with open(self.filepath, "w") as text_file:
            for run in self.runs:
                for event in run.timing_list:
                    text_file.write(f"{str(event)} ")
                text_file.write("\n")
        return


class Run:
    def __init__(self, run_number):
        self.number = run_number
        self.timing_list = list()


def generate_stimfiles(stimfile, run_files, timing_header, seperator="\t"):
    # add run timing data to stimfiles

    conditions_list = []
    for run_num, run_file in enumerate(run_files, start=1):
        # load event timing tsv files
        print(run_file)
        try:
            run_df = pd.read_csv(run_file, sep=seperator)
        except:
            stimfile.runs.append(Run(run_num))
            current_run = stimfile.runs[-1]
            if len(current_run.timing_list) == 0:
                current_run.timing_list.append("*")

        # add run to each stimfile
        stimfile.runs.append(Run(run_num))

        # append time to stimfile with associated task
        for row in zip(run_df["cue"], run_df[timing_header]):
            conditions_list.append([row[0], row[1], run_num])

            stimfile.runs[-1].timing_list.append(row[1])

        # insert * if no timing for run
        current_run = stimfile.runs[-1]
        if len(current_run.timing_list) == 0:
            current_run.timing_list.append("*")

    sub_df = pd.DataFrame(conditions_list, columns=["Cue", "Time", "Run"])
    return sub_df


# Create stimfile and generate csv for all conditions
# dir_tree = base.DirectoryTree("/data/backed_up/shared/ThalHi_MRI_2020/")
# subjects = base.get_subjects(dir_tree.deconvolve_dir, dir_tree)

# for sub in subjects:
#     run_files = base.get_ses_files(
#         sub,
#         "/mnt/cifs/rdss/rdss_kahwang/ThalHi_data/MRI_data/Behavioral_data/",
#         "*.csv",
#     )
#     print(sub.name)
#     stimfile = Stimfile("all", sub.deconvolve_dir)
#     sub_df = generate_stimfiles(
#         stimfile,
#         run_files,
#         "Time_Since_Run_Cue_Prez",
#         seperator=",",
#     )

#     stimfile.write_file()
#     sub_df.to_csv(sub.deconvolve_dir + "compiled_conditions.csv")

dir_tree = base.DirectoryTree("/data/backed_up/shared/ThalHi_MRI_2020/")
sub_10054 = base.Subject("10054", dir_tree, dir_tree.bids_dir)

cue_df = pd.read_csv(sub_10054.deconvolve_dir + "compiled_conditions.csv")

for name, group in cue_df.groupby(["Cue"]):
    os.chdir(sub_10054.deconvolve_dir)

    # sub_brick_list = group["SubBrick"].to_list()
    # print(
    #     f"3dTcat all.LSS+tlrc[{','.join(map(str, sub_brick_list))}] -output {name}.LSS.nii"
    # )
    # os.system(
    #     f"3dTcat all.LSS+tlrc[{','.join(map(str, sub_brick_list))}] -output {name}.LSS.nii"
    # )

    all_data = nib.load(f"{name}.LSS.nii").get_fdata()
    single_data = nib.load(f"{name}.LSS+tlrc.BRIK").get_fdata()
    all_data_2d = all_data.reshape(-1, all_data.shape[-1])
    single_data_2d = single_data.reshape(-1, single_data.shape[-1])

    correlations = []
    for trial in range(all_data_2d.shape[-1]):
        corr = np.corrcoef(all_data_2d[:, 0], single_data_2d[:, 0])
        correlations.append(corr[0][1])
    print(np.mean(correlations))
