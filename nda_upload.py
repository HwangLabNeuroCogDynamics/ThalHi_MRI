import os
import sys
import json
import hashlib
import argparse
from thalpy import base
from thalpy.constants import paths
import pandas as pd
import shutil
import msoffcrypto
import io
import pandas as pd
from dateutil.relativedelta import relativedelta


def load_sub_dict_from_mapping(mapping_path):
    pass


def pull_info_from_xlsx(subjects, session_date_label):
    password = input("Enter the password for the Subject Excel Sheet: ")

    file = msoffcrypto.OfficeFile(
        open("/mnt/cifs/rdss/rdss_kahwang/subject_info.xlsx", "rb")
    )

    file.load_key(password=password)

    decrypted = io.BytesIO()
    file.decrypt(decrypted)

    df = pd.read_excel(decrypted, sheet_name="Subject Info Sheet").set_index(
        "Subject #"
    )
    for sub in subjects:
        row = df.loc[int(sub.name)]
        sub.sex = row["Sex"]
        sub.interview_date = row[session_date_label]
        time_difference = relativedelta(sub.interview_date, row["DOB (age)"])
        sub.age = time_difference.years + time_difference.months / 12
        print(sub.age)


def create_participants_tsv(bids_dir, subjects):
    tsv_df = pd.DataFrame()
    for sub in subjects:
        sub_row = pd.DataFrame(
            [["sub-" + sub.name, sub.age, sub.sex, sub.interview_date]],
            columns=["participant_id", "interview_age", "sex", "age"],
        )
        tsv_df = tsv_df.append(sub_row)
    tsv_df.to_csv(os.path.join(bids_dir, paths.PARTICIPANTS_TSV), sep="\t", index=False)


if __name__ == "__main__":
    dir_tree = base.DirectoryTree(
        "/data/backed_up/shared/ThalHi_MRI_2020/",
        bids_dir="/data/backed_up/shared/ThalHi_MRI_2020/nda_bids",
    )
    os.chdir(dir_tree.bids_dir)
    sub_dict = load_sub_dict_from_mapping(dir_tree.bids_dir + "mapping.txt")

    # copy subjects we need from BIDS to nda_bids
    os.chdir("/data/backed_up/shared/ThalHi_MRI_2020/BIDS")
    for subject in sub_dict.keys():
        shutil.copytree(
            f"sub-{subject}/", dir_tree.bids_dir + f"sub-{subject}/", dirs_exist_ok=True
        )

    subjects = base.get_subjects(dir_tree.bids_dir, dir_tree)
    pull_info_from_xlsx(subjects, "session date")
    create_participants_tsv(dir_tree.bids_dir, subjects)
