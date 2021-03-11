import json
import common
import glob as glob
import sys
import basic_settings as bs
import os
import stat

sub_args = sys.argv[1:]
dataset_dir = '/data/backed_up/shared/ThalHi_MRI_2020/'
dir_tree = common.DirectoryTree(dataset_dir)
subject = common.subargs_to_subjects(sub_args, dir_tree)[0]
print(subject.name)
print(subject.bids_dir)

fieldmap_jsons = glob.glob(subject.bids_dir + bs.FMAP_DIR + '*.json')
runs = sorted(glob.glob(subject.bids_dir + bs.FUNC_DIR + '*.nii.gz'))
runs = [run.split(f'sub-{subject.name}/')[1] for run in runs]
print(runs)
for fmap in fieldmap_jsons:
    os.chmod(fmap, 0o755)

    with open(fmap, 'r') as jsonfile:
        json_content = json.load(jsonfile)

    if common.parse_run_from_file(fmap) == '001':
        intended_runs = runs[:4]
    elif common.parse_run_from_file(fmap) == '002':
        intended_runs = runs[4:]
    else:
        print('Run number must be 001 or 002')
        continue

    json_content['IntendedFor'] = intended_runs
    # correct json file intended_runs

    with open(fmap, 'w') as jsonfile:
        # you decide the indentation level
        json.dump(json_content, jsonfile, indent=4)
