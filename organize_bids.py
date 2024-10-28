# reorganize bids folder to deal with openneuro upload
import shutil
import os

df = pd.read_csv("/home/kahwang/argon/data/ThalHi/Version_Info.csv", sep=",")

# copy included subjects
for sub in df['sub']:
    source = '/mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/BIDS/sub-%s' %sub
    target = '/home/kahwang/bsh/OpenNeuro/ThalHi/BIDS/sub-%s' %sub
    shutil.copytree(source, target)

# now do the deface procedure
import os
from joblib import Parallel, delayed

def deface(sub):
    input_file = f"/home/kahwang/bsh/OpenNeuro/ThalHi/BIDS/sub-{sub}/anat/sub-{sub}_T1w.nii.gz"
    output_file = f"/home/kahwang/bsh/OpenNeuro/ThalHi/BIDS/sub-{sub}/anat/sub-{sub}_T1w.nii.gz"
    
    # Construct the command
    command = f"pydeface {input_file} --outfile {output_file} --force"
    
    # Execute the command
    os.system(command)

# Number of parallel jobs
num_jobs = 36

# Use joblib to run the deface function in parallel
Parallel(n_jobs=num_jobs)(delayed(deface)(sub) for sub in df['sub'])

#now organize tsv files
bdf = pd.read_csv("/home/kahwang/argon/data/ThalHi/ThalHi_MRI_2020_RTs.csv")
for sub in df['sub']:
    for r in np.arange(1,9):
        tdf = pd.DataFrame()
        rdf = bdf.loc[(bdf['sub']==sub) & (bdf['block']==r)].reset_index()
        tdf['onset'] = rdf['Time_Since_Run_Cue_Prez'].values
        tdf['duration'] = 3
        tdf['trial_type'] = rdf['Trial_type'].values
        tdf['response_time'] = rdf['rt'].values
        tdf['texture'] = rdf['Texture'].values
        tdf['color'] = rdf['Color'].values
        tdf['shape'] = rdf['Shape'].values
        tdf['task'] = rdf['Task'].values
        tdf['response'] = rdf['Subject_Respo'].values
        tdf['probe_pic'] = rdf['pic'].values
        tdf.fillna("n/a", inplace=True)
        fn = f"/data/backed_up/shared/OpenNeuro/ThalHi/BIDS/sub-{sub}/func/sub-{sub}_task-ThalHi_run-00{r}_events.tsv"
        tdf.to_csv(fn, sep='\t', index=False)
