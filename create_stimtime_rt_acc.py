## this is the version where we parse the RT and accurate trial stimulus timing information for 3dDeconvolve.
# Nov 29, 2023 Kai

import pandas as pd
import numpy as np

func_path = "/Shared/lss_kahwang_hpc/data/ThalHi/"
decon_path = "/Shared/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/"
df = pd.read_csv(func_path + "ThalHi_MRI_2020_RTs.csv")

subjects = df['sub'].unique().astype("int")

for sub in subjects:
    try:
        eds=pd.DataFrame()
        ids=pd.DataFrame()
        stay=pd.DataFrame()
        errors=pd.DataFrame()
        all_trials=pd.DataFrame()
        rt=pd.DataFrame()
        runs=[1,2,3,4,5,6,7,8]
        #print(sub)
        for r in runs:
            eds = eds.append(df.loc[(df['sub'] == sub) & (df['block'] == r) & (df['Trial_type'] == 'EDS') & (df['trial_Corr']==1)].Time_Since_Run_Cue_Prez.reset_index(drop=True))
            ids = ids.append(df.loc[(df['sub'] == sub) & (df['block'] == r) & (df['Trial_type'] == 'IDS') & (df['trial_Corr']==1)].Time_Since_Run_Cue_Prez.reset_index(drop=True))
            stay = stay.append(df.loc[(df['sub'] == sub) & (df['block'] == r) & (df['Trial_type'] == 'Stay') & (df['trial_Corr']==1)].Time_Since_Run_Cue_Prez.reset_index(drop=True))
            errors = errors.append(df.loc[(df['sub'] == sub) & (df['block'] == r) & (df['trial_Corr']==0)].Time_Since_Run_Cue_Prez.reset_index(drop=True)) # not sure we want to model this given the low number of error trials
            all_trials = all_trials.append(df.loc[(df['sub'] == sub) & (df['block'] == r) & (df['trial_Corr']==1)].Time_Since_Run_Cue_Prez.reset_index(drop=True))
            rt = rt.append(df.loc[(df['sub'] == sub) & (df['block'] == r) & (df['trial_Corr']==1)].rt.reset_index(drop=True))
            
        np.savetxt(decon_path + 'sub-{}/EDS.acc.1D'.format(sub),eds.to_numpy(),fmt='%1.4f') 
        np.savetxt(decon_path + 'sub-{}/IDS.acc.1D'.format(sub),ids.to_numpy(),fmt='%1.4f')
        np.savetxt(decon_path + 'sub-{}/Stay.acc.1D'.format(sub),stay.to_numpy(),fmt='%1.4f')
        np.savetxt(decon_path + 'sub-{}/Errors.1D'.format(sub),errors.to_numpy(),fmt='%1.4f')
        np.savetxt(decon_path + 'sub-{}/All_Trials.1D'.format(sub),all_trials.to_numpy(),fmt='%1.4f')
        np.savetxt(decon_path + 'sub-{}/RT.1D'.format(sub),rt.to_numpy(),fmt='%1.4f')
    except:
        print(sub)
    ## do this to get rid of nans sed -i "s/nan//g" RT.1D
    # sed -i "s/nan//g" All_Trials.1D


