import fileinput
import csv
import sys
import pandas as pd
import glob
import numpy as np
import os.path

#parse block order
#os.chdir('/home/despoB/kaihwang/TRSE/TTD')
#Subjects = glob.glob('5*')
#these are old FTTD subjects
#[601, 602, 603, 605]
#os.chdir('/home/despoB/kaihwang/TRSE/TTD/ScanLogs')


def write_stimtime(filepath, inputvec):
	''' short hand function to write AFNI style stimtime'''
	#if os.path.isfile(filepath) is False:
	f = open(filepath, 'w')
	for val in inputvec[0]:
		if val =='*':
			f.write(val + '\n')
		else:
			# this is to dealt with some weird formating issue
			#f.write(np.array2string(val))
			#print np.array2string(np.around(val,4))
			f.write(np.array2string(np.around(val,4)).replace('\n','')[3:-1] + '\n')
	f.close()



def parse_stim(s, fns, ntrials_per_run, num_runs, postfix):
	''' parse stim log from pscyhtoolbox script
	input: s, subjects
		ROIs,
		ntrials_per_rn
		num_runs 
	'''
	
	os.chdir('/data/backed_up/shared/ThalHi_MRI/ScanLogs')
	num_runs = int(num_runs)
	print "parsing stimulus timing for subject %s" %(s)	
	timing_logs = sorted(fns) #important to sort since glob seems to randomly order files.

	#to check the order is correct, in ascending order
	print('check timing logs are in ascending order:')
	print(timing_logs) 

	#load timing logs, concat depend on number of sessions
	df = pd.read_csv(timing_logs[0]).dropna()
	for i in np.arange(1, len(timing_logs), 1):
		df = df.append(pd.read_csv(timing_logs[i]).dropna())
	
	# create new column variable of "run number" for every trial (48 trials per run). 12 blocks totoal 
	df['Run'] = np.repeat(np.arange(1, num_runs+1), ntrials_per_run)


	##### irrelevant stuff!
	# extract the order of each block condition and save to a text file
	# run_order = df.groupby(['Run', 'Condition', 'MotorMapping']).sum().reset_index()[['Run', 'Condition', 'MotorMapping']]
	# fn = '/home/despoB/TRSEPPI/TTD/ScanLogs/%s_%s_run_order' %(s, site)
	# if os.path.isfile(fn) is False:			
	# 	run_order[['Condition','MotorMapping']].to_csv(fn, index=None, header=None, sep = '\t')

	# #write out target+distractor (TD) runs
	# TD_runs = run_order[(run_order['Condition'] == 'HF') | (run_order['Condition'] == 'FH')]['Run'].values	
	# fn = '/home/despoB/TRSEPPI/TTD/ScanLogs/%s_TD_runs' %s
	# if os.path.isfile(fn) is False:			
	# 	np.savetxt(fn, TD_runs, fmt='%2d')

	# #write out 2bk runs
	# TD_runs = run_order[(run_order['Condition'] == 'H2') | (run_order['Condition'] == 'F2')]['Run'].values	
	# fn = '/home/despoB/TRSEPPI/TTD/ScanLogs/%s_2B_runs' %s
	# if os.path.isfile(fn) is False:			
	# 	np.savetxt(fn, TD_runs, fmt='%2d')	

	# #write out passive runs	
	# To_runs = run_order[(run_order['Condition'] == 'Hp') | (run_order['Condition'] == 'Fp')]['Run'].values	
	# fn = '/home/despoB/TRSEPPI/TTD/ScanLogs/%s_P_runs' %s
	# if os.path.isfile(fn) is False:			
	# 	np.savetxt(fn, To_runs, fmt='%2d')


	#create FIR timing for every condition, extract trial timing for the first trial of every block
	Stay_stimtime = [['*']*num_runs] #stimtime format, * for runs with no events of this stimulus class
	IDS_stimtime = [['*']*num_runs] #create one of this for every condition
	EDS_stimtime = [['*']*num_runs]
	Switch_stimtime = [['*']*num_runs]

	for i, run in enumerate(np.arange(1,num_runs+1)):  #loop through 12 runs
		run_df = df[df['Run']==run].reset_index() #"view" of block we are curreint sorting through
		Stay_run_trials = [] #empty vector to store trial time info for the current block
		IDS_run_trials = []
		EDS_run_trials = []
		Switch_run_trials = []


		for tr in np.arange(0,ntrials_per_run):  #this is to loop through trials
			if run_df.loc[tr,'Trial_type'] in ('Stay'):
				Stay_run_trials.append(run_df.loc[tr,'Time_Since_Run_Cue_Prez']) #if a match of condition, append trial timing						
			
			if run_df.loc[tr,'Trial_type'] in ('IDS'):
				IDS_run_trials.append(run_df.loc[tr,'Time_Since_Run_Cue_Prez']) 
				
			if run_df.loc[tr,'Trial_type'] in ('EDS'):
				EDS_run_trials.append(run_df.loc[tr,'Time_Since_Run_Cue_Prez']) 

			if run_df.loc[tr,'Trial_type'] in (['EDS', 'IDS']):
				Switch_run_trials.append(run_df.loc[tr,'Time_Since_Run_Cue_Prez']) 				

		if any(Stay_run_trials):
			Stay_stimtime[0][i] = Stay_run_trials	#put trial timing of each block into the stimtime array	
		if any(IDS_run_trials):
			IDS_stimtime[0][i] = IDS_run_trials			
		if any(EDS_run_trials):
			EDS_stimtime[0][i] = EDS_run_trials	
		if any(Switch_run_trials):
			Switch_stimtime[0][i] = Switch_run_trials					


	#write out stimtime array to text file. 	
	#postfix='MB2'	
	fn = '/data/backed_up/shared/ThalHi_MRI/ScanLogs/%s_%s_Stay_stimtime.1D' %(s, postfix)
	write_stimtime(fn, Stay_stimtime)		
		
	fn = '/data/backed_up/shared/ThalHi_MRI/ScanLogs/%s_%s_IDS_stimtime.1D' %(s, postfix)
	write_stimtime(fn, IDS_stimtime)	

	fn = '/data/backed_up/shared/ThalHi_MRI/ScanLogs/%s_%s_EDS_stimtime.1D' %(s, postfix)
	write_stimtime(fn, EDS_stimtime)	

	fn = '/data/backed_up/shared/ThalHi_MRI/ScanLogs/%s_%s_Switch_stimtime.1D' %(s, postfix)
	write_stimtime(fn, Switch_stimtime)



if __name__ == "__main__":

	#Subject, ROI, nruns = raw_input().split()

	ntrials_per_run = 50
	num_runs=2
	#nruns = 12
	s='D001'
	fns = ['D01_002_Task_THHS_2019_Apr_09_1051.csv', 'D01_004_Task_THHS_2019_Apr_09_1110.csv']
	postfix='MB2'
	parse_stim(s, fns, ntrials_per_run, num_runs, postfix)

	fns = ['D01_001_Task_THHS_2019_Apr_09_1042.csv', 'D01_003_Task_THHS_2019_Apr_09_1101.csv']
	postfix='MB4'
	parse_stim(s, fns, ntrials_per_run, num_runs, postfix)

	
