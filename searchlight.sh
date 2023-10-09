# run searchlight on thalahi dataset.
#$ -N thalhi_sl
#$ -q SEASHORE
#$ -pe smp 24
#$ -t 1-73
#$ -tc 6
#$ -ckpt user
#$ -o /Users/kahwang/sge_logs/
#$ -e /Users/kahwang/sge_logs/
/bin/echo Running on compute node: `hostname`.
source activate py3.8
/bin/echo Job: $JOB_ID
/bin/echo Task: $SGE_TASK_ID
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

subjects=(10001  10006  10011  10017  10022  10027  10032  10037  10042  10057  10062  10068  10074  10169  10175 10002  10007  10012  10018  10023  10028  10033  10038  10043  10058  10063  10069  10076  10170  10176 10003  10008  10013  10019  10024  10029  10034  10039  10044  10059  10064  10071  10077  10172  10179 10004  10009  10014  10020  10025  10030  10035  10040  10054  10060  10065  10072  10080  10173 10005  10010  10016  10021  10026  10031  10036  10041  10055  10061  10066  10073  10162  10174)
echo subjects: ${subjects[@]}
echo total_subjects=${#subjects[@]}
subject="${subjects[$SGE_TASK_ID-1]}"
echo "Starting decoding on $subject"
echo $subject | python /Users/kahwang/bin/ThalHi_MRI/test_searchlight.py

