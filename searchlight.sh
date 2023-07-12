# run searchlight on thalahi dataset.
#$ -N thalhi_sl
#$ -q SEASHORE
#$ -pe smp 4
#$ -t 1-50
#$ -tc 21
#$ -ckpt user
#$ -o /Users/kahwang/sge_logs/
#$ -e /Users/kahwang/sge_logs/
/bin/echo Running on compute node: `hostname`.
source activate py3.8
/bin/echo Job: $JOB_ID
/bin/echo Task: $SGE_TASK_ID
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

subjects=(10001  10014  10027  10038  10068 10002  10017  10028  10040  10071 10003  10018  10029  10041  10073 10004  10019  10030  10042  10074 10005  10020  10031  10043  10076 10007  10021  10032  10054  10077 10008  10022  10033  10059  10080 10010  10023  10034  10060 10024  10035  10063 10012  10025  10036  10064 10013  10026  10037  10066)
echo subjects: ${subjects[@]}
echo total_subjects=${#subjects[@]}
subject="${subjects[$SGE_TASK_ID-1]}"
echo "Starting decoding on $subject"
echo $subject | python /Users/kahwang/bin/ThalHi_MRI/test_searchlight.py