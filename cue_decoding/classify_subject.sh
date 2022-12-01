#!/bin/bash
# SGE
#$ -N LDA_SEARCHLIGHT
#$ -q SEASHORE
#$ -pe smp 80
#$ -t 1-1
# OMP_NUM_THREADS=80

# subjects=(10001 10002 10003 10004 10005 10007 10008 10010 10011 10012 10013 10014 10017 10018 10019 10020 10021 10022 10023 10024 10025 10026 10027 10028 10029 10030 10031 10032 10033 10034 10035 10036 10037 10038 10040 10041 10042 10043 10054)
subjects=(10001)

classifier=LDA
echo subjects: ${subjects[@]}
echo total_subjects=${#subjects[@]}
subject="${subjects[$SGE_TASK_ID-1]}"

source /Shared/lss_kahwang_hpc/virtualenvs/base/bin/activate

python3 /Shared/lss_kahwang_hpc/scripts/thalhi/cue_decoding/decoding_subject.py $subject --searchlight --cores 80 --classifier LDA