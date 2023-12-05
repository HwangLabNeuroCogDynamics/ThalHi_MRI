subjects=(10001 10002 10003 10004 10005 \
10008 10010 10012 10013 10014 \
10017 10018 10019 10020 10022 \
10023 10024 10025 10027 10028 \
10031 10032 10033 10034 10035 \
10036 10037 10038 10039 10040 \
10041 10042 10043 10044 10054 \
10057 10058 10059 10060 10063 \
10064 10066 10068 10069 10071 \
10072 10073 10074 10076 10077 \
10080 10162 10169 10170 10173 \
10174 10175 10176 10179)

for sub in ${subjects[@]}
do
3dFWHMx -acf -input /mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/3dDeconvolve_fdpt4/sub-${sub}/sub-${sub}_FIRmodel_RT_model_errts_REML+tlrc | tail -1 >> /mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/RT_investigations/acf_parameters_59subs_RT.txt
done

3dClustSim -acf 0.928 1.605 12.492 -mask /mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/3dMEMA/EDS-IDS_59subs.nii.gz >> /mnt/nfs/lss/lss_kahwang_hpc/data/ThalHi/RT_investigations/a9subs_cluststim.txt