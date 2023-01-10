subjects = [10001,10002,10003,10004,10005,
10008,10010,10012,10013,10014,
10017,10018,10019,10020,10022,
10023,10024,10025,10027,10028,
10031,10032,10033,10034,10035,
10036,10037,10038,10039,10040,
10041,10042,10043,10044,10054,
10057,10058,10059,10060,10063,
10064,10066,10068,10069,10071,
10072,10073,10074,10076,10077,
10080,10162,10169,10170,10173,
10174,10175,10176,10179]

All_reliability = np.zeros(len(subjects))
for i, sub in enumerate(subjects):
    a = np.loadtxt("/home/kahwang/argon/data/ThalHi/reliability_test/sub-%s/All_1-4_beta.1D" %sub)
    b = np.loadtxt("/home/kahwang/argon/data/ThalHi/reliability_test/sub-%s/All_5-8_beta.1D" %sub)
    All_reliability[i] = np.corrcoef(a,b)[0,1]