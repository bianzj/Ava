import numpy as np
import time
from datetime import datetime
from normalization_temporal.fitting_DTC import *
import pandas as pd

infile = r'D:\data\AWS_HHL_2015.csv'
infile_ = r'D:\data\TVS_HHL_2015.csv'

month = 9
day = 6

aws = pd.read_csv(infile, parse_dates=['TIMESTAMP'])
low = datetime.datetime(2015, month, day, 0)
high = datetime.datetime(2015, month, day+1, 0)
aws = aws[aws['TIMESTAMP'] >= low]
aws = aws[aws['TIMESTAMP'] < high]
Time_aws = aws['TIMESTAMP']
doy = Time_aws.apply(lambda x:x.dayofyear)
hour = Time_aws.apply(lambda x:x.hour)
minute = Time_aws.apply(lambda x:x.minute)
t_aws = np.asarray(doy+ hour/24.0 + minute/60.0/24.0,dtype=np.float)
Ta_aws = np.asarray(aws['Ta_5_Avg'].values,dtype=np.float)
Ws_aws = np.asarray(aws['WS_5_Avg'].values,dtype=np.float)
Dl_aws = np.asarray(aws['DLR_CO_Avg'].values,dtype=np.float)
Ds_aws = np.asarray(aws['DR_Avg'].values,dtype=np.float)
Ul_aws = np.asarray(aws['ULR_CO_Avg'].values,dtype=np.float)
Us_aws = np.asarray(aws['UR_Avg'].values,dtype=np.float)
IRT_aws = np.asarray(aws['IRT_1_Avg'].values,dtype=np.float)

ct = pd.read_csv(infile_, parse_dates=['Time'])
ct = ct[ct['Time'] >= low]
ct = ct[ct['Time'] < high]
Time_ct = ct['Time']
doy = Time_ct.apply(lambda x:x.dayofyear)
hour = Time_ct.apply(lambda x:x.hour)
minute = Time_ct.apply(lambda x:x.minute)
t_ct = np.asarray(doy+ hour/24.0 + minute/60.0/24.0,dtype=np.float)
Tv_ct = np.asarray(ct['Tv'].values,dtype=np.float)-273.15
Ts_ct = np.asarray(ct['Ts'].values,dtype=np.float)-273.15

plt.plot(t_ct,Ts_ct - Tv_ct)
# plt.plot(t_aws,(-Us_aws+Ds_aws)/40,'r')
# plt.plot(t_aws,(-Dl_aws+Ul_aws)/10,'r')
plt.plot(t_aws,IRT_aws-Ta_aws,'r')
plt.show()