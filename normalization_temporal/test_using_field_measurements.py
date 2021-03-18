from util.myfun import *
import time
from datetime import datetime
from normalization_temporal.fitting_DTC import *


infile = r'D:\data\tvs.xlsx'

tstr = read_excel_sheet_col(infile,0,'Sheet1')[1:]
Ts = read_excel_sheet_col(infile,1,'Sheet1')[1:]
Tv = read_excel_sheet_col(infile, 2, 'Sheet1')[1:]



date = [xlrd.xldate_as_tuple(t,0) for t in tstr]
year = [perdate[0] for perdate in date]
doy = [date2DOY(perdate[0],perdate[1],perdate[2]) for perdate in date]
localt = [perdate[3]+perdate[4]/60.0 for perdate in date]

doy = np.asarray(doy)
year = np.asarray(year)
localt = np.asarray(localt)
Tv = np.asarray(Tv,dtype = np.float)
Ts = np.asarray(Ts,dtype = np.float)


ind = (year == 2014)*(doy == 310)*(localt > 9)*(localt < 20)
tt = localt[ind]
Tvv = Tv[ind]
Tss = Ts[ind]
dTvs = Tss - Tvv

x=fitting_DTC(tt,dTvs)
mod = model_DTC(x,tt)

plt.plot(tt,dTvs,'.')
plt.plot(tt,mod,'-')
plt.show()



