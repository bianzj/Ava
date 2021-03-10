import shutil
from processing_for_Aus.myfun import *

wdir = r'J:\SateTDR\\'
fileNames = ope_searchfile_accept1_rej1(wdir,'S3A','ZIP')

bt = 'emis'
for k in range(len(fileNames)):
    fileName = fileNames[k]
    UTCtime = fileName[25:27]
    if (UTCtime > '05') & (UTCtime < '20'):
        file = wdir + fileName
        print(file)
        try:
            # os.rmdir(path)
            # os.remove(file)
            shutil.rmtree(file)
        except Exception as e:
            print('Exception', e)
