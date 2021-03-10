from myfun_file import *
from myfun_image import *


###################################################
####  RENAME OR MOVE
###################################################
# indir = r'G:\s3b_tif\\'
# outdir = r'G:\s3b_tif\\'
# dirNames = search_dir_rej(indir,['SEN3'],'zip')
# num = len(dirNames)
# for k in range(num):
#     source_file = indir+dirNames[k]+'/refl_red.tif'
#     dest_file = outdir+dirNames[k]+'/red.tif'
#     move_file(source_file,dest_file)


# indir = r'G:\mod_tif\\'
# outdir = r'G:\mod_tif\\'
# dirNames = search_file(indir,['bt31'])
# num = len(dirNames)
# for k in range(num):
#     source_file = indir+dirNames[k]
#     dest_file = outdir+dirNames[k][:-8]+'bt1.tif'
#     print(source_file,dest_file)
    # move_file(source_file,dest_file)

###################################################
####  REMOVE FILES
###################################################
# indir = r'G:\s3b_raw\\'
# target_file = 'cloud.tif'
# dirNames = search_dir_rej(indir,['SEN3'],'zip')
# num = len(dirNames)
# for k in range(1,num-1):
#     target_dir = indir+dirNames[k]
#     remove_file(target_dir,target_file)


#################################################
##### WRITE WINRAR.BAT FOR UNZIP FILE
##################################################
indir = r'F:/SLSTR3A/04/'
outfile = 'F:/SLSTR3A/04/winrar.bat'
f = open(outfile,'w')
fileNames = search_file_rej(indir,'zip','SEN3')
fileNum = len(fileNames)
for k in range(fileNum):
    f.write('winrar x '+fileNames[k] +'\n')
f.close()


# indir = r'G:/s3a_raw/'
# outfile = 'G:/s3a_raw/winrar.bat'
# f = open(outfile,'w')
# fileNames = search_file_rej(indir,'zip','SEN3')
# fileNum = len(fileNames)
# for k in range(fileNum):
#     if os.path.exists(indir+fileNames[k][:-3]+'SEN3')==1:
#         continue
#     f.write('winrar x '+fileNames[k] +'\n')
# f.close()