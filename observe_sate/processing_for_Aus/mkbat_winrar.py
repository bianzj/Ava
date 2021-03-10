from processing_for_Aus.myfun import *

indir = 'J:\SLSTR3B\\'
outdir = indir

outfile = outdir+'winrar5.bat'
file_start = 280
file_end = 358
outfilenum = 0

infiles = ope_search1(indir,'S3B')
f = open(outfile,'w')

filenum = 0
for k in range(len(infiles)):
    infile = infiles[k]
    passdate = infile[16:8 + 16]
    passtime = infile[25:25 + 6]

    if (passtime > '090000') and (passtime < '200000'):
        continue
    print(filenum,passdate,passtime)
    filenum = filenum+1
    if (filenum > file_start) and (filenum<= file_end):
        f.write('winrar x '+infile+'\n')


f.close()




