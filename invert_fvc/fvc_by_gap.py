from util.myfun import *



indir_lai = r'D:\data\lst\lai_AREA\\'
indir_ci = r'D:\data\lst\CI_AREA\\'
outdir_fvc = r'd:\data\lst\fvc_area\\'

G = 0.5
days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
kdoy = 0
for kmonth in range(12):
    daynum = days_of_month[kmonth]

    infile0 = indir_ci + '2019_month%02d' % (kmonth+1) + '.tif'
    [ci, ns, nl, nb, goeg, proj] = read_image_gdal(infile0)

    for kday in range(daynum):

        kdoyby5 = (kdoy//4)

        infile1 = indir_lai + '2019_doy%03d'%(kdoyby5*4+1)+'.tif'
        infile2 = indir_lai + '2019_doy%03d'%(kdoyby5*4+5)+'.tif'
        if kdoyby5*4+5 >= 365: infile2 = indir_lai + '2019_doy365.tif'
        prop = (kdoy % 4) / 4.0
        print(infile1,infile2, prop)
        [lai1,ns,nl,nb,geog,proj] = read_image_gdal(infile1)
        [lai2,ns,nl,nb,geog,proj] = read_image_gdal(infile2)

        lai = lai1 * (1-prop) + lai2 *prop

        fvc = 1-np.exp(-G*lai*ci)
        fvc[fvc>0.995] = 0

        kdoy = kdoy + 1

        outfile = outdir_fvc + '2019_doy%03d'%kdoy+'_fvc.tif'
        write_image_gdal(fvc,ns,nl,nb,geog,proj,outfile)