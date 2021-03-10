clc;
clear;

indir1 = 'J:\SLSTR3A_day_lst\\';
indir2 = 'J:\AHI8_day_selected\\';
indir3 = 'J:\SLSTR3a_day\\';
indir4 = 'J:\AHI8_day\\'
outdir = 'J:\multiple_result\\';

lstmin = 250;
lstmax = 350;
ndvimin = 0;
ndvimax = 1.0;

rmse_ = [];
for k = 2:2
   
    k
    datestr(now)
    
    fileName = ['Aus_201812' num2str(k,'%02d') '_'];
    
    infile1 = [indir1 fileName 'LST_n.tif'];
    infile2 = [indir1 fileName 'LST_o.tif'];
    infile3 = [indir2 fileName 'LST.tif'];
    
    infile4 = [indir3 fileName 'refl_red_n.tif'];
    infile5 = [indir3 fileName 'refl_red_o.tif'];
    infile6 = [indir2 fileName 'red.tif'];
    
    infile7 = [indir3 fileName 'refl_nir_n.tif'];
    infile8 = [indir3 fileName 'refl_nir_o.tif'];
    infile9 = [indir2 fileName 'nir.tif'];
    
    infile10 = [indir3 fileName 'cloud_n.tif'];
    infile11 = [indir3 fileName 'cloud_o.tif'];
    
    %% lst
     mydata = Tiff(infile1,'r+');
    mydata.setDirectory(1);
    lst_slstr_n = mydata.read();
%     lst_slstr_n = imresize(lst_slstr_n);
    
    mydata = Tiff(infile2,'r+');
    mydata.setDirectory(1);
    lst_slstr_o = mydata.read();
%     lst_slstr_o = imresize(lst_slstr_o);
        
    mydata = Tiff(infile3,'r+');
    mydata.setDirectory(1);
    lst_ahi = mydata.read();
    lst_ahi = imresize(lst_ahi,2);
        
    %% red and nir
    mydata = Tiff(infile4,'r+');
    mydata.setDirectory(1);
    red_slstr_n = mydata.read();
%     red_slstr_n = imresize(red_slstr_n,0.5);
        
    mydata = Tiff(infile5,'r+');
    mydata.setDirectory(1);
    red_slstr_o = mydata.read();
%     red_slstr_o = imresize(red_slstr_o,0.5);
    
        mydata = Tiff(infile6,'r+');
    mydata.setDirectory(1);
    red_ahi = mydata.read();
    red_ahi = imresize(red_ahi,2);
    
        mydata = Tiff(infile7,'r+');
    mydata.setDirectory(1);
    nir_slstr_n = mydata.read();
%     nir_slstr_n = imresize(nir_slstr_n,0.5);
    
        mydata = Tiff(infile8,'r+');
    mydata.setDirectory(1);
    nir_slstr_o = mydata.read();
%      nir_slstr_o = imresize(nir_slstr_o,0.5);
     
        mydata = Tiff(infile9,'r+');
    mydata.setDirectory(1);
    nir_ahi = mydata.read();
    nir_ahi = imresize(nir_ahi,2);

        mydata = Tiff(infile10,'r+');
    mydata.setDirectory(1);
    cloud_n = mydata.read();
%     cloud_n = imresize(cloud_n,0.5);
    
    mydata = Tiff(infile11,'r+');
    mydata.setDirectory(1);
    cloud_o = mydata.read();
%     cloud_o = imresize(cloud_o,0.5);
    
    ndvi_slstr_n = (nir_slstr_n - red_slstr_n)./(nir_slstr_n+red_slstr_n);
     ndvi_slstr_o = (nir_slstr_o - red_slstr_o)./(nir_slstr_o+red_slstr_o);
      ndvi_ahi = (nir_ahi - red_ahi)./(nir_ahi+red_ahi);
    
    bf_slstr_n = (nir_slstr_n + red_slstr_n)/2;
     bf_slstr_o = (nir_slstr_o + red_slstr_o)/2;
      bf_ahi = (nir_ahi + red_ahi)/2;

    difndvi = 0.5;
    ind1 = (ndvi_slstr_n < ndvimax) .* (ndvi_slstr_n > ndvimin) .* ...
          (ndvi_slstr_o < ndvimax) .* (ndvi_slstr_o > ndvimin) .* ...
          (ndvi_ahi < ndvimax) .* (ndvi_ahi > ndvimin) .* ...
          (cloud_n==0) .* (cloud_o==0) .* (abs(ndvi_slstr_n - ndvi_slstr_o) < difndvi) .* ...
          (abs(ndvi_ahi - ndvi_slstr_n) < difndvi) .* (abs(ndvi_ahi - ndvi_slstr_o) < difndvi);

    diflst = 15;
    ind2 = (lst_slstr_n<lstmax).*(lst_slstr_n>lstmin).* ...
          (lst_slstr_o<lstmax).*(lst_slstr_o>lstmin).* ...
          (lst_ahi<lstmax).*(lst_ahi>lstmin).*...
          (cloud_n==0).*(cloud_o==0).*(abs(lst_slstr_n-lst_slstr_o)<diflst).*...
          (abs(lst_ahi-lst_slstr_n)<diflst).*(abs(lst_ahi-lst_slstr_o)<diflst);
    ind0 = ind1.*ind2;
    
    
%     ind = ind0==1;     
%     lst_slstr_n = lst_slstr_n(ind);
%     lst_slstr_o = lst_slstr_o(ind);
%     lst_ahi = lst_ahi(ind);
%     ndvi_slstr_n = ndvi_slstr_n(ind);
%     ndvi_slstr_o = ndvi_slstr_o(ind);
%     ndvi_ahi = ndvi_ahi(ind);
%     bf_slstr_n = bf_slstr_n(ind);
%     bf_slstr_o = bf_slstr_o(ind);
%     bf_ahi = bf_ahi(ind);
    
    
    [nl,ns] = size(ind0);
%     yy_ = zeros(num,3);
%     y_ = [lst_slstr_n, lst_slstr_o, lst_ahi];
    
    yy_ = [];
    y_ = [];
    nld2 = nl/2;
    nsd2 = ns/2;
    coeffa = zeros(nld2,nsd2);
    coeffb = zeros(nld2,nsd2);
    coeffc = zeros(nld2,nsd2);
    
    
    for ks =1:ns/2
        for kl = 1:nl/2
        
        offl = kl*2-1;
        endl = kl*2;
        offs = ks*2-1;
        ends = ks*2;
        
        ind = ind0(offl:endl,offs:ends);
        ind = ind == 1;
        
        num = sum(ind(:));
        if(num < 1) continue
        end
        lst3 = lst_ahi(offl,offs);
        ndvi3 = ndvi_ahi(offl,offs);
        ind3 = (lst3 < lstmax)*(lst3 > lstmin)*(ndvi3 < ndvimax)*(ndvi3 > ndvimin);
        if(ind3 < 1) continue
        end
        
        
        lst1 = lst_slstr_n(offl:endl,offs:ends);
        lst2 = lst_slstr_o(offl:endl,offs:ends);        
        ndvi1 = ndvi_slstr_n(offl:endl,offs:ends);
        ndvi2 = ndvi_slstr_o(offl:endl,offs:ends);       
        bf1 = bf_slstr_n(offl:endl,offs:ends);
        bf2 = bf_slstr_o(offl:endl,offs:ends);
        bf3 = bf_ahi(offl,offs);
        
        ndvi = [ndvi1(ind); ndvi2(ind); ndvi3];
        lst = [lst1(ind); lst2(ind); lst3];
        bf = [bf1(ind); bf2(ind); bf3];
        iso = ones(1,length(ndvi));
        
        x = [ndvi bf iso'];
        y = lst;
        
        if(length(y)<3)
            continue
        end
        
        [b,bint,r,rint,stats]=regress(y,x);
        
        coeffa(kl,ks) = b(1);
        coeffb(kl,ks) = b(2);
        coeffc(kl,ks) = b(3);
        
%         yy = x*b;
%         yy_ = [yy_;yy];
%         y_ = [y_;y];
        
        
                
        end
    end
    
    coeffa(coeffa<-300) = 0;
    coeffa(coeffa>360) = 0;
    coeffb(coeffb<-300) = 0;
    coeffb(coeffb>360) = 0;
     coeffc(coeffc<-300) = 0;
    coeffc(coeffc>360) = 0;
    
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_9VTa.tif'];
    saveTif(coeffa,outfile);
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_9VTb.tif'];
    saveTif(coeffb,outfile);
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_9VTc.tif'];
    saveTif(coeffc,outfile);
    
    
%     y_ = y_ -273.15;
%     yy_ = yy_ -273.15;
%     
%     data1 = y_;
%     data2 = yy_;
%     rmse = sqrt(mean((data1-data2).^2));
%     bias = mean(data2-data1);
%     coef = corrcoef(data1,data2);
%     r2 = coef(2)^2;
%     sigma = std(data2-data1);
%     
%     rmse_ = [rmse_ rmse];
%     
%     figure
%     plot(data1,data2,'k.');
%     set(gcf,'position',[400,400,450,400]);
%     [min1, min2, max1, max2, step1, step2] = deal(0, 0, 60, 60, 10, 10);
%     axis([min1 max1 min2 max2])
%     grid on;
%     set(gca,'XTick',[min1:step1:max1]) 
%     set(gca,'XTickLabel',[min1:step1:max1])
%     set(gca,'YTick',[min2:step2:max2]) 
%     set(gca,'YTickLabel',[min2:step2:max2])
%     xlabel('Observed LSTs (\circC)');
%     ylabel('Simulated LSTs (\circC)');
%     text(35,12,{['RMSE =' num2str(rmse,'%.2f') '°„C' ] ['Bias = ' num2str(bias,'%.2f') '°„C'] 
%         ['R^2 = ' num2str(r2,'%.2f')]  ['\sigma = ' num2str(sigma,'%.2f') '°„C']},...
%          'background','no','fontsize',15);
%     outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_VT.tif'];
%     f=getframe(gcf);
%     imwrite(f.cdata,outfile);
        
end
 
% datestr(now)
% save('rmse_VT.mat','rmse');
