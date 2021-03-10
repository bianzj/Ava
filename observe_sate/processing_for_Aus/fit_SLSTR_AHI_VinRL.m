clc;
clear;

indir1 = r'J:\SLSTR_day_lst\\';
indir2 = r'J:\AHI8_day_selected\\';
indir3 = r'J:\SLSTR_day\\';
outdir = r'J:\\figure\\';

lstmin = 250;
lstmax = 350;
ndvimin = 0;
ndvimax = 1.0;

rmse_ = [];
for k = 1:31
   
    k
    datestr(now)
    
    fileName = ['Aus_201801' num2str(k,'%02d') '_'];
    
    infile1 = [indir1 fileName 'LST_n.tif'];
    infile2 = [indir1 fileName 'LST_o.tif'];
    infile3 = [indir2 fileName 'LST.tif'];
    
    infile4 = [indir3 fileName 'psi_n.tif'];
    infile5 = [indir3 fileName 'psi_o.tif'];
    infile6 = [indir2 fileName 'psi.tif'];
    
    infile7 = [indir3 fileName 'vza_n.tif'];
    infile8 = [indir3 fileName 'vza_o.tif'];
    infile9 = [indir2 fileName 'vza.tif'];
    
    infile10 = [indir3 fileName 'cloud_n.tif'];
    infile11 = [indir3 fileName 'cloud_o.tif'];
    
    infile12 = [indir3 fileName 'sza_n.tif'];
    
    mydata = Tiff(infile1,'r+');
    mydata.setDirectory(1);
    lst_slstr_n = mydata.read();
    lst_slstr_n = imresize(lst_slstr_n,0.5);
    
    mydata = Tiff(infile2,'r+');
    mydata.setDirectory(1);
    lst_slstr_o = mydata.read();
    lst_slstr_o = imresize(lst_slstr_o,0.5);
        
    mydata = Tiff(infile3,'r+');
    mydata.setDirectory(1);
    lst_ahi = mydata.read();
//     lst_ahi = imresize(lst_ahi);
        
    mydata = Tiff(infile4,'r+');
    mydata.setDirectory(1);
    psi_n = mydata.read();
    psi_n = imresize(psi_n,0.5);
        
    mydata = Tiff(infile5,'r+');
    mydata.setDirectory(1);
    psi_o = mydata.read();
    psi_o = imresize(psi_o,0.5);
    
    mydata = Tiff(infile6,'r+');
    mydata.setDirectory(1);
    psi_ahi = mydata.read();
%     psi_ahi = imresize(psi_ahi,2);
    
    mydata = Tiff(infile7,'r+');
    mydata.setDirectory(1);
    vza_n = mydata.read();
    vza_n = imresize(vza_n,0.5);
    
    mydata = Tiff(infile8,'r+');
    mydata.setDirectory(1);
    vza_o = mydata.read();
    vza_o = imresize(vza_o,0.5);
     
        mydata = Tiff(infile9,'r+');
    mydata.setDirectory(1);
    vza_ahi = mydata.read();
%     vza_ahi = imresize(vza_ahi,2);

    mydata = Tiff(infile10,'r+');
    mydata.setDirectory(1);
    cloud_n = mydata.read();
    cloud_n = imresize(cloud_n,0.5);
    
    mydata = Tiff(infile11,'r+');
    mydata.setDirectory(1);
    cloud_o = mydata.read();
    cloud_o = imresize(cloud_o,0.5);
    
    mydata = Tiff(infile12,'r+');
    mydata.setDirectory(1);
    sza_n = mydata.read();
    sza_n = imresize(sza_n,0.5);

    diflst = 15;
    ind0 = (lst_slstr_n<lstmax).*(lst_slstr_n>lstmin).* ...
          (lst_slstr_o<lstmax).*(lst_slstr_o>lstmin).* ...
          (lst_ahi<lstmax).*(lst_ahi>lstmin).*...
          (cloud_n==0).*(cloud_o==0).*(abs(lst_slstr_n-lst_slstr_o)<diflst).*...
          (abs(lst_ahi-lst_slstr_n)<diflst).*(abs(lst_ahi-lst_slstr_o)<diflst);
    
    
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
    
    coeffa = zeros(nl,ns);
    coeffb = zeros(nl,ns);
    coeffc = zeros(nl,ns);
    coeffd = zeros(nl,ns);
    
    for ks =1:ns
        for kl = 1:nl
        
%         offl = kl*2-1;
%         endl = kl*2;
%         offs = ks*2-1;
%         ends = ks*2;
        
        offl = kl;
        endl = kl;
        offs = ks;
        ends = ks;
        
        ind = ind0(offl:endl,offs:ends);
        ind = ind == 1;
        
        num = sum(ind(:));
        if(num < 1) continue
        end
        lst3 = lst_ahi(offl,offs);
        psi3 = psi_ahi(offl,offs);
        ind3 = (lst3 < lstmax)*(lst3 > lstmin);
        if(ind3 < 1) continue
        end
        
        
        lst1 = lst_slstr_n(offl:endl,offs:ends);
        lst2 = lst_slstr_o(offl:endl,offs:ends);        
        psi1 = psi_n(offl:endl,offs:ends);
        psi2 = psi_o(offl:endl,offs:ends);       
        vza1 = vza_n(offl:endl,offs:ends);
        vza2 = vza_o(offl:endl,offs:ends);
        vza3 = vza_ahi(offl,offs);
        sza1 = sza_n(offl:endl,offs:ends);
        sza2 = sza1;
        sza3 = sza1(1);
        
                
        vza = [vza1(ind); vza2(ind); vza3];
        lst = [lst1(ind); lst2(ind); lst3];
        psi = [psi1(ind); psi2(ind); psi3];
        sza = [sza1(ind);sza2(ind); sza3];
        
        f = sqrt((tand(sza).^2+tand(vza).^2-2.*tand(vza).*tand(sza).*cosd(psi)));
        kemis = 1-cosd(vza);
        X = [f kemis sza];
        Y = lst;
        [Estimates,R2]=fitting_VinRL(@Kernel_VinRL,X,Y,4);
        
        coeffa(kl,ks) = Estimates(1);
        coeffb(kl,ks) = Estimates(2);
        coeffc(kl,ks) = Estimates(3);
        coeffd(kl,ks) = Estimates(4);
        
%         predicted =Model_VinRL(Estimates,X);
%         dif = predicted - Y;
%         if (sum(isnan(dif))>0) continue
%         end
        
        
        
%         yy_ = [yy_;predicted];
%         y_ = [y_;Y];
        
                     
        end
    end
    
    
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_VinRLa.tif'];
    saveTif(coeffa,outfile);
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_VinRLb.tif'];
    saveTif(coeffb,outfile);
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_VinRLc.tif'];
    saveTif(coeffc,outfile);
    outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_VinRLd.tif'];
    saveTif(coeffd,outfile);
    
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
%     outfile = [outdir 'Aus_201801' num2str(k,'%02d') '_VinRL.tif'];
%     f=getframe(gcf);
%     imwrite(f.cdata,outfile);
        
end

% datestr(now)
% save('rmse_VinRL.mat','rmse');

function newdbt = Model_VinRL(Estimates,Input)

    in1=Input(:,1);
    in2=Input(:,2);
    in3=Input(:,3);
    C=Estimates(1);
    A=Estimates(2);
    B=Estimates(3);
    D=Estimates(4);
    newdbt = C.*in1 +  A.*(exp(-B.*in2)-exp(-B.*tan(in3)))./(1-exp(-B.*tan(in3)))+D;

end
