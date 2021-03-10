function [Estimates,R2]=fitting_VinRL(myfit,x,y,n)

starting=[1.1 1.1 1.1 15];
options=optimset('Display','off','TolX',1e-6,'TolFun',1e-6);
Estimates=fminsearch(myfit,starting,options,x,y);
R2=1-myfit(Estimates,x,y);

end