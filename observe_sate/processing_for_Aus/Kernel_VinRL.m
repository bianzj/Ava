function sse=Kernel_VinRL(params,Input,Actual_Output)

C=params(1);
A=params(2);
B=params(3);
D=params(4);

in1=Input(:,1);
in2=Input(:,2);
in3=Input(:,3);
Fitted_Curve=C.*in1 +  A.*(exp(-B.*in2)-exp(-B.*tan(in3)))./(1-exp(-B.*tan(in3)))+D;
Error_Vector=Fitted_Curve - Actual_Output;
Var_Data=Actual_Output-mean(Actual_Output);
sse=(sum(Error_Vector.^2))./(sum(Var_Data.^2));

end