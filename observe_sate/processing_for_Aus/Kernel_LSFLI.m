function sse=Kernel_LSFLI(params,Input,Actual_Output)

A=params(1);
B=params(2);
C=params(3);


in1=Input(:,1);
in2=Input(:,2);

Fitted_Curve=A.*in1 + B.*in2 +C;
Error_Vector=Fitted_Curve - Actual_Output;
Var_Data=Actual_Output-mean(Actual_Output);
sse=(sum(Error_Vector.^2))./(sum(Var_Data.^2));

end