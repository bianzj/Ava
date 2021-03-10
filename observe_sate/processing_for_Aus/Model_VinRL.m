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