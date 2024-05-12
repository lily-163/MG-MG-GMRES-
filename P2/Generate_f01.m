function int_fphi = Generate_f1(K1,K2,x,y)
        f11 = 0*x;
        for i = 1:length(f11)
            if x(i)>0.5 
                f11(i) = 0.05*(exp(5*abs(x(i)-0.5))-1);
            elseif x(i)<0.5
                f11(i) = -0.05*(exp(5*abs(x(i)-0.5))-1);
            end
        end
        int_fphi=zeros(size(x,1),1);
        int_fphi(1)=[0.5 0 0.5]*f11;
        int_fphi(2)=[0.5 0.5 0]*f11;
        int_fphi(3)=[0 0.5 0.5]*f11;
end
