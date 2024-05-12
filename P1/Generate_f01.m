function int_fphi = Generate_f01(K1,K2,x,y)
        f11 = 0*x;
        int_fphi=zeros(size(x,1),1);
        int_fphi(1)=[0.5 0 0.5]*f11;
        int_fphi(2)=[0.5 0.5 0]*f11;
        int_fphi(3)=[0 0.5 0.5]*f11;
end
