function int_fphi = Generate_f02(K1,K2,x,y)
        f12 = -0.5*10^(-1)*ones(size(x));
        int_fphi=zeros(size(x,1),1);
        int_fphi(1)=[0.5 0 0.5]*f12;
        int_fphi(2)=[0.5 0.5 0]*f12;
        int_fphi(3)=[0 0.5 0.5]*f12;
end 
