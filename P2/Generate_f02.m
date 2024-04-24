function int_fphi = Generate_f2(K1,K2,x,y)
%         z = (K1+K2)*(6*y.*cos(pi*x))...
%         + (K1+0.5*K2)*(-y.*exp(x).*(3*y - 2))...
%         + 0.5*K2*(-pi^2*y.*cos(pi*x).*(y.^2 - 1));      %f0(2)²¿·Ö
        f12 = 0*x;
        int_fphi=zeros(size(x,1),1);
        int_fphi(1)=[0.5 0 0.5]*f12;
        int_fphi(2)=[0.5 0.5 0]*f12;
        int_fphi(3)=[0 0.5 0.5]*f12;
end 