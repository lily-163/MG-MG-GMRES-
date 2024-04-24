function int_fphi = Generate_f1(K1,K2,x,y)
%         z = (K1+K2)*(-y.^2.*exp(x).*(y - 1))...
%         + (K1+0.5*K2)*(-pi*sin(pi*x).*(3*y.^2 - 1))...
%         + 0.5*K2*(-2*exp(x).*(3*y - 1));     % f0(1)����
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