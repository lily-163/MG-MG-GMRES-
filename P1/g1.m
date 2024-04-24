function z = g1(K1,K2,x,y)
%         z = (K1+K2)*(-y.^2.*exp(x).*(y - 1))...
%             +K1*(cos(pi*x).*(3*y.^2 - 1));
    z = -0.3*sin(6*pi*y-pi/2)-0.3-2*y*(1-y);
%     if y>=0.2 && y<1
%         z = +5*(y-0.2)*(y-1);
%     else
%         z = 0;
%     end
end