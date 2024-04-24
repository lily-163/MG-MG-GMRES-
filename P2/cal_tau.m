function z = cal_tau(lambda,u)
    a=0.05;b=0.03;beta=100;
    up=(a-b)*exp(-beta*abs(u))+b;
    if norm(lambda)<= norm(up)
        z=lambda;
%         disp('!!!!!!')
    else
        z=lambda*up/(abs(lambda));
%         disp('спак')
    end
end