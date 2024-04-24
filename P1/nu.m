function y = nu(x)
    u_nu1 = 0.04;u_nu2 = 0.06;
    c_nu1 = 1;c_nu2 = -0.5;c_nu3 = 4;
    y = zeros(length(x),1); 
    for i = 1:1:length(x)
        if x(i) < 0
            y(i) = 0;
        elseif x(i) < u_nu1
            y(i) = c_nu1*x(i);
%             fprintf('--');
        elseif x(i) < u_nu2
            y(i) = c_nu1*u_nu1+c_nu2*(x(i)-u_nu1);
%             fprintf('?');
        else
            y(i) = c_nu1*u_nu1+c_nu2*(u_nu2-u_nu1)+c_nu3*(x(i)-u_nu2);
%             fprintf('!');
        end
    end
    y=y/1.2;
end