function Y = Generate_tau(lambda)
    % lambda=[lambda(p1),lambda(p2)];
    Y(1) = lambda(1)*1+lambda(2)*0 + 4*mean(lambda)*0.5;
    Y(2) = lambda(1)*0+lambda(2)*1 + 4*mean(lambda)*0.5;
end