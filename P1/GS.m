function [x0,j,res] = GS(A,b,i,x0)
    if isempty(x0) 
        x0 = zeros(size(b));
    end
    res = [];
    D = diag(diag(A));
    L = D-tril(A);
    U = D-triu(A);
    flag = 1;
    for j =1:1:i
        x0 = (D-L)\(b+U*x0);
        res = [res;norm(b-A*x0)/norm(b)];
        if res(end) <1e-6
            flag= 0;
            break;
        end
    end
    if flag ==1 && i>10
        fprintf('%d 步内没有迭代收敛',i);
    end
end