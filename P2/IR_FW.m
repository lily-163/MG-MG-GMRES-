function [I,R]=IR_FW(H)
    % H:为粗网格层的划分H，插值到h=H/2的网格上
    n=1/H+1;
    N=2*n-1;
    M1=sparse(1:2:N,1:n,1,N,n)+sparse(2:2:N-1,1:n-1,0.5,N,n)+sparse(2:2:N-1,2:n,0.5,N,n);
    M2=M1/2;
    T1=sparse(1:2:N,1:n,1,N,n);
    T2=sparse(2:2:N-1,1:n-1,1,N,n);
    T3=sparse(2:2:N-1,2:n,1,N,n);
    I = kron(T1,M1)+kron(T2,M2)+kron(T3,M2);
    R = I';
    R=R./sum(R,2);
    egg=sparse(N^2,n^2);
    I = [I,egg;egg,I];
    R = [R,egg';egg',R];
end
