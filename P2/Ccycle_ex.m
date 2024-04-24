clear;
% clc;
% 从最粗网格层迭代
load('A-B.mat');    load('Uref.mat');
% 参数
E =1;  % Young modulus
k = 0.3;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operato r i ~= j;
tol = 1e-6; g = 0.06;   rho=2;   i1=2;  i2=1;
uerror = zeros(8,3); href=1/(length(U1ref)-1);
[noderef,elemref,Db,Nb_L,Nb_R,Cb]=generation_elem_node(href);
u0=[]; Lambda=[]; restart=[]; M=[];
maxN = 8;
minN = 3;
for l=minN:maxN
    h=1/(2^(l+1));
    A=allA{l};  B=allB{l}; % 其中B={B{1},B{2}};
    if l > minN
        u0 = I2htoh{l}*u_;
        [m,n]=size (I2htoh{l});
        I=I2htoh{l}(1:m/2,1:n/2);
        Lambda=I*Lambda;
    end
    % u0=[]; Lambda=[];   %注意
    tic;
    [u_,Lambda,resvec,res,iter,Niter] = Friction_Solve(K1,K2,allA,B,h,u0,l,Lambda,restart,M,tol,g,rho,I2htoh,Rhto2h,i1,i2);
    t = toc;
    fprintf("网格1/ %d的求解时间:%.4f(s)\n",1/h,t);
    % 
    u=u_;
%     if l>minN
%          [node,elem,Db,Nb_L,Nb_R,Cb]=generation_elem_node(h);
%          Cbpoint=unique(Cb);
%          ux = u_(Cbpoint);
%          subplot(1,2,2);
%          x = (0:1/(length(Cbpoint)-1):1)';
%          plot(x,ux); hold on;
%     end
    while h>href
        [I,R] = IR_FW(h);
        u = I*u;
        h=h/2;
    end
    % 返回二维数据
    N = sqrt(length(u)/2); 
    u1=u(1:length(u)/2);u2=u(length(u)/2+1:end);
    U1=reshape(full(u1),N,N);U2=reshape(full(u2),N,N);
    err = Q_error(K1,K2,U1,U2,U1ref,U2ref,elemref,noderef);
    sss=[iter,err,Niter];
    uerror(l,:) = sss;
end


%% 计算Q范数误差
for i = 3:1:maxN-1
    orderu = uerror(i-1,2)/uerror(i,2);
    fprintf('1/%d层与1/%d层间收敛阶:%.4f\n',2^(i),2^(i+1),log2(orderu));
end   