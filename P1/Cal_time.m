clear;
% clc;
% �������������
load('A-B.mat');    load('Uref.mat');

% ����
E = 2;  % Young modulus
k = 0.4;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operato r i ~= j;
tol = 1e-6; g = 0.07;   rho=2.5;   i1=2;  i2=2;
uerror = zeros(8,3); href=1/(length(U1ref)-1);
[noderef,elemref,Db,Nb_L,Nb_R,Cb]=generation_elem_node(href);
maxN = 8;
minN = 3;
for kk=5:8
    u0=[]; Lambda=[]; restart=[]; M=[];
    tic;
    for l=minN:kk
        h=1/(2^(l+1));
        A=allA{l};  B=allB{l}; % ����B={B{1},B{2}};
        if l > minN
            u0 = I2htoh{l}*u_;
            [m,n]=size (I2htoh{l});
            I=I2htoh{l}(1:m/2,1:n/2);
            Lambda=I*Lambda;
        end
        u0=[]; Lambda=[];   %ע�����޳�ʼֵ
        [u_,Lambda,resvec,res,iter,Niter,B] = Friction_Solve(K1,K2,allA,B,h,u0,l,Lambda,restart,M,tol,g,rho,I2htoh,Rhto2h,i1,i2);
    end
    t = toc;
    fprintf("����1/ %d�����ʱ��:%.4f(s),ƽ���ڵ�������:%d\n",1/h,t,Niter);
end
