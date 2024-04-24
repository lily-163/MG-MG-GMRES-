clear;
clc;
% �������������
load('A-B.mat');  
% ����
E = 2;  % Young modulus
k = 0.4;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operato r i ~= j;
tol = 1e-6; g = 0.07;   rho=2.5;   i1=2;  i2=2;
t = 8;
h = 1/(2^(t+1));
A=allA{t};  B=allB{t}; % ����B={B{1},B{2}};
[u,Lambda,resvec,res,iter,Niter,B] = Friction_Solve(K1,K2,allA,B,h,[],t,[],[],[],tol,g,rho,I2htoh,Rhto2h,i1,i2);
N = sqrt(length(u)/2); 
u1=u(1:length(u)/2);u2=u(length(u)/2+1:end);
U1ref=reshape(full(u1),N,N);U2ref=reshape(full(u2),N,N);
save('Uref.mat','U1ref','U2ref');