clear;
% clc;
% 从最粗网格层迭代
load('A-B.mat');    load('Uref.mat');
% 参数
E =1;  % Young modulus
k = 0.3;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operato r i ~= j;
tol = 1e-6; g = 0.06;   rho=2;   i1=2;  i2=2;
uerror = zeros(8,3); href=1/(length(U1ref)-1);
[noderef,elemref,Db,Nb_L,Nb_R,Cb]=generation_elem_node(href);
u0=[]; Lambda=[]; restart=[]; M=[]; allres1=cell(6,4);
maxN = 8;
minN = 3;
ff = ["Vcycle","Fcycle","Wcycle"];
for nn=1:1
    for j=1:1
        way = ff(j);
for kk=5:5
    u0=[]; Lambda=[]; restart=[]; M=[];
    tic;
    for l=minN:kk
        h=1/(2^(l+1));
        A=allA{l};  B=allB{l}; % 其中B={B{1},B{2}};
        if l > minN
            u0 = I2htoh{l}*u;
            [m,n]=size (I2htoh{l});
            I=I2htoh{l}(1:m/2,1:n/2);
            Lambda=I*Lambda;
        end
%         u0=[]; Lambda=[];   %注意有无初始值
        [u,Lambda,resvec,res,iter,Niter] = Friction_Solve(K1,K2,allA,B,h,u0,l,Lambda,restart,M,tol,g,rho,I2htoh,Rhto2h,i1,i2,way,nn);
    end
    t = toc;
    loc=(nn-1)*3+j; loc2=kk-4;
    allres1{loc,loc2}=resvec;
    fprintf("网格1/ %d的求解时间:%.4f(s),平均内迭代次数：%d\n",1/h,t,Niter);
%     [node,elem,Db,Nb_U,Nb_R,Cb]=generation_elem_node(h);
%     subplot(1,2,2);
%     Cbpoint=unique(Cb);
%     x = 0:h:1;
%     y = u(Cbpoint);
%     plot(x',y);hold on;
%     
%     if kk==5
%         [X,Y]=meshgrid(0:h:1);
%         N = sqrt(length(u)/2); 
%         u1=u(1:length(u)/2);u2=u(length(u)/2+1:end);
%         U1=reshape(full(u1),N,N);U2=reshape(full(u2),N,N);
%         subplot(1,2,1);
%         trisurf(elem,X+U1,Y+U2,0*U1);
%         hold on;plot([-0.1 1.1],[-g -g],'r');
%         hold on;plot([-0.1 1.1],[-0.3 -0.3],'r');
%         hold on;plot([-0.1 1.1],[0 0],'r');
%         hold on;plot([-0.1 -0.1],[0 -0.3],'r');
%         hold on;plot([1.1 1.1],[0 -0.3],'r');
%         axis off;
%     end
end
    end
end

