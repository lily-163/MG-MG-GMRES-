function [u,Lambda,resvec,res,iter,Niter] = Friction_Solve(K1,K2,allA,B,h,u0,l,Lambda,restart,M,tol,g,rho,I2htoh,Rhto2h,i1,i2,way,nn)
%% Input:
%   K1,K2:  方程系数
%   A:      当前层总纲矩阵
%   B:      当前层右端向量
%   u0:     当前层初始值
%   Lambda: 当前层的初始Lagrange乘子
%   restart:        gmres算法的重启限制
%   M:              gmrees预处理子
%   tol:            Lagrange迭代的收敛阈值
%   g:              法线约束的限制, 如g=0.02m
%   rho:            Lagrange乘子更新的步长控制

%% Output:
%   u:      当前网格层问题数值解
%   Lambda: 当前网格层收敛的Lagrange乘子
%   resvec: 最后一次Lagrange迭代过程中求解Au=B时的残差向量
%   res:    相邻两次Lagrange迭代的收敛解的相对误差

%% 参数
A = allA{l};
maxiter=1000;resvec=[];res=[];last_res = [];
[node,elem,Db,Nb_U,Nb_R,Cb]=generation_elem_node(h);
N = size(node,1);   % 节点个数
NT = size(elem,1);  % 单元个数
Niter=[];

Cbpoint=unique(Cb); nCb=length(Cbpoint);
DbNodes = unique(Db);  %Dirichlet points is certain.
nDb = length(DbNodes);
B1=B{1};B2=B{2};
z = unique(Cb);  B2(z)=2*B2(z);
z = unique(Nb_U);B2(z)=2*B2(z);
z = unique(Nb_R);B2(z)=2*B2(z);
B2 = B2./h^2;
% gamma3摩擦边界;
if size(Cb,1) ~= 0   
    for iter=1:1:maxiter 
        % 更新lambda乘子;
        if iter==1 
            if isempty(Lambda)
                Lambda=sparse(N,1);
                Lambda(Cbpoint) = sparse(nCb,1);
            elseif ~isempty(u0)
                for i=1:1:length(Cbpoint)
                    loc = Cbpoint(i);
                    lambda2=Lambda(loc)+rho*u0(loc);
                    Lambda(loc)=cal_tau(lambda2,u0(loc)); % 表示u2
                end
            end
        else
            Lambda0 = Lambda;
            for i=1:1:length(Cbpoint)
                loc = Cbpoint(i);
                lambda2=Lambda(loc)+rho*u(loc);
                lambda2=cal_tau(lambda2,u(loc)); % 表示u2
                Lambda(loc) = 0.5*Lambda(loc)+0.5*lambda2;
            end
        end
        % 计算(lambda,v_nu);
        CB1=B1;
        for i=1:1:size(Cb,1)
            curCb=Cb(i,:);
            length4curedge = norm(node(curCb(1),:)-node(curCb(2),:));
            F = Generate_tau(Lambda(curCb));
            CB1(curCb) = CB1(curCb) - (1/6)*length4curedge*F';
        end
        % 组装B
        z = unique(Cb);  CB1(z)=2*CB1(z);
        z = unique(Nb_U);CB1(z)=2*CB1(z);
        z = unique(Nb_R);CB1(z)=2*CB1(z);
        CB1 = CB1./h^2;
        if size(Db,1) ~= 0
            CB1(DbNodes)=sparse(nDb,1);B2(DbNodes)=sparse(nDb,1);
            B = [CB1;B2];
        end
        % 迭代收敛B1
        if iter==1
%             u = A\B;
            if h ==1/8
                u=A\B;
            else
                [u,niter,~] = GS(A,B,100000,u0);
%                 [u,~,~,niter,~] = MG(allA,B,tol,u0,I2htoh,Rhto2h,l,'Wcycle',i1,i2);
%                 [u,~,~,niter,~] = MGgmres(allA,B,[],tol,length(B),u0,I2htoh,Rhto2h,l,'Vcycle',i1,i2);
%                 [u,~,~,niter,~] = Way(allA,B,[],tol,length(B),u0,I2htoh,Rhto2h,l,way,i1,i2,nn);
                Niter=[Niter;niter(end)];
            end 
        else
%             u_ = A\B;
            if h==1/8
                u_=A\B;
            else
                [u_,niter,~] = GS(A,B,100000,u0);
%                 [u_,~,~,niter,last_res] = MG(allA,B,tol,u0,I2htoh,Rhto2h,l,'Wcycle',i1,i2);
%                 [u_,~,~,niter,last_res] = MGgmres(allA,B,[],tol,length(B),u0,I2htoh,Rhto2h,l,'Vcycle',i1,i2);
%                 [u_,~,~,niter,last_res] = Way(allA,B,[],tol,length(B),u0,I2htoh,Rhto2h,l,way,i1,i2,nn);
                Niter=[Niter;niter(end)];
            end
            temp = norm(u_-u)/norm(u);
            temp1 = norm(Lambda-Lambda0)/norm(Lambda0);
            res=[res;temp];
            if temp < tol  %&& (temp1 < tol)
%                 [~,nnn,~] = GS(A,B,100000,u0);
%                 fprintf('%d次',nnn);
                break;
            else
                u=u_;
            end
        end
    end
end  
u=u_;
if isempty(Niter)
    Niter=1;
else
    Niter = mean(Niter);
end
resvec = last_res;
end
