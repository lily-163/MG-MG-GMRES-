function [u,Lambda,resvec,res,iter,Niter] = Friction_Solve(K1,K2,allA,B,h,u0,l,Lambda,restart,M,tol,g,rho,I2htoh,Rhto2h,i1,i2,way,nn)
%% Input:
%   K1,K2:  ����ϵ��
%   A:      ��ǰ���ܸپ���
%   B:      ��ǰ���Ҷ�����
%   u0:     ��ǰ���ʼֵ
%   Lambda: ��ǰ��ĳ�ʼLagrange����
%   restart:        gmres�㷨����������
%   M:              gmreesԤ������
%   tol:            Lagrange������������ֵ
%   g:              ����Լ��������, ��g=0.02m
%   rho:            Lagrange���Ӹ��µĲ�������

%% Output:
%   u:      ��ǰ�����������ֵ��
%   Lambda: ��ǰ�����������Lagrange����
%   resvec: ���һ��Lagrange�������������Au=Bʱ�Ĳв�����
%   res:    ��������Lagrange�������������������

%% ����
A = allA{l};
maxiter=1000;resvec=[];res=[];last_res = [];
[node,elem,Db,Nb_U,Nb_R,Cb]=generation_elem_node(h);
N = size(node,1);   % �ڵ����
NT = size(elem,1);  % ��Ԫ����
Niter=[];

Cbpoint=unique(Cb); nCb=length(Cbpoint);
DbNodes = unique(Db);  %Dirichlet points is certain.
nDb = length(DbNodes);
B1=B{1};B2=B{2};
z = unique(Cb);  B2(z)=2*B2(z);
z = unique(Nb_U);B2(z)=2*B2(z);
z = unique(Nb_R);B2(z)=2*B2(z);
B2 = B2./h^2;
% gamma3Ħ���߽�;
if size(Cb,1) ~= 0   
    for iter=1:1:maxiter 
        % ����lambda����;
        if iter==1 
            if isempty(Lambda)
                Lambda=sparse(N,1);
                Lambda(Cbpoint) = sparse(nCb,1);
            elseif ~isempty(u0)
                for i=1:1:length(Cbpoint)
                    loc = Cbpoint(i);
                    lambda2=Lambda(loc)+rho*u0(loc);
                    Lambda(loc)=cal_tau(lambda2,u0(loc)); % ��ʾu2
                end
            end
        else
            Lambda0 = Lambda;
            for i=1:1:length(Cbpoint)
                loc = Cbpoint(i);
                lambda2=Lambda(loc)+rho*u(loc);
                lambda2=cal_tau(lambda2,u(loc)); % ��ʾu2
                Lambda(loc) = 0.5*Lambda(loc)+0.5*lambda2;
            end
        end
        % ����(lambda,v_nu);
        CB1=B1;
        for i=1:1:size(Cb,1)
            curCb=Cb(i,:);
            length4curedge = norm(node(curCb(1),:)-node(curCb(2),:));
            F = Generate_tau(Lambda(curCb));
            CB1(curCb) = CB1(curCb) - (1/6)*length4curedge*F';
        end
        % ��װB
        z = unique(Cb);  CB1(z)=2*CB1(z);
        z = unique(Nb_U);CB1(z)=2*CB1(z);
        z = unique(Nb_R);CB1(z)=2*CB1(z);
        CB1 = CB1./h^2;
        if size(Db,1) ~= 0
            CB1(DbNodes)=sparse(nDb,1);B2(DbNodes)=sparse(nDb,1);
            B = [CB1;B2];
        end
        % ��������B1
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
%                 fprintf('%d��',nnn);
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
