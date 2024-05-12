clear;
clc;
close all;

%% -------------------
E =2;  % Young modulus
k = 0.4;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operator i ~= j; 
maxiter=1000;
rho=2.5;   g=0.07;  
tol=1e-6;
uerror = []; res=[];Niter=[];
maxgrid = 8;   href=1/(2^(maxgrid+1));
% load('A-B.mat');load('Uref.mat');
h = 1/(2^(maxgrid+1));
[noderef,elemref,Db,Nb_L,Nb_R,Cb]=generation_elem_node(h);

for t =5:-1:5
    h=1/(2^(t+1));
    [node,elem,Db,Nb_L,Nb_R,Cb]=generation_elem_node(h);

    N = size(node,1);  %节点个数
    NT = size(elem,1); %单元个数
    freeNodes = setdiff(1:N,unique(Db));  %Dirichlet points is certain.
    Cbpoint=unique(Cb);nCb=length(Cbpoint);
    A11 = sparse(N,N);   %总刚矩阵
    A12 = A11;A21 = A11;A22 = A11;
    B1 = sparse(N,1);   %右端项
    B2 = sparse(N,1); 

    % 总刚矩阵----标准单元上
    a1=0; b1=1; c1=0;
    a2=-1;b2=-1;c2=1;
    a3=1; b3=0; c3=0;coef=[a1 b1 c1;a2 b2 c2;a3 b3 c3];
    dxx = 0.5*[a1*a1 a2*a1 a3*a1;a1*a2 a2*a2 a3*a2;a1*a3 a2*a3 a3*a3];
    dyy = 0.5*[b1*b1 b2*b1 b3*b1;b1*b2 b2*b2 b3*b2;b1*b3 b2*b3 b3*b3];
    dxy = 0.5*[a1*b1 a2*b1 a3*b1;a1*b2 a2*b2 a3*b2;a1*b3 a2*b3 a3*b3];
    dyx = dxy';
    S=h*h/2;

    cA11=(K1+K2)*dxx+0.5*K2*dyy;
    cA12=K1*dyx+0.5*K2*dxy;
    cA21=K1*dxy+0.5*K2*dyx;
    cA22=(K1+K2)*dyy+0.5*K2*dxx;
    for i=1:NT
        curelem = elem(i,:);
        v = node(curelem,:);  
        A11(curelem,curelem)=A11(curelem,curelem)+cA11;
        A12(curelem,curelem)=A12(curelem,curelem)+cA12;
        A21(curelem,curelem)=A21(curelem,curelem)+cA21;
        A22(curelem,curelem)=A22(curelem,curelem)+cA22;

        B1(curelem)= B1(curelem)+Generate_f01(K1,K2,v(:,1),v(:,2))*S/3;
        B2(curelem)= B2(curelem)+Generate_f02(K1,K2,v(:,1),v(:,2))*S/3;
    end
% Gamma_2应力条件
if size(Nb_L,1) ~= 0
    for i = 1:1:size(Nb_L,1)
        curNb_L = Nb_L(i,:);
        length4curedge = norm(node(curNb_L(1),:)-node(curNb_L(2),:));
        v = node(curNb_L,:); 
        [GL1,GL2] = Generate_gL(K1,K2,v); % 左边边界
        B1(curNb_L) = B1(curNb_L)+(1/6)*length4curedge*GL1';
        B2(curNb_L) = B2(curNb_L)+(1/6)*length4curedge*GL2'; 
    end
end
if size(Nb_R,1) ~= 0
    for i = 1:1:size(Nb_R,1)
        curNb_R = Nb_R(i,:);
        length4curedge = norm(node(curNb_R(1),:)-node(curNb_R(2),:));
        v = node(curNb_R,:); 
        [GR1,GR2] = Generate_gR(K1,K2,v); % 右边边界
        B1(curNb_R) = B1(curNb_R)+(1/6)*length4curedge*GR1';
        B2(curNb_R) = B2(curNb_R)+(1/6)*length4curedge*GR2'; 
    end
end
% Gamma_1边界条件
if size(Db,1) ~= 0
    DbNodes = unique(Db);  %Dirichlet points is certain.
    nDb = length(DbNodes);
    A11(DbNodes,:)=0;A12(DbNodes,:)=0;
    A21(DbNodes,:)=0;A22(DbNodes,:)=0;
    A11(DbNodes,DbNodes)=speye(nDb);A22(DbNodes,DbNodes)=speye(nDb);
    A  = [A11,A12;A21,A22];
    B1(DbNodes)=sparse(nDb,1);B2(DbNodes)=sparse(nDb,1);
    B = [B1;B2];
end
% gamma3摩擦边界;
if size(Cb,1) ~= 0   
    for iter=1:1:maxiter 
        % 更新lambda乘子;
        if iter==1
            % tem=1;
            Lambda=sparse(N,1);
            Lambda(Cbpoint) = sparse(nCb,1) ;
        else
            Lambda0 = Lambda;
            for i=1:1:length(Cbpoint)
                loc = Cbpoint(i);
                flag1=nu(-u(N+loc)); % 表示u_{nu} = -u2
                flag2=Lambda(loc)+rho*((-u(N+loc))-g);
                % Lambda(loc) = 0.5*Lambda(loc)+0.5*flag1;     % pro2
                Lambda(loc)=max(flag1,flag2);   % pro3
            end
        end
        % 计算(lambda,v_nu);
        CB2=B2;
        for i=1:1:size(Cb,1)
            curCb=Cb(i,:);
            length4curedge = norm(node(curCb(1),:)-node(curCb(2),:));
            v=node(curCb,:);
            F = Generate_nu(Lambda(curCb));
            CB2(curCb) = CB2(curCb) - (1/6)*length4curedge*F';
        end
        % 组装B
        if size(Db,1) ~= 0
            B1(DbNodes)=sparse(nDb,1);CB2(DbNodes)=sparse(nDb,1);
            B = [B1;CB2];
        end
        % 迭代收敛B1
        if iter==1
            u = A\B;
            % [u,flag,relres,it0,resvec] = gmres(A,B,[],10^(-6),length(B),L,U,[]);
        else
            u_ = A\B;
            % [u_,flag,relres,it0,resvec] = gmres(A,B,[],10^(-6),length(B),L,U,[]);
            temp = norm(u_-u)/norm(u);
            temp1 = norm(Lambda(Cbpoint)-Lambda0(Cbpoint))/norm(Lambda0(Cbpoint));
            temp
            res=[res;temp];
            if (temp < tol) % && (temp1 < tol)
                break;
            else
                u=u_;  
            end
        end
    end
else
    u_ = A\B; 
end  % gamma3摩擦边界;

u=u_;
while h>href
    [I,R] = IR_FW(h);
    u = I*u;
    h=h/2;
end
% 返回二维数据

N = sqrt(length(u)/2); 
u1=u(1:length(u)/2);u2=u(length(u)/2+1:end);
U1=reshape(full(u1),N,N);U2=reshape(full(u2),N,N);
disp(iter)

if t==maxgrid
    U1ref=U1;
    U2ref=U2;
    ss = 1;
    elemref=elem;
    noderef=node;
else
    err = Q_error(K1,K2,U1,U2,U1ref,U2ref,elemref,noderef);
    uerror = [uerror,err];
    Niter=[Niter,iter];
    ss = ss+1;
end
end

%% 计算Q范数误差
t = length(uerror);
for i = 1:1:t-1
    orderu = uerror(i+1)/uerror(i);
    fprintf('1/%d层与1/%d层间收敛阶:%.4f\n',2^(t-i+2),2^(t-i+1),log2(orderu));
end   
U1ref=U1;
U2ref=U2;
plot_sigmanu
