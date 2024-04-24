clear;clc;
%----------------------------------------------
% ���룺����ģ�� E �� ���ɱ� k
% ��������������ܸپ���A����Omega�ϵ��Ҷ�����
% �������� h=1/4 ������1������
% �ڵ�Ԫ�ϵ��ܸپ���A������ͬ��
% �������� ע�⣺Dirichlet condition �Ƿ�����ȥ��������������
%----------------------------------------------

E = 1;  % Young modulus
k = 0.3;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operator i ~= j;
allA=cell(8,1); allB=cell(8,1);  
I2htoh=cell(8,1); Rhto2h=cell(8,1);
for t=1:1:8
    fprintf("�ڼ���1/%d�����ϵ������;",2^(t+1));
    h=1/(2^(t+1));
    [node,elem,Db,Nb_U,Nb_R,Cb]=generation_elem_node(h);

    N = size(node,1);  %�ڵ����
    NT = size(elem,1); %��Ԫ����
    freeNodes = setdiff(1:N,unique(Db));  %Dirichlet points is certain.
    Cbpoint=unique(Cb);nCb=length(Cbpoint);
    A11 = sparse(N,N);   %�ܸվ���
    A12 = A11;A21 = A11;A22 = A11;
    B1 = full(sparse(N,1));   %�Ҷ���
    B2 = sparse(N,1); 
    % �ܸվ���
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
    % Gamma_2Ӧ������
    if size(Nb_U,1) ~= 0
        for i = 1:1:size(Nb_U,1)
            curNb_L = Nb_U(i,:);
            length4curedge = norm(node(curNb_L(1),:)-node(curNb_L(2),:));
            v = node(curNb_L,:); 
            [GL1,GL2] = Generate_gR(K1,K2,v); % ��߽߱�
            B1(curNb_L) = B1(curNb_L)+(1/6)*length4curedge*GL1';
            B2(curNb_L) = B2(curNb_L)+(1/6)*length4curedge*GL2'; 
        end
    end
%     if size(Nb_R,1) ~= 0
%         for i = 1:1:size(Nb_R,1)
%             curNb_R = Nb_R(i,:);
%             length4curedge = norm(node(curNb_R(1),:)-node(curNb_R(2),:));
%             v = node(curNb_R,:); 
%             [GR1,GR2] = Generate_gR(K1,K2,v); % �ұ߽߱�
%             B1(curNb_R) = B1(curNb_R)+(1/6)*length4curedge*GR1';
%             B2(curNb_R) = B2(curNb_R)+(1/6)*length4curedge*GR2'; 
%         end
%     end
    z = unique(Cb);  A11(z,:)=2*A11(z,:);A12(z,:)=2*A12(z,:);A21(z,:)=2*A21(z,:);A22(z,:)=2*A22(z,:);
    z = unique(Nb_U);A11(z,:)=2*A11(z,:);A12(z,:)=2*A12(z,:);A21(z,:)=2*A21(z,:);A22(z,:)=2*A22(z,:);
    z = unique(Nb_R);A11(z,:)=2*A11(z,:);A12(z,:)=2*A12(z,:);A21(z,:)=2*A21(z,:);A22(z,:)=2*A22(z,:);
    % Gamma_1�߽�����
    if size(Db,1) ~= 0
        DbNodes = unique(Db);  %Dirichlet points is certain.
        nDb = length(DbNodes);
        A11(DbNodes,:)=0;A12(DbNodes,:)=0;
        A21(DbNodes,:)=0;A22(DbNodes,:)=0;
        A11(DbNodes,DbNodes)=speye(nDb);A22(DbNodes,DbNodes)=speye(nDb);
        A21(Cbpoint,:)=0;A22(Cbpoint,:)=0;  % ��֤gamma3��u_{nu}=0
        A22(Cbpoint,Cbpoint)=speye(nCb);    % ��֤gamma3��u_{nu}=0
        A  = [A11,A12;A21,A22];
        B1(DbNodes)=sparse(nDb,1);B2(DbNodes)=sparse(nDb,1);
        B2(Cbpoint)=sparse(nCb,1); % ��֤gamma3��u_{nu}=0
        B = [B1;B2];
    end
    A = A./(h^2);
    % �洢�նȾ���
    allA{t}=A;  allB{t}={B1,B2};
    
    % ��ֵ���Ӻ���������
    [I,R]=IR_FW(2*h);
    I2htoh{t}=I;    Rhto2h{t}=R;
    fprintf("1/%d�����ϵ�����������ϡ�\n",2^(t+1));
    pause(0.5);
end
save('A-B.mat','allA','allB','I2htoh','Rhto2h');

clear;
clc;
% �������������
load('A-B.mat');  
% ����
E =1;  % Young modulus
k = 0.3;   % Poisson ratio of the material
K1 = E*k/(1-k^2);    % Coefficient of operator i = j ;
K2 = E/(1+k);        % Coefficient of operator i ~= j; 
tol = 1e-6; g = 0.06;   rho=2;   i1=2;  i2=2;
t = 8;
h = 1/(2^(t+1));
A=allA{t};  B=allB{t}; % ����B={B{1},B{2}};
[u,Lambda,resvec,res,iter,Niter] = Friction_Solve(K1,K2,allA,B,h,[],t,[],[],[],tol,g,rho,I2htoh,Rhto2h,i1,i2);
N = sqrt(length(u)/2); 
u1=u(1:length(u)/2);u2=u(length(u)/2+1:end);
U1ref=reshape(full(u1),N,N);U2ref=reshape(full(u2),N,N);
save('Uref.mat','U1ref','U2ref');