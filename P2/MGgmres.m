function [xm,flag,relres,iter,resvec] = MGgmres(allA,b,restart,tol,maxit,x0,I2htoh,Rhto2h,l,MGpre,i1,i2)
%% Input & Ouput parameter
%           ��GMRES�㷨������Ԥ����
% allA:     ��������ϵ��ܸپ���ϵ������--- ���ϡ�财��
% b:        ��ϸ����㣨�����㣩�ϵ��Ҷ�����
% restart:  restart��gmres������Ҳ������ν��GMRES(m)�㷨
% tol:      �������ȣ�||b-Ax_n||/||b||��
% maxit:    �������������������Ļ�һ��Ϊlength(b);��
% x0:       ��ʼֵ
% I2htoh:   �������Ĳ�ֵ����
% Ihto2h:   �����������ƾ���
% -------------------------------------------
% xm:        ���������ֵ��
% flag:     �ж��Ƿ��������������tol��ָ�ꣻflag = 0������������tol������û������
% relres:   ���һ�ε�������Բв�
% iter:     ������������ʽΪ[outiter,initer]
% resvec:   ÿ�ε�������Բв�
% MGpre:    �������͵Ķ�������ȥ����Ԥ������������������

if isempty(restart)
    restart = length(b);
end

A = allA{l};
if isempty(x0)
    x0 = zeros(size(b));
end
r0 = b-A*x0;   
beta = norm(r0);    
resvec = [norm(r0)];
normb = norm(b);
if beta/normb < tol     % ��ʼֵ�ͺܺã����㾫��
    flag = 0;
    iter = [0,0];
    relres =  beta/normb ;
    resvec = beta;
    xm = x0;
    return;
end

v(:,1) = r0/beta;
kx(1,:) = beta*1;    % beta*e1;
for outit = 1:maxit
    r0 = b-A*x0;    
    beta = norm(r0);  
    v(:,1) = r0/beta;
    kx(1,:) = beta*1;    % beta*e1;
    h = []; c=[];s =[]; z= [];
    for j = 1:restart
        % w = M\v(:,j);
        if MGpre == 'Vcycle'
            [z(:,j)] = V_cycle(allA,v(:,j),[],I2htoh,Rhto2h,i1,i2,l,1e-6);
        elseif MGpre == 'Fcycle'
            [z(:,j)] = F_cycle(allA,v(:,j),[],I2htoh,Rhto2h,i1,i2,l,1e-6,l);
        elseif MGpre == 'Wcycle'
            [z(:,j)] = W_cycle(allA,v(:,j),[],I2htoh,Rhto2h,i1,i2,l,1e-6,l);
        else
            fprintf('����Ԥ����������������Vcycle,Fcycle,Wcycle');
            quit;
        end
        w = A*z(:,j);
    
        for i = 1:j
            h(i,j) = v(:,i)'*w;
            w = w-h(i,j)*v(:,i);
        end
        h(j+1,j) = norm(w);
    
        if h(j+1,j) == 0
            break;
        end
    
        v(:,j+1) = w/h(j+1,j);  % ����µ�������
    
        for i = 1:j-1
            t = h(i,j);
            h(i,j) = c(i,:)* h(i,j) + s(i,:)*h(i+1,j);
            h(i+1,j) = -s(i,:)*t + c(i,:)*h(i+1,j);
        end
    
        if abs(h(j,j)) > abs(h(j+1,j))
            tau = h(j+1,j)/h(j,j); c(j,:) = 1/sqrt(1+tau^2); s(j,:) = c(j,:)*tau;
        else
            tau = h(j,j)/h(j+1,j); s(j,:) = 1/sqrt(1+tau^2); c(j,:) = s(j,:)*tau;
        end
    
        h(j,j) = c(j,:)*h(j,j) + s(j,:)*h(j+1,j);
        h(j+1,j) = 0;
        
        t = kx(j,:);
        kx(j,:) = c(j,:) * kx(j,:);
        kx(j+1,:) = -s(j,:) * t;
        
        relres = norm(kx(j+1,:)) / normb;
        resvec = [resvec;norm(kx(j+1,:))];
        if relres < tol
            break;
        end
    end
    m = j;
    ym = h(1:m,1:m)\kx(1:m,:);
    xm = x0+z(:,1:m)*ym;

    if relres < tol
        break;
    end
    x0 = xm;
end
iter = [outit,m];
if relres < tol
%     fprintf(['use ' MGpre ' get a convergent solution.\n']);
    flag = 0;
else
    fprintf(['use ' MGpre ' not get a convergent solution.\n']);
    flag = 1;
end
resvec = resvec/norm(b);
end 