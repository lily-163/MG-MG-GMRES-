function [x0,flag,relres,iter,resvec] = MG(allA,b,tol,x0,I2htoh,Rhto2h,l,MGpre,i1,i2)
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
relres=[];
resvec = [];
if isempty(x0)
    r = b;
else
    r=  b - allA{l}*x0;
end
relres = norm(r)/norm(b);
resvec=[resvec;relres];
for iter = 1:1000
    if MGpre == 'Vcycle'
        x0 = V_cycle(allA,b,x0,I2htoh,Rhto2h,i1,i2,l,tol);
    elseif MGpre == 'Fcycle'
        x0 = F_cycle(allA,b,x0,I2htoh,Rhto2h,i1,i2,l,tol,l);
    elseif MGpre == 'Wcycle'
        x0 = W_cycle(allA,b,x0,I2htoh,Rhto2h,i1,i2,l,tol,l);
    else
        fprintf('����MG������ʽ����������Vcycle,Fcycle,Wcycle');
    end
    r = b - allA{l}*x0;
    
    relres = norm(r)/norm(b);
    resvec = [resvec;relres];
    if relres < tol
        flag=0;
        break
    end
    flag=1;
end















