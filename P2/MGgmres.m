function [xm,flag,relres,iter,resvec] = MGgmres(allA,b,restart,tol,maxit,x0,I2htoh,Rhto2h,l,MGpre,i1,i2)
%% Input & Ouput parameter
%           对GMRES算法进行右预处理！
% allA:     各网格层上的总纲矩阵（系数矩阵）--- 最好稀疏储存
% b:        最细网格层（需求解层）上的右端向量
% restart:  restart次gmres后重启也就是所谓的GMRES(m)算法
% tol:      收敛精度（||b-Ax_n||/||b||）
% maxit:    最大迭代次数（不重启的话一般为length(b);）
% x0:       初始值
% I2htoh:   各网格层的插值矩阵
% Ihto2h:   各网格层的限制矩阵
% -------------------------------------------
% xm:        迭代后的数值解
% flag:     判断是否迭代收敛到精度tol的指标；flag = 0，代表收敛到tol，否则没有收敛
% relres:   最后一次迭代的相对残差
% iter:     迭代次数，格式为[outiter,initer]
% resvec:   每次迭代的相对残差
% MGpre:    哪种类型的多重网格去进行预处理，这里设置了三类

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
if beta/normb < tol     % 初始值就很好，满足精度
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
            fprintf('输入预处理方法有误，请输入Vcycle,Fcycle,Wcycle');
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
    
        v(:,j+1) = w/h(j+1,j);  % 获得新的正交基
    
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