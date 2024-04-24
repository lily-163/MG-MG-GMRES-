function [x0,flag,relres,iter,resvec] = MG(allA,b,tol,x0,I2htoh,Rhto2h,l,MGpre,i1,i2)
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
        fprintf('输入MG方法格式有误，请输入Vcycle,Fcycle,Wcycle');
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















