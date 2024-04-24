function [uh,res] = C_cycle(allA,allb,I2htoh,l,restart,tol)
%-------------------------------------------------------
    % 瀑布型多重网格 Cascadic Multigrid
    % allA:     各网格层的总纲矩阵Ah
    % allb:     各网格层Omega上积分的右端向量
    % u0:       当前网格层上的初始值
    % I2htoh:   各网格层的插值算子
    % restart:  几步完成一个重启
    % l:        当前层 默认最粗网格为h=1/4，此时l=1
    % tol:      相对精度，默认1e-6
%--------------------------------------------------------
    uh = [];
    for i=1:1:5
        [uh,~,~,~,res] = gmres(allA{i},[allb{i}{1};allb{i}{2}],[],tol,length(allb{i}{1})*2,[],[],uh);
        % u=allA{i}\allb{i};
        % norm(u-uh)/norm(u)
        if i==l
            break;
        else 
            uh=I2htoh{i+1} * uh;
            % u=allA{i+1}\allb{i+1};
            % norm(u-uh)/norm(u)
        end
    end
end