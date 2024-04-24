function [uh] = W_cycle(allA,b,u0,I2htoh,Rhto2h,i1,i2,l,tol,lmax)
%-------------------------------------------------------
    % allA:     各网格层的总纲矩阵Ah
    % b:        当前网格层Omega上积分的右端向量
    % u0:       当前网格层上的初始值
    % I2htoh:   各网格层的插值算子
    % Rhto2h:   各网格层的限制算子
    % i1:       前光滑次数
    % i2:       后光滑次数
    % l:        当前层 默认最粗网格为h=1/4，此时l=1
    % tol:      相对精度，默认1e-6
%--------------------------------------------------------
%     [uh,~,~,~,~] = gmres(allA{l},b,[],tol,i1,[],[],u0);
    [uh,~,~] = GS(allA{l},b,i1,u0);
%     fprintf("level:%d\n",l);
    r2h = Rhto2h{l} * (b - allA{l} * uh);  % 限制到粗网格上
   
    if l-1==1
        e2h = allA{l-1}\r2h;
%         fprintf("level:%d\n",l-1);
    else
        e2h = W_cycle(allA,r2h,[],I2htoh,Rhto2h,i1,i2,l-1,tol,lmax);
    end
    uh=uh+I2htoh{l} * e2h;
%     [uh,~,~,~,~] = gmres(allA{l},b,[],tol,i2,[],[],uh);
    [uh,~,~] = GS(allA{l},b,i2,uh);
%     fprintf("level:%d\n",l);
    if l ~= lmax
    r2h = Rhto2h{l} * (b - allA{l} * uh);  % 限制到粗网格上
   
    if l-1==1
        e2h = allA{l-1}\r2h;
%         fprintf("level:%d\n",l-1);
    else
        e2h = W_cycle(allA,r2h,[],I2htoh,Rhto2h,i1,i2,l-1,tol,lmax);
    end

    uh=uh+I2htoh{l} * e2h;
%     [uh,~,~,~,~] = gmres(allA{l},b,[],tol,i2,[],[],uh);
    [uh,~,~] = GS(allA{l},b,i1,uh);
%     fprintf("level:%d\n",l);
    end 
end