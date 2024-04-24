function [uh] = W_cycle(allA,b,u0,I2htoh,Rhto2h,i1,i2,l,tol,lmax)
%-------------------------------------------------------
    % allA:     ���������ܸپ���Ah
    % b:        ��ǰ�����Omega�ϻ��ֵ��Ҷ�����
    % u0:       ��ǰ������ϵĳ�ʼֵ
    % I2htoh:   �������Ĳ�ֵ����
    % Rhto2h:   ����������������
    % i1:       ǰ�⻬����
    % i2:       ��⻬����
    % l:        ��ǰ�� Ĭ���������Ϊh=1/4����ʱl=1
    % tol:      ��Ծ��ȣ�Ĭ��1e-6
%--------------------------------------------------------
%     [uh,~,~,~,~] = gmres(allA{l},b,[],tol,i1,[],[],u0);
    [uh,~,~] = GS(allA{l},b,i1,u0);
%     fprintf("level:%d\n",l);
    r2h = Rhto2h{l} * (b - allA{l} * uh);  % ���Ƶ���������
   
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
    r2h = Rhto2h{l} * (b - allA{l} * uh);  % ���Ƶ���������
   
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