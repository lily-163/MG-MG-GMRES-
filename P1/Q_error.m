function SS = Q_error(K1,K2,U1,U2,U1ref,U2ref,elem,node)
% ------------------------------------
% 已修改为最细网格上做误差估计，偏导都采用系数
%------------------------------
    N = sqrt(size(node,1));
    h = 1/(N-1);href=1/(size(U1ref,1)-1);
    Qref = 0;  
    EQ = 0;
%     [X,Y]=meshgrid(0:href:1);
%     F1=scatteredInterpolant(X(:),Y(:),U1ref(:),'linear');
%     F2=scatteredInterpolant(X(:),Y(:),U2ref(:),'linear');
    S = 1; % 分子分母都有S，就默认1；
    
    parfor i = 1:1:size(elem,1)
        v = node(elem(i,:),:);
        coef =[0 1 0;-1 -1 1;1 0 0]*(-1)^(i+1)/h;
%         coef=[a1 b1 c1;a2 b2 c2;a3 b3 c3];
        loc_yx = (N-1)*v+1;
        u1 = zeros(3,1); u2 = zeros(3,1);
        u1ref = zeros(3,1); u2ref = zeros(3,1);
        for j = 1:1:3
            u1ref(j)=U1ref(loc_yx(j,2),loc_yx(j,1));
            u2ref(j) = U2ref(loc_yx(j,2),loc_yx(j,1));
            u1(j) = U1(loc_yx(j,2),loc_yx(j,1));
            u2(j) = U2(loc_yx(j,2),loc_yx(j,1));
        end   
        
%         x = v(:,1); y = v(:,2);
%         x=[0.5*(x(1)+x(2));0.5*(x(2)+x(3));0.5*(x(3)+x(1))];
%         y=[0.5*(y(1)+y(2));0.5*(y(2)+y(3));0.5*(y(3)+y(1))];
%         du1refdx = ( F1(x+href,y)-F1(x-href,y) )/(2*href);
%         du1refdy = ( F1(x,y+href)-F1(x,y-href) )/(2*href);
%         du2refdx = ( F2(x+href,y)-F2(x-href,y) )/(2*href);
%         du2refdy = ( F2(x,y+href)-F2(x,y-href) )/(2*href);
        
        du1refdx = coef(:,1)'*u1ref;
        du1refdy = coef(:,2)'*u1ref;
        du2refdx = coef(:,1)'*u2ref;
        du2refdy = coef(:,2)'*u2ref;
        
        qref = ((K1+K2)*du1refdx+K1*du2refdy).*du1refdx +...
            ((K1+K2)*du2refdy+K1*du1refdx).*du2refdy+...
            2*K2*(du1refdy/2+du2refdx/2).^2;
        Qref =Qref + mean(qref)*S;
        %%%%%%%%%%%%%%
        edu1dx = du1refdx - coef(:,1)'*u1;
        edu1dy = du1refdy - coef(:,2)'*u1;
        edu2dx = du2refdx - coef(:,1)'*u2;
        edu2dy = du2refdy - coef(:,2)'*u2;
        
        eq = ((K1+K2)*edu1dx+K1*edu2dy).*edu1dx +...
            ((K1+K2)*edu2dy+K1*edu1dx).*edu2dy+...
            2*K2*(edu1dy/2+edu2dx/2).^2;
        EQ = EQ + mean(eq)*S;
        
    end
    Qref = sqrt(Qref/2);
    EQ = sqrt(EQ/2);
    SS = EQ/Qref;

end