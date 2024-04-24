N = size(U1ref,1);
h=1/(N-1);
[node,elem,Db,Nb_L,Nb_R,Cb]=generation_elem_node(h);
g=0.06;
[X,Y]=meshgrid(0:h:1);

F1=scatteredInterpolant(X(:),Y(:),U1ref(:),'linear');
F2=scatteredInterpolant(X(:),Y(:),U2ref(:),'linear');

x=X(1,:)';
y=Y(:,1);

N=size(U1ref,1); 
sigmanu_left=sparse(N,2);
sigmanu_right=sparse(N,2);
sigmanu_up=sparse(N,2);
sigmanu_down=sparse(N,2);

% sigmanu_left(:,1)=-(K1+K2)*(F1(0*x+h,y)-F1(0*x-h,y))/(2*h)-K1*(F2(0*x,y+h)-F2(0*x,y-h))/(2*h);
% sigmanu_left(:,2)=-0.5*K2*(F1(0*x,y+h)-F1(0*x,y-h))/(2*h)-0.5*K2*(F2(0*x+h,y)-F2(0*x-h,y))/(2*h);
% figure(1);
% quiver(0*x,y,sigmanu_left(:,1),sigmanu_left(:,2),0);
% title('sigma_{nu} left')
% 
% one = ones(N,1);
% sigmanu_right(:,1)=(K1+K2)*(F1(one+h,y)-F1(one-h,y))/(2*h)+K1*(F2(one,y+h)-F2(one,y-h))/(2*h);
% sigmanu_right(:,2)=0.5*K2*(F1(one,y+h)-F1(one,y-h))/(2*h)+0.5*K2*(F2(one+h,y)-F2(one-h,y))/(2*h);
% figure(2);
% quiver(one,y,sigmanu_right(:,1),sigmanu_right(:,2),0);
% title('sigma_{nu} right')

% sigmanu_up(:,1)=0.5*K2*(F1(x,one+h)-F1(x,one-h))/(2*h)+0.5*K2*(F2(x+h,one)-F2(x-h,one))/(2*h);
% sigmanu_up(:,2)=(K1+K2)*(F2(x,one+h)-F2(x,one-h))/(2*h)+K1*(F1(x+h,one)-F1(x-h,one))/(2*h);
% figure(3);
% quiver(x',one',sigmanu_up(:,1)',sigmanu_up(:,2)',0);
% title('sigma_{nu} up')

sigmanu_down(:,1)=-0.5*K2*(F1(x,0*y+h)-F1(x,0*y-h))/(2*h)-0.5*K2*(F2(x+h,0*y)-F2(x-h,0*y))/(2*h);
sigmanu_down(:,2)=-(K1+K2)*(F2(x,0*y+h)-F2(x,0*y-h))/(2*h)-K1*(F1(x+h,0*y)-F1(x-h,0*y))/(2*h);
figure(4);
quiver(x',0*y',-sigmanu_down(:,1)',-0*sigmanu_down(:,2)');
title('sigma_{nu} down(等比例放缩后)')

% lll=Lambda(Cbpoint);
figure(1);
trisurf(elem,X+U1ref,Y+U2ref,0*U1ref);
hold on;plot([-0.1 1.1],[-g -g],'r');
hold on;plot([-0.1 1.1],[-0.3 -0.3],'r');
hold on;plot([-0.1 1.1],[0 0],'r');
hold on;plot([-0.1 -0.1],[0 -0.3],'r');
hold on;plot([1.1 1.1],[0 -0.3],'r');
% hold on;quiver(U1ref(1,:)+x',0*y',-sigmanu_down(:,1)',-sigmanu_down(:,2)',0,'b');
% axis([-0.2 1.2 -0.4 1.1]);
axis off;




