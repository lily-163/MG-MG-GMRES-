clear;
syms x y;
% u1 = sin(pi*x)*exp(x)*y*(1-y);
% u2 = sin(2*pi*y)*log(1+x)*(1-x);

u1 = y^2*(1 - y)*exp(x);
u2 = y*(y^2-1)*cos(pi*x);

du1dxx = simplify(diff(diff(u1,x),x));
du1dxy = simplify(diff(diff(u1,x),y));
du1dyy = simplify(diff(diff(u1,y),y));
du1dyx = simplify(diff(diff(u1,y),x));

du2dxx = simplify(diff(diff(u2,x),x));
du2dxy = simplify(diff(diff(u2,x),y));
du2dyy = simplify(diff(diff(u2,y),y));
du2dyx = simplify(diff(diff(u2,y),x));

du1dx = simplify(diff(u1,x));
du1dy = simplify(diff(u1,y));
du2dx = simplify(diff(u2,x));
du2dy = simplify(diff(u2,y));