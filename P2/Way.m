function [xm,flag,relres,iter,resvec] = Way(allA,b,restart,tol,maxit,x0,I2htoh,Rhto2h,l,MGpre,i1,i2,nn)
    if nn==1
        [xm,flag,relres,iter,resvec] = MG(allA,b,tol,x0,I2htoh,Rhto2h,l,MGpre,i1,i2);
    elseif nn==2
        [xm,flag,relres,iter,resvec] = MGgmres(allA,b,restart,tol,maxit,x0,I2htoh,Rhto2h,l,MGpre,i1,i2);
    end
end