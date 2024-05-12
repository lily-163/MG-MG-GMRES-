function [G1,G2] = Generate_gL(K1,K2,v)
    x1=v(1,1);y1=v(1,2);x2=v(2,1);y2=v(2,2);xmid=(x1+x2)/2;ymid=(y1+y2)/2;
    G1(1) = - g1(K1,K2,x1,y1)*1 - g1(K1,K2,x2,y2)*0 - 4*g1(K1,K2,xmid,ymid)*0.5 ; %f2(1)²¿·Ö
    G1(2) = - g1(K1,K2,x1,y1)*0 - g1(K1,K2,x2,y2)*1 - 4*g1(K1,K2,xmid,ymid)*0.5 ; 

    G2(1) = - g2(K1,K2,x1,y1)*1 - g2(K1,K2,x2,y2)*0 - 4*g2(K1,K2,xmid,ymid)*0.5 ; %f2(2)²¿·Ö
    G2(2) = - g2(K1,K2,x1,y1)*0 - g2(K1,K2,x2,y2)*1 - 4*g2(K1,K2,xmid,ymid)*0.5 ;
    if y1==0
        G1(1)=1*G1(1);G2(1)=1*G2(1);  % 与Gamma1接触的
    elseif y2==1
        G1(2)=0*G1(2);G2(2)=0*G2(2);
    end
end
