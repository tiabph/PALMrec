function f=sigmoid2(fp,T)

f=fp(1)+fp(2)./(1+exp((fp(3)-T)./fp(4)))+fp(5)./(1+exp((fp(6)-T)./fp(7)));

