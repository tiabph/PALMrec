function f=sigmoid(fp,T)

f=fp(1)+fp(2)./(1+exp((fp(3)-T)./fp(4)));

