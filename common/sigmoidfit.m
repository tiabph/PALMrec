function [Tnew,fitX,R]=sigmoidfit(T,X,flag)

% T=a(:,1);
% X=a(:,2);
% flag=1;
len=length(X);
base=min(X);
maxx=max(X)-min(X);
% if flag==1
%    maxx=100;
% else
%    maxx=-100;
% end
%  xhalf=T(round(len/2),1);
 xhalf=70;
 rate=abs(maxx)/len*2;
 
 fp=[base,maxx,xhalf,rate,-maxx,2*xhalf,rate];
%  fp=[base,maxx,xhalf,rate];
 options = optimset('Display','off','MaxIter',1000,'TolFun',1e-9,'LargeScale','off');
 b = lsqcurvefit(@sigmoid2,fp,T,X,[],[],options);
%  options=statset('Display','off','MaxIter',1000,'TolFun',1e-9);
%  b = nlinfit(T,X,@sigmoid,fp,options);
 
 fpp=b;
 fitX=fpp(1)+fpp(2)./(1+exp((fpp(3)-T)./fpp(4)))+fpp(5)./(1+exp((fpp(6)-T)./fpp(7)));
% fitX=fpp(1)+fpp(2)./(1+exp((fpp(3)-T)./fpp(4)));
Tnew=1;
R=abs(fpp(2)/fpp(1));
%  if fpp(2)<0   
%       ninetieth=fpp(1)+0.1*fpp(2);
%  else
%       ninetieth=fpp(1)+0.9*fpp(2);
%  end
% 
%  for i=1:len-1
%      if (fitX(i,1)-ninetieth)*(fitX(i+1,1)-ninetieth)<0
%          ninetiethI=i;
%          break;
%      else
%          ninetiethI=0;
%      end
%  end
%  
%  if ninetiethI>0
%  Tnew=T(ninetiethI,1);
%  else
%      Tnew=NaN;
%  end
%  
%  if flag==3  &&  fpp(2)>0
%      Tnew=NaN;
%  end
%  if flag==1  &&  fpp(2)<0
%      Tnew=NaN;
%  end
%  R=abs(fpp(2)/fpp(1));
%  Tnew
%  plot(T,X,T,fitX);
 
 