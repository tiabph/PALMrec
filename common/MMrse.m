function [beta,sigma,SE1,SE]=MMrse(data,bdp)
% 
% Computes MM-estimator of regression with Robust standard errors, as in
%       Croux, C., Dhaene, G., Hoorelbeke, D. (2003), 
%       ``Robust Standard Errors for Robust Estimators'' 
% Starting values of the MM-estimator is fast-S-estimator 
%       Salibian-Barrera, M. and Yohai, V. (2005),
%       ``A fast algorithm for S-regression estimates.
%       (Translated to Matlab by Kristel Joossens, K.U.Leuven, Belgium)
%
% INPUT
%      data : [Y,X] matrix of size n x (K+1), containing the response 
%             vector in the first column and the design matrix in the 
%             others. Add a column of ones if an intercept is in the 
%             regression model.
%      bdp  : breakdown point of the MM-estimator (for example 0.25 or 0.50)
%             Efficiency of MM-estimator is equal to 95%
%
% OUTPUT
%      beta  : estimated regression coefficients
%      sigma : estimated regression scale
%      SE1   : standard errors around beta, robust for heteroskedasticity
%      SE    : standard errors aroumd beta, robust for heteroskedasticity 
%              and autocorrelation
%
% Direct calls: auxrobcovmm_ng0, mmregres, fasts
% Indirect calls: dpsibi, fw, gint, lossS, oursolve, psibi, ress, rhobi,
%                 scale1, Tbsb, Tbsc

warning off;
[n,p]=size(data);
Y=data(:,1);
X=data(:,2:p);
resfs=fastsreg(X,Y,bdp);
bs=resfs.beta;
ss=resfs.scale;
bm=mmregres(X,Y,bs,ss,0);
em=(Y-X*bm)/ss;
es=(Y-X*bs)/ss;    
[cov3,cov2s,cov1s,cov1,cov]=auxrobcovmm_ng0(es,em,ss,bdp,X);
SE1=diag(sqrt(cov1));
SE=diag(sqrt(cov));
beta=bm;
sigma=ss;

%--------------------------------------------------------------------------
% Direct calls
%-------------------------------------------------------------------------- 
function [cov3,cov2s,cov1s,cov1,cov]=auxrobcovmm_ng0(es,em,s,bdp,X);

%INPUT
%     es : S-residuals 
%     em : MM-residuals
%     s  : scale estimate
%     c  : constant TB
%
% Direct calls: rhobi, psipi, dpsibi, Tbsc, Tbsb
% Indirect calls: gint

XX=X'*X;
n=max(size(X));
k=min(size(X)); % # parameters is k+1  
c0=Tbsc(bdp,1); %initial S-estimator  
b0=Tbsb(c0,1)*c0/6; 
c1=4.685; %MM-estimator

rho0=rhobi(es,c0);  
rho0m=rhobi(em,c0);
drho0=psibi(es,c0);  
drho0m=psibi(em,c0);
ddrho0=dpsibi(es,c0);
psi=psibi(em,c1);  
dpsi=dpsibi(em,c1);  

Epsi2=(1/n)*psi'*psi; 
EXX=(1/n)*XX;  
Edpsi=(1/n)*sum(dpsi);  
cov2s=(1/n)*s^2*Epsi2*inv(EXX)/Edpsi^2;
cov3=(1/n)*s^2*1.0526*inv(EXX);

%A-matrix (E dg/dtheta) 
h1=dpsi*ones(1,k);  
h1=h1.*X;  
EdpsiXX=(1/n)*X'*h1;  
h1=ddrho0*ones(1,k); 
h1=h1.*X;
Eddrho0XX=(1/n)*X'*h1; 
EdpsiXem=(1/n)*(dpsi.*em)'*X;  
Eddrho0Xes=(1/n)*(ddrho0.*es)'*X;  
Edrho0es=(1/n)*drho0'*es; 
A=(-1/s)*[EdpsiXX    EdpsiXem'   zeros(k,k); 
          zeros(1,k) Edrho0es    zeros(1,k);    
          zeros(k,k) Eddrho0Xes' Eddrho0XX]; 
  
%B-matrix (E g_t g'_t)
h1=psi*ones(1,k); 
h0=drho0*ones(1,k); 
G=[h1.*X rho0-b0 h0.*X]; % elements g_t 
B=(1/n)*G'*G; 

iA=inv(A);
cov1=(1/n)*iA*B*iA'; 
cov1=cov1(1:k,1:k);

h3=inv(EdpsiXX);
cov1s=s^2*(1/n)*h3*B(1:k,1:k)*h3;

q=floor(4*(n/100)^(2/9)); % truncation parameter Bartlett kernel 
S=B; 
for i=1:q 
    Si=(1/n)*G(i+1:n,:)'*G(1:n-i,:);
    S=S+(1-i/(q+1))*(Si+Si'); 
end
cov=(1/n)*iA*S*iA';
cov=cov(1:k,1:k);

% -------------------------------------------------------------------------

function result = fastsreg(x, y, bdp, N) 

% This is the Matlab program implemented by Kristel Joossens (K.U.Leuven, 
% Belgium), which is a translation from the R program of 
% Salibian-Barrera, M. and Yohai, V.J. (2005)
% Reference paper: "A fast algorithm for S-regression estimates" (2005)
% The R software can also be found on 
%           http://hajek.stat.ubc.ca/~matias/soft.html
%
% INPUT
%      x   : data matrix         (n x p)
%            add a column of ones if you which an intercept
%      y   : response matrix     (n x 1)
%      bdp : breakdownpoint      (e.g. .5)
%      N   : cant de sub-samples (if not given, default = 20)
% Default arguments, redefineable in the beginning of the program
%      k : number of refining iterations in each subsample (default = 2)
%          (k = 0 means "raw-subsampling")
%      bestr : number of "best betas" to remember from the subsamples. 
%              These will be later iterated until convergence (default=5)
%
% OUTPUT is a list with components
%      beta : robust estimate of the regression coefficients (p x 1)
%      scale: values of the objective function               (1 x 1)
%
% Direct calls: lossS, oursolve, ress, scale1, Tbsb, Tbsc
% Indirect calls: fw, gint, rhobi

if nargin == 3
   N= 20; 
end
k=2; bestr=5; seed=0; 

[n,p] =size(x);
c = Tbsc(bdp,1);
kp = (c/6) * Tbsb(c,1);

bestbetas = zeros(bestr,p);
bestscales = 1e20 * ones(bestr,1);
sworst = 1e20;
nref = 1;
    
for i=1:N
%    i    %<- print this index to see how fast the algorithm works     
    % get a subsample
    singular =1; itertest=1;
    while (singular==1 && itertest<100)
         [dummysort,index]=sort(rand(n,1));
         index=index(1:p);
         xs = x(index,:);
        ys = y(index);
        beta = oursolve(xs,ys);
        singular = any(isnan(beta));  
        itertest=itertest+1;
    end
    if itertest==100
        error('too many degenerate subsamples')
    end 
    if k>0  
       % do the refining
       tmp = ress(x,y,beta,k,0,kp,c);
       betarw = tmp.betarw;
       scalerw = tmp.scalerw;
       resrw = y - x * betarw;
   else
       % k = 0 means "no refining"
       betarw = beta;
       resrw = y - x * betarw;
       scalerw = median(abs(resrw))/.6745;
   end
   if i > 1  
      % if this isn't the first iteration....
      scaletest = lossS(resrw,sworst,c);
      if scaletest < kp
          sbest = scale1(resrw,kp,c,scalerw);
          [yss,yi]=sort(bestscales);
          ind=yi(bestr);
          bestscales(ind) = sbest;
          bestbetas(ind,:) = betarw';
          sworst = max(bestscales);
      end
    else 
        % if this is the first iteration, then this is the best beta...
        bestscales(bestr) = scale1(resrw,kp,c,scalerw);
        bestbetas(bestr,:) = betarw';
    end
end

% do the complete refining step until convergence (conv=1) starting
% from the best subsampling candidate (possibly refined)
superbestscale = 1e20;
% magic number alert
for i=bestr:-1:1
    tmp = ress(x,y,bestbetas(i,:)',0,1,kp,c,bestscales(i));
    if tmp.scalerw < superbestscale
      superbestscale = tmp.scalerw;
      superbestbeta = tmp.betarw;
    end
end

result.beta=superbestbeta;
result.scale=superbestscale;

%--------------------------------------------------------------------------

function res=mmregres(X,Y,b0,s,bdp)
% INPUT
%      X   : data matrix         (n x p)
%      Y   : response matrix     (n x 1)
%      b0,s: initial S-estimate with high breakdownpoint (eg .5)
%      bdp : bdp for MM-estimation, choose bdp==zero for .95 efficiency at
%      normal  %
% Direct calls: Tbsc, psibi
% Indirect calls: Tbsb, gint
k=min(size(X));

if bdp == 0
    c=4.685;    
else
    c=Tbsc(bdp,1);
end
maxit=100;tol=10^(-10);eps=10^(-200);
iter=0;crit=1000;
b1=b0;
while (iter <= maxit) && (crit > tol)    
    r1=(Y-X*b1)/s;
    tmp = find(abs(r1) <= eps);  
    [n1 n2] = size(tmp);
    if n1 ~= 0  
        r1(tmp) = eps;
    end;   
    w=psibi(r1,c)./r1;
    W=w*ones(1,k);
    XW=X'.*W';    
    b2=inv(XW*X)*XW*Y;
    d=b2-b1;
    crit=max(abs(d));
    iter=iter+1;
    b1=b2;
end
res=b2;

%--------------------------------------------------------------------------
% Indirect calls
%-------------------------------------------------------------------------- 

function z=dpsibi(x,c)
% computes derivative of tukey's biweight psi function with constant c for 
% all values in the vector x.
%
% Direct calls: none
z = (abs(x) < c) .* (1 - x.^2 .*(6/c^2 - 5*x.^2/c^4));

%--------------------------------------------------------------------------

function tmp = fw(u,c)
% weight function = psi(u)/u
%
% Direct calls: none
tmp = (1 - (u/c).^2).^2;
tmp = tmp .* (c^2/6);
tmp( abs(u/c) > 1 )= 0;

%--------------------------------------------------------------------------  

function res=gint(k,c,p)

% this procedures computes the integral from zero to c of
% the function r^k g(r^2), where g(||x||^2) is the density function
% of a p-dimensional standardnormal distribution 
% 
% Direct calls: none
e=(k-p-1)/2;
numerator=(2^e)*gamcdf((c^2)/2,(k+1)/2)*gamma((k+1)/2);
res=(numerator/(pi^(p/2)));

%--------------------------------------------------------------------------

function res = lossS(u,s,c)
% the objective function, we solve loss.S(u,s,c)=b for "s"
%
% Direct calls: rhobi
res = mean(rhobi(u/s,c));

%--------------------------------------------------------------------------  

function res = oursolve(a,b) 
% Direct calls: none
[Q,R]=qr(a);
p = size(Q,2);
if rank(a) < p
    res = NaN;
else 
    res = inv(a) * b;
end

%--------------------------------------------------------------------------  

function z=psibi(x,c)
% psi function for biweight.  Use for robust estimation.
% c : tuning parameter
%
% Direct calls: none
z = (abs(x) < c) .* x .* ( 1 - (x./c).^2 ).^2 ;

%--------------------------------------------------------------------------

function v = rhobi(u,c) 
% rho function for biweight.  Used for robust scale estimation.
% c : tuning parameter
%
% Direct calls: none
w = (abs(u)<=c);
v = (u.^2/(2).*(1-(u.^2/(c^2))+(u.^4/(3*c^4)))).*w +(1-w)*(c^2/6);

%--------------------------------------------------------------------------

function result = ress(x,y,initialbeta,k,conv,kp,c,initialscale) 
% does "k" IRWLS refining steps from "initial.beta"
%
% if "initial.scale" is present, it's used, o/w the MAD is used
% k = number of refining steps
% conv = 0 means "do k steps and don't check for convergence"
% conv = 1 means "stop when convergence is detected, or the maximum number
%                 of iterations is achieved"
% kp and c = tuning constants of the equation
% 
% Direct calls: fw, oursolve, rhobi
[n,p]=size(x);    
res = y - x * initialbeta;
if (nargin < 8)
    scale = median(abs(res))/.6745;
    initialscale = scale;
else
    scale = initialscale;
end

if (conv == 1) 
    k = 50;
end
% if conv == 1 then set the max no. of iterations to 50 magic number alert!

beta = initialbeta;
scale = initialscale;
lowerbound = median(abs(res))/c;
for i=1:k
    % do one step of the iterations to solve for the scale
    scalesuperold = scale;
    %lower.bound <- median(abs(rdis))/1.56
    scale = sqrt( scale^2 * mean( rhobi(res/scale,c) ) / kp );
    % now do one step of IRWLS with the "improved scale"
    weights = fw(res/scale,c);
    sqweights = weights.^(1/2);
    sqW = sqweights * ones(1,p);
    xw = x .* sqW;
    yw = y .* sqweights;
    beta1 = oursolve(xw'*xw,xw'*yw);
    if (any(isnan(beta1))) 
        beta1 = initialbeta;
        scale = initial.scale;
        break
    end
    
    if  (conv==1)
        % check for convergence
        if ( norm( beta - beta1 ) / norm(beta) < 1e-20 ) 
            break
            % magic number alert!!!
        end 
    end
    res = y - x * beta1;
    beta = beta1;
end
res = y - x * beta;
% get the residuals from the last beta
result.betarw = beta1;
result.scalerw = scale;

%--------------------------------------------------------------------------  

function sc = scale1(u, kp, c, initialsc) 
% Direct calls: rhobi
if nargin<3
    initialsc = median(abs(u))/.6745;
end
% find the scale, full iterations
max.it = 200;
% magic number alert
% sc = median(abs(u))/.6745;
sc = initialsc;
i = 0; 
eps = 1e-20;
% magic number alert
err = 1;
while  (( i < max.it ) && (err > eps))
    sc2 = sqrt( sc^2 * mean(rhobi(u/sc,c)) / kp);
    err =abs(sc2/sc - 1);
    sc = sc2;
    i=i+1;
end

%-------------------------------------------------------------------------- 

function res=Tbsb(c,p)
% Direct calls: gint
y1=gint(p+1,c,p)/2-gint(p+3,c,p)/(2*c^2)+gint(p+5,c,p)/(6*c^4);
y2=(6/c)*2*(pi^(p/2))/gamma(p/2);
y3=c*(1-chi2cdf(c^2,p));
res=(y1*y2+y3);

%-------------------------------------------------------------------------- 

function res = Tbsc(alpha,p)
% constant for Tukey Biweight S 
%
% Direct cals: Tbsb
talpha = sqrt(chi2inv(1-alpha,p));
maxit = 1000; 
eps = 10^(-8);
diff = 10^6;
ctest = talpha;
iter = 1;
while ((diff>eps) && iter<maxit)
    cold = ctest;
    ctest = Tbsb(cold,p)/alpha;
    diff = abs(cold-ctest);
    iter = iter+1;
end
res = (ctest);