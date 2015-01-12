function bfactors=bleachFactors(A,n)
minI=min(min(A(:,:,1)));
maxI=max(max(A(:,:,1)));
m=[];
for i=1:n
   a=A(:,:,i);
   a=a(:);
   a=(a-minI)/(maxI-minI);
   m(i)=mean(a);
end
maxM(1:n)=max(m);
bfactors=maxM./m;