a=length(V);
V(:,4)=0;
n=0;
for i=1:a-1
    i
    for j=i+1:a
        if sqrt((V(i,2)-V(j,2))^2+(V(i,1)-V(j,1))^2)<0.3 %&& i~=j
            if V(i,4)==0 && V(j,4)==0
                n=n+1;
                V(i,4)=n;
                V(j,4)=n;                
            end
            if V(i,4)>0 && V(j,4)==0
                V(j,4)=V(i,4);
            end
            if V(j,4)>0 && V(i,4)==0
                V(i,4)=V(j,4);
            end
            if V(j,4)>0 && V(i,4)>0
                [IX II]=find(V(:,4)==V(j,4));
                V(IX,4)=V(i,4);
            end
        end
    end
end
N=zeros(n,1);
T=[];
t=[];
m=0;
for i=1:n
    N(i,1)=sum(V(:,4)==i);
    t=V(V(:,4)==i,3);
    t=sort(t);
    for j=1:length(t)-1
        m=m+1;
        T(m,1)=t(j+1)-t(j);
    end
end
T=T*0.1;
N=N(N>0);
bin=(2:max(N))';
out=(hist(N,bin))';
R(1,1)=sum(V(:,4)==0);
R(2:10,1)=out(1:9);
R=R./(sum(out)+R(1,1))*100;
