function [tx1 result flag order tfval] = GetPointsReg(p1,p2)
    l1 = size(p1,1);
    l2 = size(p2,1);
%     if(l1>=l2)
%         x1_ = p2;
%         x2_ = p1;
%     else
%         x1_ = p1;
%         x2_ = p2;
%     end
    x1_ = p1;
    x2_ = p2;
    if(min(l1,l2)<5)
        warning('need more points');
    end
    
    initpar = [1 0 0 1 0 0];
    %tol = 0.0000;
    options = optimset('Display ','off','MaxFunEvals',3000,'MaxIter',2000,'TolFun',1e-9,'TolX',1e-7);
    %[rzlt] = lsqnonlin(@EvalPointSet,initpar,[],[],options,x1_,x2_);
    f=@(par)EvalPointSet(par,x1_,x2_);
%     tic
    [result tfval exitflag output] = fminsearch(f,initpar,options);
    if(tfval>1e-1)%结果不好
        flag=0;
    else
        flag=1;
    end
%     toc
%     rx1=zeros(0,2);
%     rx2=zeros(0,2);
    tx1 = x1_*[result(1) result(2)
               result(3) result(4)];
    tx1(:,1)=tx1(:,1)+result(5);
    tx1(:,2)=tx1(:,2)+result(6);

    [order] = FindPointOrder(tx1,x2_,0.02);
%     rx1 = x1_(order(:,1),:);
%     rx2 = x2_(order(:,2),:);
end