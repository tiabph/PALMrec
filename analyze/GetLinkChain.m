function result = GetLinkChain(linkData_chain, linkData_sp, chainlen)
    if nargin<3
        chainlen=10000;
    end
    result = zeros(chainlen,1);
    currPos = linkData_sp;
    chainCnt = 0;
    while(1)
        chainCnt = chainCnt+1;
        result(chainCnt) = currPos;
        currPos = linkData_chain(currPos);
        if currPos<0
            break;
        else
            currPos = currPos +1;
        end
    end
    
    result = result(1:chainCnt);
end