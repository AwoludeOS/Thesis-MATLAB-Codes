% Code to generate G(n,p).
% Input: number of nodes n and edge probability p.
% Output: 
% (i) Matrix with n rows, row i contains neighbours of node i.
% (ii) Vectorized form of this matrix, containing all neighbours of 1, then
%      all  neighbours of 2 and so on. Length: sum of degrees.
% (iii) Vector of length n giving degrees of all nodes.

function [NMat,NVec,Degr]=Erdos_Renyi(n,p)


    
    m=nchoosek(n,2);%maximal possible number of edges
    
    X = binornd(m,p);%(varying) number of edges in the graph 
    
    
    NMat=zeros(n,ceil(n*p+4*sqrt(n*p)));
    Degr=zeros(n,1);
    
    for i = 1:X
        success = 0;
        while success == 0
            x=ceil(rand*n);
            y=ceil(rand*n);
            if x ~= y && sum(NMat(x,1:Degr(x)) == y) == 0
                success =1;
                Degr(x)= Degr(x)+1;
                Degr(y)= Degr(y)+1;
                NMat(x,Degr(x)) = y;
                NMat(y,Degr(y)) = x;
            end
        end
    end
    
    NVec = zeros(2*X,2);
    index = 1;
    
    for i = 1:n
        for d = 1:Degr(i)
            NVec(index,1) = i;
            NVec(index,2) = NMat(i,d);
            index = index+1;
        end
    end
end
        