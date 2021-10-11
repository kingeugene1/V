function [K,M]= Assemble(K,M,k,idx,me)

edof= length(idx);
% K = zeros (DoF,DoF);
% M=zeros (DoF,DoF);
for i=1:edof
%     for c=1:n
    ii=idx(i);
    for j=1:edof
        jj=idx(j);
        K(ii,jj)=K(ii,jj)+k(i,j);
        M(ii,jj)=M(ii,jj)+me(i,j);
    end
end
end