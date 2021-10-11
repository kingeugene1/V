function M=M_vector(NPE,N)
%m= N.*N;
for i=1:3*NPE
    for j= 1:3*NPE
        M(i,j)= N(i)*N(j);
    end
end
end