function [kk,ff] = constraints(kk,ff,bcdof,bcval)

% apply constraint to matrix [K]{x} = {F}

%k stiffness matrix before applying constrain
% f force vector befor apply constrain

% n=length(BC_DoF);
% sdof=size(k);
% 
% for i=1:n
%     c=BC_DoF(i);
%     for j=1:sdof
%         k(c,j)=0;
%     end
%     
%     k(c,c)=1;
%     f(c)=BC_Value(i);
n=length(bcdof);
 sdof=size(kk);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
    end

    kk(c,c)=1;
    ff(c)=bcval(i);
 end
end