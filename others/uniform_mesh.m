function [NL,EL]=uniform_mesh(d1,d2,p,m)
% This function returns the Node coordinates and Elements nodes as matrices
% Written by Eugene Mudahemuka 2021, 09/06
% NL node list with coordinates
% EL element list with nodes
 pd=2;                                          % dimension of the problem

 NoN= (p+1)*(m+1);                   % Total number of nodes
NL= zeros(NoN,pd);                     % Matrix holding node list
 dx=d1/p;                                     % increment along x direction
 dy=d2/m;                                    % increment along x direction
 
%  N=zeros(121,2);
jj=repmat((1:1:p+1)',p+1,1);
ii=repelem((1:1:m+1)',m+1,1);
NL(:,1)=(jj-1)*dx;
NL(:,2)=(ii-1)*dy;
% Elements
 
EL=[];
for yy=0:m-1
     for xx=0:p-1
        
            X1=     xx * dx;
            X2= (xx+1) * dx;
            
            Y1=     yy * dy;
            Y2= (yy+1) * dy;
           
            idx = zeros(4,1);
            
            idx(1)  = find(sum(NL(:,1:2) == [X1,Y1],2)==2);
            idx(2)  = find(sum(NL(:,1:2) == [X2,Y1],2)==2);
            idx(3)  = find(sum(NL(:,1:2) == [X2,Y2],2)==2);
            idx(4)  = find(sum(NL(:,1:2) == [X1,Y2],2)==2);
           
            
            EL  = [EL; idx'];
        
    end
end
end