function [index]=elementdof(node,nnel,ndof)
     
 k=0;
   for i=1:nnel
     start = (node(i)-1)*ndof;
       for j=1:ndof
         k=k+1;
         index(k)=start+j;
       end
   end