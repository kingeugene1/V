function [N]= dNfunction(X)
 
  N  = zeros(4,1);
  
  N(1) =  3/4*(X^2-1);
  N(2) =  1/4*(3*X^2-2*X-1);
  N(3) = -3/4*(X^2-1);
  N(4) =  1/4*(3*X^2+2*X-1);
  
end