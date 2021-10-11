function [N] =ddNfunction (X)
N  = zeros(4,1);
  
  N(1) =  3/2*(X);
  N(2) =  1/2*(3*X-1);
  N(3) = -3/2*(X);
  N(4) =  1/2*(3*X+1);
end