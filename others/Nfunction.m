function [N]=Nfunction(X)
  N= zeros(4,1);
  N(1) =  1/4*(X-1).^2.*(X+2); 
  N(2) =  1/4*(X-1).^2.*(X+1); 
  N(3) = -1/4*(X+1).^2.*(X-2); 
  N(4) =  1/4*(X+1).^2.*(X-1); 
end