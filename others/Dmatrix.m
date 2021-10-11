function [D]=Dmatrix(E,N,t,nu)
for i=1:N
 D0= t^3/(12*(1-nu^2));
 
 D= E(i)*D0*[1 nu 0;
              nu 1 0;
              0 0 0.5*(1-nu)];
end
end