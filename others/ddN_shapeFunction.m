function [ddN1,ddN2,ddN12]= ddN_shapeFunction(xi,eta)
 ddN1  = zeros(4,1);
  
      ddN1(1) =3/8*(xi)*(eta-1).^2.*(eta+2) ;
      ddN1(2) = 1/8*(3*xi-1)*(eta-1).^2.*(eta+2);
      ddN1(3) = -3/8*(xi)*(eta-1).^2.*(eta+1);
      ddN1(4) = 1/8*(3*xi+1)*(eta-1).^2.*(eta+1);
      
 ddN2  = zeros(4,1);
      
      ddN2(1) = 3/8*(xi-1).^2.*(xi+2)*(eta);
      ddN2(2) = 3/8*(xi-1).^2.*(xi+1)*(eta);
      ddN2(3) = -1/8*(xi+1).^2.*(xi-2)*(3*eta-1);
      ddN2(4) = 1/8*(xi+1).^2.*(xi-1)*(3*eta-1);
      
      ddN12  = zeros(4,1);
  
      ddN12(1) =  3/4*(xi)*(eta);
      ddN12(2) = 1/4*(3*xi-1)*(eta);
      ddN12(3) =  -3/4*(xi)*(3*eta-1);
      ddN12(4) = 1/4*(3*xi+1)*(3*eta-1);
end