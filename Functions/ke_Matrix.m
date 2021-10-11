function [A,G,B]= ke_Matrix(xi,eta,x,y)
[~,dNdxi, dNdeta] = shape_function1(xi,eta);

[ddN1,ddN2,ddN12]= ddN_shapeFunction1();

% derivative of global cordinate wrt local cordinate(xi,eta)
dxdX = x * dNdxi; 
dxdY = x* dNdeta; 
dydX = y * dNdxi;
dydY = y * dNdeta; 

% Second derivative of global cordinate wrt to local cordinates (xi,eta)

ddxdX2 = x * ddN1;
ddxdXdY = x * ddN12;
ddxdY2 = x * ddN2; 
ddydX2 = y* ddN1; 
ddydXdY = y * ddN12;
ddydY2 = y * ddN2; 

%j=dxdX*dydY-dydX*dxdY; % jacobian Matrix
% call for Herimite shape functions
% H= Hfunction(xi,eta);


% First differential of the hermite function H 
% in the local coordinate system (X,Y)
% (of the hermite element)
% [N,dHdxi,dHdeta]= dH_function(xi,eta,a)
dHdX = dH_function(xi,eta,0);
dHdY = dH_function(xi,eta,1);

% Second differential of the herimite function H 
% in the local coordinate system (X,Y)
% (of the hermite element)

ddHdX2  = ddH_function(xi,eta,0);
ddHdY2  = ddH_function(xi,eta,1);
ddHdXdY = ddH_function(xi,eta,2);


j=[ 1    0    0;
    0 dxdX dxdY;
    0 dydX dydY;];

JJ(1:3,1:3) = j;
JJ(4:6,4:6) = j;
JJ(7:9,7:9) = j;
JJ(10:12,10:12) = j;
dNdX = JJ*dHdX;
dNdY = JJ*dHdY;

ddNdX2  = JJ*ddHdX2;
ddNdY2  = JJ*ddHdY2;
ddNdXdY = JJ*ddHdXdY;

G = [dNdX.'; dNdY.'; ddNdX2.'; ddNdY2.'; ddNdXdY.'];


C = [ dxdX    dydX    0           0           0;
      dxdY    dydY    0           0           0;
      ddxdX2  ddydX2  (dxdX^2)    (dydX^2)    (2*dxdX*dydX);
      ddxdY2  ddydY2  (dxdY^2)    (dydY^2)    (2*dxdY*dydY);
      ddxdXdY ddydXdY (dxdX*dxdY) (dydX*dydY) (dxdX*dydY+dxdY*dydX)];
 AA  = inv(C);  
A = [AA(3,:);
       AA(4,:);
       2*AA(5,:)];
   B = [AA(3,:);
       AA(4,:);
       AA(5,:)];

end