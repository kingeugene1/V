function [det_J, inv_J]= Jacobian(dNdX,dNdY,xx,yy)
 
Jacobian= zeros(2,2);


    Jacobian(1,1)= Jacobian(1,1)+ xx*dNdX;
    Jacobian(1,2)= Jacobian(1,2)+ xx*dNdY;
    Jacobian(2,1)= Jacobian(2,1)+ yy*dNdX;
    Jacobian(2,2)= Jacobian(2,2)+ yy*dNdY;
    

det_J= det(Jacobian); % Determinant of Jacobian Matrix
inv_J= inv(Jacobian); % inverse of jacobian Matrix


end