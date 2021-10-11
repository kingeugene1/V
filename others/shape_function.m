function [N,dNdx, dNdy] = shape_function (xi,eta)
    N=zeros(4,1);
    N(1) =  1/16*(xi-1).^2.*(xi+2)*(eta-1).^2.*(eta+2); 
  N(2) =  1/16*(xi-1).^2.*(xi+1)*(eta-1).^2.*(eta+2); 
  N(3) = -1/16*(xi+1).^2.*(xi-2)*(eta-1).^2.*(eta+1); 
  N(4) =  1/16*(xi+1).^2.*(xi-1)*(eta-1).^2.*(eta+1);
    
    dNdx=zeros(4,1);
    dNdx(1,1) = 3/16*(xi^2-1)*(eta-1).^2.*(eta+2);
    dNdx(2,1) = 1/16*(3*xi^2-2*xi-1)*(eta-1).^2.*(eta+2);
    dNdx(3,1) = -3/16*(xi^2-1)*(eta-1).^2.*(eta+1);
    dNdx(4,1) = 1/16*(3*xi^2+2*xi-1)*(eta-1).^2.*(eta+1);
    
        dNdy=zeros(4,1);
        dNdy(1,1) = 3/16*(xi-1).^2.*(xi+2)*(eta^2-1);
        dNdy(2,1) = 3/16*(xi-1).^2.*(xi+1)*(eta^2-1);
        dNdy(3,1) = -1/16*(xi+1).^2.*(xi-2)*(3*eta^2-2*eta-1);
        dNdy(4,1) = 1/16*(xi+1).^2.*(xi-1)*(3*eta^2-2*eta-1);
end