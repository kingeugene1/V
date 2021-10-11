function [N,dNdx, dNdy] = shape_function1(xi,eta)
    N=zeros(4,1);
    N(1,1) = 1/4*(1-xi)*(1-eta);
    N(2,1) = 1/4*(1+xi)*(1-eta);
    N(3,1) = 1/4*(1+xi)*(1+eta);
    N(4,1) = 1/4*(1-xi)*(1+eta);
    
    dNdx=zeros(4,1);
    dNdx(1,1) = -1/4*(1-eta);
    dNdx(2,1) = 1/4*(1-eta);
    dNdx(3,1) = 1/4*(1+eta);
    dNdx(4,1) = -1/4*(1+eta);
    
        dNdy=zeros(4,1);
        dNdy(1,1) = -1/4*(1-xi);
        dNdy(2,1) = -1/4*(1+xi);
        dNdy(3,1) = 1/4*(1+xi);
        dNdy(4,1) = 1/4*(1-xi);
end