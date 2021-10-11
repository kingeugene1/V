function [dNdx, dNdy]= dN_global(NPE,dNdX,dNdY,inv_J)
for i= 1:NPE
    dNdx(i)= inv_J(1,1)*dNdX(i)+inv_J(1,2)*dNdY(i);
    dNdy(i)= inv_J(2,1)*dNdX(i)+inv_J(2,2)*dNdY(i);
end
end