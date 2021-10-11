function   [ke]= Element_stiffness(NPE,dNdx, dNdy)
 for i= 1:NPE
   index1=(i-1)*3+1;
   index2= (i-1)*3+2;
   index3= (i-1)*3+3;
   
   ke(2,index2)= dNdy(i);
   ke(3,index2)= dNdx(i);
   ke(1,index3)=-dNdx(i);
   ke(3,index3)=-dNdy(i);
   ke(3,index1)= 0;
end
   
end