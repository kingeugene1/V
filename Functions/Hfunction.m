function [H]= Hfunction(xi,eta)
H= zeros(12,1);

N1= Nfunction(xi);
N2= Nfunction(eta);
  H(1) = N1(1)*N2(1);        
  H(2) = N1(2)*N2(1);           
  H(3) = N1(1)*N2(2);           
  
  
  %(X,Y)=(1,-1)
  
  H(4) = N1(3)*N2(1);        
  H(5) = N1(4)*N2(1);           
  H(6) = N1(3)*N2(2);          
  
  
  %(X,Y)=(1,1)
  
  H(7) = N1(3)*N2(3);        
  H(8) = N1(4)*N2(3);           
  H(9) = N1(3)*N2(4);           
 
  %(X,Y)=(-1,1)
  
  H(10) = N1(1)*N2(3);        
  H(11) = N1(2)*N2(3);          
  H(12) = N1(1)*N2(4);           
end