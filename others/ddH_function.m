function [N]= ddH_function(xi,eta,a)

  N= zeros(12,1);
  P1=Nfunction(xi);
  P2=Nfunction(eta);
  dP1=dNfunction(xi);
  dP2=dNfunction(eta);
  dN1=ddNfunction (xi);
  dN2=ddNfunction (eta);
  
  if a==0 % X-direction second differential --- ddN/dX^2
  
      N(1) = dN1(1)*P2(1);        %  w(0,0)
      N(2) = dN1(2)*P2(1);           %  0x(0,0)
      N(3) = dN1(1)*P2(2);                   %  0y(0,0)
      
  
      N(4) = dN1(3)*P2(1);        %  w(L,0)
      N(5) = dN1(4)*P2(1);           %  0x(L,0)
      N(6) = dN1(3)*P2(2);                   %  0y(L,0)
     
      
      N(7) = dN1(3)*P2(3);        %  w(L,L)
      N(8) = dN1(4)*P2(3);           %  0x(L,L)
      N(9) = dN1(3)*P2(4);                   %  0y(L,L)
     
      
      N(10) = dN1(1)*P2(3);        %  w(0,L)
      N(11) = dN1(2)*P2(3);           %  0x(0,L)
      N(12) = dN1(1)*P2(4);                   %  0y(0,L)
      
      
      
  elseif a==1  % Y-direction second differential --- ddN/dY^2
      
       N(1) = P1(1)*dN2(1);        %  w(0,0)
       N(2) = P1(2)*dN2(1);                    %  0x(0,0)
       N(3) = P1(1)*dN2(2);           %  0y(0,0)
       
       N(4) = P1(3)*dN2(1);        %  w(L,0)
       N(5) = P1(4)*dN2(1);                    %  0x(L,0)
       N(6) = P1(3)*dN2(2);           %  0y(L,0)
  
       N(7) = P1(3)*dN2(3);        %  w(L,L)
       N(8) = P1(4)*dN2(3);                   %  0x(L,L)
       N(9) = P1(3)*dN2(4);           %  0y(L,L)
       
       N(10) = P1(1)*dN2(3);        %  w(0,L)
       N(11) = P1(2)*dN2(3);                   %  0x(0,L)
       N(12) = P1(1)*dN2(4);           %  0y(0,L)
       
       
       
  else  % X,Y-direction differential --- ddN/dXdY
      
       N(1) = dP1(1)*dP2(1);        %  w(0,0)
       N(2) = dP1(2)*dP2(1);           %  0x(0,0)
       N(3) = dP1(1)*dP2(2);           %  0y(0,0)
  
       N(4) = dP1(3)*dP2(1);        %  w(L,0)
       N(5) = dP1(4)*dP2(1);           %  0x(L,0)
       N(6) = dP1(3)*dP2(2);           %  0y(L,0)
  
       N(7) = dP1(3)*dP2(3);        %  w(L,L)
       N(8) = dP1(4)*dP2(3);           %  0x(L,L)
       N(9) = dP1(3)*dP2(4);           %  0y(L,L)
      
       N(10) = dP1(1)*dP2(3);        %  w(0,L)
       N(11) = dP1(2)*dP2(3);           %  0x(0,L)
       N(12) = dP1(1)*dP2(4);           %  0y(0,L)
       
  end
end