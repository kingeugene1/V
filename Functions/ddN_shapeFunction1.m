function [ddN1,ddN2,ddN12]= ddN_shapeFunction1()
 ddN1  = zeros(4,1);
  
      ddN1(1) =0;
      ddN1(2) =0;
      ddN1(3) =0;
      ddN1(4) = 0;
      
 ddN2  = zeros(4,1);
      
      ddN2(1) = 0;
      ddN2(2) = 0;
      ddN2(3) =0;
      ddN2(4) = 0;
      
      ddN12  = zeros(4,1);
  
      ddN12(1) =  1/4;
      ddN12(2) = -1/4;
      ddN12(3) =  1/4;
      ddN12(4) = -1/4;
end