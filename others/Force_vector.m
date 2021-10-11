function F=Force_vector(NPE,N,P)
    f= N*P;
    for i=1:NPE
       index1= (i-1)*3 +1;
       index2= (i-1)*3 +2;
       index3= (i-1)*3 +3;
       F(index1,1)= f(1);
       F(index2,1)= 0;
       F(index3,1)= 0;
%F(i)=N(i)*P;
    end
end