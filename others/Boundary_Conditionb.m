function BC_DoF= Boundary_Conditionb(BC_Type,NL)

L1= find(NL(:,2)==min(NL(:,2))); % y=0
L2= find(NL(:,1)==max(NL(:,1))); %x= d1
L3=find(NL(:,2)==max(NL(:,2)));  % y=d2
L4=find(NL(:,1)==min(NL(:,1)));  % x=0
n= length(L1);
switch BC_Type
    case 'ss4'     % simply supported on all edges
        
        dofl1= zeros(1,2*n);
        dofl2= zeros(1,2*n);
        dofl3= zeros(1,2*n);
        dofl4= zeros(1,2*n);
        
        for i=1:n
            index1= 2*(i-1)+2;
            index2= index1-1;
            
            dofl1(index1)= 3*L1(i);
            dofl1(index2)= 3*L1(i)-2;
            dofl3(index1)= 3*L3(i);
            dofl3(index2)= 3*L3(i)-2;
            
        end
        for i=1:n
            index1= 2*(i-1)+2;
            index2= index1-1;
            
            dofl2(index1)= 3*L2(i)-1;
            dofl2(index2)= 3*L2(i)-2;
            dofl4(index1)= 3*L4(i)-1;
            dofl4(index2)= 3*L4(i)-2;
        end
        L1UL3 = union(dofl1,dofl3) ;
        L2UL4 = union(dofl2,dofl4) ;
        BC_DoF = union(L1UL3,L2UL4) ;
        
    case 'cc4'
        
        dofL1 = zeros(1,2*n) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*n) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:n
            index1 = 2*(i-1)+2 ;
            index2 = index1-1 ;
            dofL1(index1) = 3*L1(i)-1;
            dofL1(index2) = 3*L1(i)-2 ;
            dofL3(index1) = 3*L3(i)-1 ;
            dofL3(index2) = 3*L3(i)-2 ;
        end
        for i = 1:n
            index1 = 2*(i-1)+2 ;
            index2 = index1-1 ;
            dofL2(index1) = 3*L2(i) ;
            dofL2(index2) = 3*L2(i)-2 ;
            dofL4(index1) = 3*L4(i) ;
            dofL4(index2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        BC_DoF = union(L1UL3,L2UL4) ;
    case 'fp'
        dofL1 = zeros(1,2*n) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*n) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:n
            index1 = 2*(i-1)+2 ;
            index2 = index1-1 ;
            dofL1(index1) = 3*L1(i)-1;
            dofL1(index2) = 3*L1(i)-2 ;
            dofL3(index1) = 3*L3(i)-1 ;
            dofL3(index2) = 3*L3(i)-2 ;
        end
        for i = 1:n
            index1 = 2*(i-1)+2 ;
            index2 = index1-1 ;
            dofL2(index1) = 3*L2(i) ;
            dofL2(index2) = 3*L2(i)-2 ;
            dofL4(index1) = 3*L4(i) ;
            dofL4(index2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        BC_DoF = union(L1UL3,L2UL4) ;
end

end