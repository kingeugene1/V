function [Kb,Mb]= M_K_mat(NoE,NPE,PD,NL,EL,D,rhA_data)
Dof_E= 3*size(EL,2);
N_Dof=3;
DoF=N_Dof*size(NL,1);
ke= zeros(Dof_E,Dof_E);
me= zeros(Dof_E,Dof_E);
Kb= zeros(DoF,DoF);                                     % initialization of stiffness, Force and Mass matrices
Mb= zeros(DoF,DoF);
for i= 1:NoE
    index= EL(i,1:NPE);
    x= zeros(NPE,PD);
    
    for j= 1:NPE
        
        x(j,:)= NL(index(j),1:PD);
    end
    
    x=x';
    xx=x(1,:);                                            %    x-coordinates of elements
    yy=x(2,:);                                            % y-coordinates of elements
    %gaussian intergral
    %    GP= 2;
    GP=4;
    [X,weight]=lgwt(GP,-1,1);
    %    [X,weight]=gausspoint();
    for ii= 1:GP
        xi= X(ii);
        wx= weight(ii);
        for j=1:GP
            eta= X(j);
            wy=weight(j);
            [~,dNdxi, dNdeta] = shape_function1(xi,eta);   % Shape function and their first derivatives
            [det_J, ~]= Jacobian(dNdxi,dNdeta,xx,yy);   % inverse jacobian and determinant jacobian
            %           [dNdx, dNdy]= dN_global(NPE,dNdxi,dNdeta,inv_J);
            %calculate first derivatives of shape function with respect to global coordinates (x,y)
            
            H= Hfunction(xi,eta);
            m1=M_vector(NPE,H);
            % element stiffness matrix
            [A,G,B]= ke_Matrix(xi,eta,xx,yy);
            ke=ke+det_J*wx*wy*(G')*(A')*D*B*G;
            
            me=me+ rhA_data(i)*m1*wx*wy*det_J;
        end
        
    end
    
    % ASSEMBLE STIFFNESS AND FORCE MATRIX
    idx= elementdof(index,NPE,N_Dof);                            % extract dofs for a given element as a list
    
    [Kb,Mb]= Assemble(Kb,Mb,ke,idx,me);
    
end
end