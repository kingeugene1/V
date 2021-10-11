function [L]= L_distribution(DoF,p,m,EL,NL,da,db,tt2,x_x,y_y)
L21 = zeros( DoF, 1 ); %-- distribtuion matrix
L22 = zeros( DoF, 1 );
L23 = zeros( DoF, 1 );
L24 = zeros( DoF, 1);
for i2=1:p
    for j2=1:m
        l2=EL(((j2-1)*p)+i2,:);
        x_c_2=NL(l2,1)/da;
        y_c_2=NL(l2,2)/db;
        z_c_2=[x_c_2 y_c_2];
        b12=z_c_2(2,:)-z_c_2(1,:);
        b22=z_c_2(3,:)-z_c_2(2,:);
        b32=z_c_2(4,:)-z_c_2(3,:);
        b42=z_c_2(1,:)-z_c_2(4,:);
        B_c2=[b12;b22;b32;b42];
        
        sx_1=(x_x(1,tt2))/da;
        c1 = x_x(1,tt2) > da*(i2-1) & x_x(1,tt2) <= da*i2;
        sx_2=(x_x(2,tt2))/da;
        c2 = x_x(2,tt2) > da*(i2-1) & x_x(2,tt2) <= da*i2;
        sx_3=(x_x(3,tt2))/da;
        c3 = x_x(3,tt2) > da*(i2-1) & x_x(3,tt2) <= da*i2;
        sx_4=(x_x(4,tt2))/da;
        c4 = x_x(4,tt2) > da*(i2-1) & x_x(4,tt2) <= da*i2;
        sy_1=y_y(1,tt2)/db;
        sy_2=y_y(2,tt2)/db;
        sy_3=y_y(3,tt2)/db;
        sy_4=y_y(4,tt2)/db;
        k_c1=[sx_1 sy_1];
        k_c2=[sx_2 sy_2];
        k_c3=[sx_3 sy_3];
        k_c4=[sx_4 sy_4];
        %first point
        a12_1=k_c1-z_c_2(1,:);
        a22_1=k_c1-z_c_2(2,:);
        a32_1=k_c1-z_c_2(3,:);
        a42_1=k_c1-z_c_2(4,:);
        %second point
        a12_2=k_c2-z_c_2(1,:);
        a22_2=k_c2-z_c_2(2,:);
        a32_2=k_c2-z_c_2(3,:);
        a42_2=k_c2-z_c_2(4,:);
        %third point
        a12_3=k_c3-z_c_2(1,:);
        a22_3=k_c3-z_c_2(2,:);
        a32_3=k_c3-z_c_2(3,:);
        a42_3=k_c3-z_c_2(4,:);
        %fourth point
        a12_4=k_c4-z_c_2(1,:);
        a22_4=k_c4-z_c_2(2,:);
        a32_4=k_c4-z_c_2(3,:);
        a42_4=k_c4-z_c_2(4,:);
        
        A_c_21=[a12_1;a22_1;a32_1;a42_1];
        A_c_22=[a12_2;a22_2;a32_2;a42_2];
        A_c_23=[a12_3;a22_3;a32_3;a42_3];
        A_c_24=[a12_4;a22_4;a32_4;a42_4];
        BA_21=B_c2.*A_c_21;
        BA_22=B_c2.*A_c_22;
        BA_23=B_c2.*A_c_23;
        BA_24=B_c2.*A_c_24;
        s_21=BA_21(BA_21~=0);
        s_22=BA_22(BA_22~=0);
        s_23=BA_23(BA_23~=0);
        s_24=BA_24(BA_24~=0);
        c21=all(s_21>0);
        c22=all(s_22>0);
        c23=all(s_23>0);
        c24=all(s_24>0);
        
        %L1
        L21(l2(1)*3-2,:)= L21(l2(1)*3-2,:)+1/4*( c1*c21.*((-sx_1+1)) .*((-sy_1+1)))';
        L21(l2(1)*3-1,:)=L21(l2(1)*3-1,:)+1/4* ( c1*c21.*( (-sx_1+1)).*((-sy_1+1)))';
        L21(l2(1)*3,:)=L21(l2(1)*3,:)+1/4* ( c1*c21.*((sx_1+1)).*((sy_1+1)))';
        L21(l2(2)*3-2,:)= L21(l2(2)*3-2,:)+1/4*( c1*c21.*((-sx_1+1)) .*((-sy_1+1)))';
        L21(l2(2)*3-1,:)=L21(l2(2)*3-1,:)+1/4* ( c1*c21.*( (-sx_1+1)).*((-sy_1+1)))';
        L21(l2(2)*3,:)=L21(l2(2)*3,:)+1/4* ( c1*c21.*((sx_1+1)).*((sy_1+1)))';
        L21(l2(3)*3-2,:)=L21(l2(3)*3-2,:)+1/4*( c1*c21.*((-sx_1+1)) .*((-sy_1+1)))';
        L21(l2(3)*3-1,:)=L21(l2(3)*3-1,:)+1/4* ( c1*c21.*( (-sx_1+1)).*((-sy_1+1)))';
        L21(l2(3)*3,:)=L21(l2(3)*3,:)+1/4* ( c1*c21.*((sx_1+1)).*((sy_1+1)))';
        L21(l2(4)*3-2,:)=L21(l2(4)*3-2,:)+1/4*( c1*c21.*((-sx_1+1)) .*((-sy_1+1)))';
        L21(l2(4)*3-1,:)=L21(l2(4)*3-1,:)+1/4* ( c1*c21.*( (-sx_1+1)).*((-sy_1+1)))';
        L21(l2(4)*3,:)= L21(l2(4)*3,:)+1/4* ( c1*c21.*((sx_1+1)).*((sy_1+1)))';
        
        %L2
        L22(l2(1)*3-2,:)= L22(l2(1)*3-2,:)+1/4*( c2*c22.*((-sx_2+1)) .*((-sy_2+1)))';
        L22(l2(1)*3-1,:)=L22(l2(1)*3-1,:)+1/4* ( c2*c22.*( (-sx_2+1)).*((-sy_2+1)))';
        L22(l2(1)*3,:)=L22(l2(1)*3,:)+1/4* ( c2*c22.*((sx_2+1)).*((sy_2+1)))';
        L22(l2(2)*3-2,:)= L22(l2(2)*3-2,:)+1/4*( c2*c22.*((-sx_2+1)) .*((-sy_2+1)))';
        L22(l2(2)*3-1,:)=L22(l2(2)*3-1,:)+1/4* ( c2*c22.*( (-sx_2+1)).*((-sy_2+1)))';
        L22(l2(2)*3,:)=L22(l2(2)*3,:)+1/4* ( c2*c22.*((sx_2+1)).*((sy_2+1)))';
        L22(l2(3)*3-2,:)=L22(l2(3)*3-2,:)+1/4*( c2*c22.*((-sx_2+1)) .*((-sy_2+1)))';
        L22(l2(3)*3-1,:)=L22(l2(3)*3-1,:)+1/4* (c2*c22.*( (-sx_2+1)).*((-sy_2+1)))';
        L22(l2(3)*3,:)=L22(l2(3)*3,:)+1/4* ( c2*c22.*((sx_2+1)).*((sy_2+1)))';
        L22(l2(4)*3-2,:)=L22(l2(4)*3-2,:)+1/4*( c2*c22.*((-sx_2+1)) .*((-sy_2+1)))';
        L22(l2(4)*3-1,:)=L22(l2(4)*3-1,:)+1/4* ( c2*c22.*( (-sx_2+1)).*((-sy_2+1)))';
        L22(l2(4)*3,:)= L22(l2(4)*3,:)+1/4* ( c2*c22.*((sx_2+1)).*((sy_2+1)))';
        %L3
     L23(l2(1)*3-2,:)= L23(l2(1)*3-2,:)+1/4*( c3*c23.*((-sx_3+1)) .*((-sy_3+1)))';
        L23(l2(1)*3-1,:)=L23(l2(1)*3-1,:)+1/4* ( c3*c23.*( (-sx_3+1)).*((-sy_3+1)))';
        L23(l2(1)*3,:)=L23(l2(1)*3,:)+1/4* ( c3*c23.*((sx_3+1)).*((sy_3+1)))';
        L23(l2(2)*3-2,:)= L23(l2(2)*3-2,:)+1/4*( c3*c23.*((-sx_3+1)) .*((-sy_3+1)))';
       L23(l2(2)*3-1,:)=L23(l2(2)*3-1,:)+1/4* ( c3*c23.*( (-sx_3+1)).*((-sy_3+1)))';
        L23(l2(2)*3,:)=L23(l2(2)*3,:)+1/4* ( c3*c23.*((sx_3+1)).*((sy_3+1)))';
        L23(l2(3)*3-2,:)=L23(l2(3)*3-2,:)+1/4*( c3*c23.*((-sx_3+1)) .*((-sy_3+1)))';
       L23(l2(3)*3-1,:)=L23(l2(3)*3-1,:)+1/4* (c3*c23.*( (-sx_3+1)).*((-sy_3+1)))';
       L23(l2(3)*3,:)=L23(l2(3)*3,:)+1/4* ( c3*c23.*((sx_3+1)).*((sy_3+1)))';
       L23(l2(4)*3-2,:)=L23(l2(4)*3-2,:)+1/4*( c3*c23.*((-sx_3+1)) .*((-sy_3+1)))';
       L23(l2(4)*3-1,:)=L23(l2(4)*3-1,:)+1/4* ( c3*c23.*( (-sx_3+1)).*((-sy_3+1)))';
        L23(l2(4)*3,:)= L23(l2(4)*3,:)+1/4* ( c3*c23.*((sx_3+1)).*((sy_3+1)))';
        %L4
      L24(l2(1)*3-2,:)=L24(l2(1)*3-2,:)+1/4*( c4*c24.*((-sx_4+1)) .*((-sy_4+1)))';
      L24(l2(1)*3-1,:)=L24(l2(1)*3-1,:)+1/4* ( c4*c24.*( (-sx_4+1)).*((-sy_4+1)))';
      L24(l2(1)*3,:)=L24(l2(1)*3,:)+1/4* ( c4*c24.*((sx_4+1)).*((sy_4+1)))';
      L24(l2(2)*3-2,:)= L24(l2(2)*3-2,:)+1/4*( c4*c24.*((-sx_4+1)) .*((-sy_4+1)))';
       L24(l2(2)*3-1,:)=L24(l2(2)*3-1,:)+1/4* ( c4*c24.*( (-sx_4+1)).*((-sy_4+1)))';
      L24(l2(2)*3,:)=L24(l2(2)*3,:)+1/4* ( c4*c24.*((sx_4+1)).*((sy_4+1)))';
       L24(l2(3)*3-2,:)=L24(l2(3)*3-2,:)+1/4*( c4*c24.*((-sx_4+1)) .*((-sy_4+1)))';
     L24(l2(3)*3-1,:)=L24(l2(3)*3-1,:)+1/4* (c4*c24.*( (-sx_4+1)).*((-sy_4+1)))';
      L24(l2(3)*3,:)=L24(l2(3)*3,:)+1/4* ( c4*c24.*((sx_4+1)).*((sy_4+1)))';
      L24(l2(4)*3-2,:)=L24(l2(4)*3-2,:)+1/4*( c4*c24.*((-sx_4+1)) .*((-sy_4+1)))';
      L24(l2(4)*3-1,:)=L24(l2(4)*3-1,:)+1/4* ( c4*c24.*( (-sx_4+1)).*((-sy_4+1)))';
        L24(l2(4)*3,:)= L24(l2(4)*3,:)+1/4* ( c4*c24.*((sx_4+1)).*((sy_4+1)))';
    end
end
L=[L21 L22 L23 L24];
end