tic
clear;  close;    clc;
d1= 5;                                                          % x-length [m]
d2= 4;                                                          % y- length [m]
p= 10;                                                          % number of elements in x-direction
m=10;                                                          % number of elements in y-direction
da= d1/p;                                                     % increment size in x-direction
db=d2/m;
% dy= d2/m;                                                 % increment size in y-direction
NoE= p*m;                                                   % Total number of elements
NoN= (p+1)*(m+1);                                     % Total number of nodes
NPE= 4;                                                       % Number of elements per node
PD= 2;                                                         % Problem dimension
N_Dof = 3;                                                   % number of degree of freedom per node
DoF = N_Dof*NoN;                                      % Total degree of freedom of the system
Dof_E= N_Dof*NPE;                                    % degree of freedom per element
 E= 2.56e15;                                                   % Plate young's modulus
% E= 2.5e6;
% E = 2.5e8; %-- flexural rigidity 
thich= 0.3;                                                        % Plate thickness
nu=0.3;
D0= E*thich^3/(12*(1-nu^2));                           % Flexural rigidity of plate
rhA=4400;
g=-9.81;                                           %[m/s^2]
%%  vehicle Parameters
fprintf('setting vehicle parameters...\n')
ms         =     9000;                                             %-- sprung-mass (kg)
mu1      =      400;                                                 %-- unsprung-mass (kg)
mu2      =      500;                                                 %-- unsprung-mass (kg)
mu3      =      400;
mu4      =      500;
cs1 = 2000;  %-- sprung-damping (kg/s)
cs2 = 2000;  %-- sprung-damping (kg/s)
cs3 = 2000;  %-- sprung-damping (kg/s)
cs4 = 2000;  %-- sprung-damping (kg/s)

ks1= 456000;  %-- sprung-stifness (N/m)=(kg/s/s)
ks2= 4500;  %-- sprung-stifness (N/m)=(kg/s/s)
ks3= 4500;  %-- sprung-stifness (N/m)=(kg/s/s)
ks4= 4500;  %-- sprung-stifness (N/m)=(kg/s/s)

ku1 = 60000; %-- unsprung-stifness (N/m)=(kg/s/s)
ku2 = 60000; %-- unsprung-stifness (N/m)=(kg/s/s)
ku3 = 60000; %-- unsprung-stifness (N/m)=(kg/s/s)
ku4 = 60000; %-- unsprung-stifness (N/m)=(kg/s/s)

% cs          =      2000;                                                %-- sprung-damping (kg/s)
% ks          =      4500;                                                %-- sprung-stifness (N/m)=(kg/s/s)
% ku          =      20000;                                             %-- unsprung-stifness (N/m)=(kg/s/s)
% d            =      1.75;                                                 %-- distance from G to s1 (m)
% L= 0.9;
dd1            =      1.75;                                                 %-- distance from G to s1 (m)
dd2= 1.75;                                                         %[m]
L= 1.8;  

v   = 10;                                                               %-- speed (m/s)
%v=18 km/h;                                                                  % change (m/s)
dt  = 0.01;                                                             %-- time increment
% T   = 15;                                                                  %-- end time
x0  = 10;                                                                  %-- start position of G(m)
xb0= 20;    % Entrance position
% yb0= 0.6;
%-- calculate input
fprintf('calculating vehicle input...\n')
R = load('profile-1.txt');                                 %-- load: R(x)
T=length(R);

x_1 = (0:0.10:(length(R)-1)*0.10)-xb0;
% y_1 = (0:0.10:(length(R)-1)*0.10)-yb0;
t1=60;

t2 = 0:dt:t1;
n= length(t2);
t = 0:dt:T;                                                                  %-- time, dt=0.001(s)
% y=zeros(length(t),1);
y_y= [ 0.68*ones(1,n); 2.48*ones(1,n);0.68*ones(1,n);2.48*ones(1,n)]; 

x_x= [x0+v*t2+dd1 ;x0+v*t2+dd1; x0+v*t2-dd2 ;x0+v*t2-dd2]-xb0;    %-- vehicle position: x(t)

rx = spline(x_1,R,x_x);                                                     %-- spline interpolation: R(x) to r(t)=R(x(t))
% plot(t2,rx(1,:),'r')
% xlim([0 10])
% ylim([-0.01 0.001])
%-- simulate output
fprintf('vehicle vibration simulation: newmark-beta method...\n')
Gama = 1/2;                                                                        %-- gammma coefficient for newmark-beta method
Beta = 1/4;                                                                         %-- beta coefficient for newmark-beta method
% Iy = 0.5*d*ms;
% Ix =0.5* L*ms;

Iy =( dd1*dd2*ms)/(dd1+dd2);
Ix =0.5* L*ms;
% Mv = [0.25*ms          0.25*ms       0.25*ms        0.25*ms          0         0        0          0;
%          0.25*Iy/d       -0.25*Iy/d     0.25*Iy/d    -0.25*Iy/d        0         0        0          0;
%         -0.25*Ix/L      -0.25*Ix/L     0.25*Ix/L     0.25*Ix/L        0         0        0          0;
%          0                      0                    0                   0                      0         0        0         0;
%          0                      0                    0                   0                       mu1   0        0         0;
%          0                      0                    0                   0                       0         mu2  0         0; 
%          0                      0                    0                   0                       0         0        mu3   0;
%          0                      0                    0                   0                       0         0        0        mu4];
Mv = [(L*dd2)/((dd1+dd2 )*(2*L))*ms           (L*dd1)/((dd1+dd2 )*(2*L))*ms       (L*dd2)/((dd1+dd2 )*(2*L))*ms       (L*dd1)/((dd1+dd2 )*(2*L))*ms      0      0    0     0;
         0.5*Iy/(dd1+dd2)     -0.5*Iy/(dd1+dd2)     0.5*Iy/(dd1+dd2)    -0.5*Iy/(dd1+dd2)    0      0    0     0;
        -dd2*Ix/(L*(dd1+dd2))     -dd1*Ix/(L*(dd1+dd2))     dd2*Ix/(L*(dd1+dd2))     dd1*Ix/(L*(dd1+dd2))    0      0    0     0;
         0              0             0             0            0      0    0     0;
         0              0             0             0            mu1    0    0     0;
         0              0             0             0            0      mu2  0     0; 
         0              0             0             0            0      0    mu3   0;
         0              0             0             0            0      0    0     mu4];       

% Cv = [cs           cs           cs             cs         -cs       -cs       -cs          -cs;
%          d* cs     -d* cs      d* cs       -d*cs     -d*cs    d*cs    -d*cs       d* cs;
%          -L*cs    -L*cs       L*cs         L*cs      L*cs    L*cs    -L*cs     -L*cs  ;
%           0           0             0               0           0         0           0             0;
%          -cs         0             0               0           cs        0           0             0;
%           0         -cs            0               0            0         cs          0            0;
%           0          0            -cs              0            0          0           cs          0;
%          0       0       0             -cs           0          0           0            cs];
 Cv = [cs1       cs2      cs3        cs4      -cs1      -cs2     -cs3    -cs4;
      dd1* cs1   -dd2* cs2   dd1* cs3    -dd2*cs4    -dd1*cs1    dd2*cs2   -dd1*cs3    dd2* cs4;
     -L*cs1    -L*cs2    L*cs3      L*cs4     L*cs1    L*cs2    -L*cs3    -L*cs4;
      0        0       0         0        0       0        0        0;
     -cs1       0       0         0        cs1      0        0        0;
      0       -cs2      0         0        0       cs2       0        0;
      0        0      -cs3        0        0       0        cs3       0;
      0        0       0        -cs4       0       0        0        cs4]; 

% Kv =[ ks          ks        ks          ks        -ks       -ks         -ks         -ks;
%         d*ks     -d*ks     d*ks    -d*ks     -d*ks    d*ks     -d*ks       d*ks;
%         -L*ks    -L*ks    L*ks     L*ks       L*ks    L*ks    -L*ks      -L*ks;
%          1            -1          -1           1             0         0           0             0;
%        -ks            0         0           0             ks+ku  0           0            0;
%         0             -ks       0           0             0           ks+ku   0            0;
%         0              0        -ks         0             0           0           ks+ku    0;
%         0              0         0          -ks           0           0           0             ks+ku];
 Kv =[ ks1       ks2     ks3       ks4     -ks1     -ks2     -ks3     -ks4;
      dd1*ks1    -dd2*ks2   dd1*ks3    -dd2*ks4   -dd1*ks1    dd2*ks2   -dd1*ks3    dd2*ks4;
     -L*ks1    -L*ks2   L*ks3     L*ks4    L*ks1    L*ks2   -L*ks3   -L*ks4;
      1       -1     -1        1       0       0       0       0;
     -ks1       0      0        0       ks1+ku1   0       0       0;
      0       -ks2     0        0       0       ks2+ku2   0       0;
      0        0     -ks3       0       0       0       ks3+ku3   0;
      0        0      0       -ks4      0       0       0       ks4+ku4];
  
%     Mp = [0.25*ms          0                    0                  0                      mu1      0            0          0;
%               0                      0.25*ms         0                  0                      0           mu2       0          0;
%               0                      0                    0.25*ms        0                     0           0            mu3      0;
%               0                      0                    0                   0.25*ms          0            0           0           mu4]; 
 Mp = [(L*dd2)/((dd1+dd2 )*(2*L))*ms          0                    0                  0                      mu1      0            0          0;
              0                      (L*dd1)/((dd1+dd2 )*(2*L))*ms         0                  0                      0           mu2       0          0;
              0                      0                    (L*dd2)/((dd1+dd2 )*(2*L))*ms        0                     0           0            mu3      0;
              0                      0                    0                   (L*dd1)/((dd1+dd2 )*(2*L))*ms          0            0           0           mu4];                    
%variables for newmark************************************
a=0.200; b= 0.005200;

% Fv= [zeros(4,n); [ku1 0 0 0; 0 ku2 0 0; 0 0 ku3 0;0 0 0 ku4]*rx];
%     [ ~, z, dz, ddz ] = newmark_beta( Mv,Cv,Kv,Fv,1/1000,zeros(8,3), Gama, Beta );
%     %-- bridge system
%     P = Mp*(repmat(g,8,n)-ddz);
%  [f]= moving_load(DoF,n,P,x_x,y_y);
 
%%%%%*******************Bridge System*************************************************
fprintf('bridge matrices...\n')
[NL,EL]=uniform_mesh(d1,d2,p,m);             % this function return the node list as coordinates (x,y) and element list
[D]=Dmatrix(E,thich,nu);                                   % this function return D matrix of isoparametric material

Kb= zeros(DoF,DoF);                                     % initialization of stiffness, Force and Mass matrices
% Force= zeros(DoF,n);
Mb= zeros(DoF,DoF);

% extract x and y coordinate

 for i= 1:NoE
      index= EL(i,1:NPE);
     x= zeros(NPE,PD);

     for j= 1:NPE
       
         x(j,:)= NL(index(j),1:PD);
     end
     
      x=x';
     xx=x(1,:);                                            %    x-coordinates of elements
     yy=x(2,:);                                            % y-coordinates of elements
     

  % Initialize element Matrices
   ke= zeros(Dof_E,Dof_E);
   kke= zeros(3,Dof_E);
   me= zeros(Dof_E,Dof_E);
%    f= zeros(Dof_E,1);
%    p0=-1;
   %gaussian intergral
   
%    GP= 2;
   GP=4;
   
   [X,weight]=gausspoint();
  for ii= 1:GP
      xi= X(ii,1);
      wx= weight(ii,1);
      for j=1:GP
          eta= X(j,2);
          wy=weight(j,2);
          [N,dNdxi, dNdeta] = shape_function (xi,eta);   % Shape function and their first derivatives
          [det_J, inv_J]= Jacobian(dNdxi,dNdeta,xx,yy);   % inverse jacobian and determinant jacobian
          
          %calculate first derivatives of shape function with respect to global coordinates (x,y)
          [dNdx, dNdy]= dN_global(NPE,dNdxi,dNdeta,inv_J);
          H= Hfunction(xi,eta);
          m1=M_vector(NPE,H);
          % element stiffness matrix
%           fe= Force_vector(NPE,N,p0);
      
          [A,B]= AB_Matrix(xi,eta,xx,yy);
          ke=ke+det_J*wx*wy*A'*B'*D*B*A;
%           f=f+fe*wx*wy*det_J;
          me=me+ rhA*m1*wx*wy*det_J;
      end
      
  end
  
  % ASSEMBLE STIFFNESS AND FORCE MATRIX
     idx= elementdof(index,NPE,N_Dof);                            % extract dofs for a given element as a list
     
     [Kb,Mb]= Assemble(Kb,Mb,ke,idx,me);
     
 end
 Cb= a*Mb+ b*Kb;    % Rayleigh damping
 
 
 %%
 
 
 ddy=zeros(DoF,n);
 dy=zeros(DoF,n);
 y=zeros(DoF,n);
%  mem= zeros(DoF,n); %-- memory
 
 fprintf('VBI-system simulation: newmark-beta method...\n')
 for count=1:5
     yt = zeros(4,length(t2));
     for tt=1:length(t2)
         L11 = zeros( DoF, 1 ); %-- distribtuion matrix
         L12 = zeros( DoF, 1 );
         L13 = zeros( DoF, 1 );
         L14 = zeros( DoF, 1);
        for il=1:p
            for jl=1:m
                %%%%
                  l2=EL(((jl-1)*p)+il,:);
                 x_c_2=NL(l2,1)/da;
                 y_c_2=NL(l2,2)/db;
                 z_c_2=[x_c_2 y_c_2];
                 b12=z_c_2(2,:)-z_c_2(1,:);
                 b22=z_c_2(3,:)-z_c_2(2,:);
                 b32=z_c_2(4,:)-z_c_2(3,:);
                 b42=z_c_2(1,:)-z_c_2(4,:);
                 B_c2=[b12;b22;b32;b42];
%                 B_c2= repmat(B_c2,4,1);
                 sx_1=(x_x(1,tt)-da*(il-1))/da;
                 sx_2=(x_x(2,tt)-da*(il-1))/da;
                 sx_3=(x_x(3,tt)-da*(il-1))/da;
                 sx_4=(x_x(4,tt)-da*(il-1))/da;
                 sy_1=y_y(1,tt)/db;
                 sy_2=y_y(2,tt)/db;
                 sy_3=y_y(3,tt)/db;
                 sy_4=y_y(4,tt)/db;
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
                L11(l2(1)*3-2,:)= L11(l2(1)*3-2,:)+1/16*( c21.*(sx_1.^3-3*sx_1+2) .*(sy_1.^3-3*sy_1+2))';
                L11(l2(1)*3-1,:)=L11(l2(1)*3-1,:)+1/16* ( c21.*( sx_1.^3-sx_1.^2-sx_1+1).*(sy_1.^3-3*sy_1+2))';
                L11(l2(1)*3,:)=L11(l2(1)*3,:)+1/16* ( c21.*(sx_1.^3-3*sx_1+2).*(sy_1.^3-sy_1.^2-sy_1+1 ))';
                L11(l2(2)*3-2,:)= L11(l2(2)*3-2,:)-1/16* ( c21.*(sx_1.^3-3*sx_1-2).*(sy_1.^3-3*sy_1+2) )';
                L11(l2(2)*3-1,:)=L11(l2(2)*3-1,:)+1/16*(c21.*(sx_1.^3-3*sx_1.^2-sx_1-1).*(sy_1.^3-3*sy_1+2))';
                L11(l2(2)*3,:)=L11(l2(2)*3,:)-1/16*(c21.*((sx_1.^3-3*sx_1-2).*(sy_1.^3-sy_1.^2-sy_1+1 )))';
                L11(l2(3)*3-2,:)=L11(l2(3)*3-2,:)+1/16*(c21.*(sx_1.^3-3*sx_1-2).*(sy_1.^3-3*sy_1-2))';
                L11(l2(3)*3-1,:)=L11(l2(3)*3-1,:)-1/16*(c21.*(sx_1.^3-3*sx_1.^2-sx_1-1).*(sy_1.^3-3*sy_1-2))';
                L11(l2(3)*3,:)=L11(l2(3)*3,:)-1/16*(c21.*(sx_1.^3-3*sx_1-2).*(sy_1.^3-3*sy_1.^2-sy_1-1))';
                L11(l2(4)*3-2,:)=L11(l2(4)*3-2,:)-1/16*(c21.*(sx_1.^3-3*sx_1+2).*(sy_1.^3-3*sy_1-2))';
                L11(l2(4)*3-1,:)=L11(l2(4)*3-1,:)-1/16*(c21.*(sx_1.^3-sx_1.^2-sx_1+1 ).*(sy_1.^3-3*sy_1-2))';
                L11(l2(4)*3,:)= L11(l2(4)*3,:)+1/16*(c21.*(sx_1.^3-3*sx_1+2).*(sy_1.^3-3*sy_1.^2-sy_1-1))';
                %L2
                L12(l2(1)*3-2,:)= L12(l2(1)*3-2,:)+1/16*( c22.*(sx_2.^3-3*sx_2+2) .*(sy_2.^3-3*sy_2+2))';
                L12(l2(1)*3-1,:)=L12(l2(1)*3-1,:)+1/16* ( c22.*( sx_2.^3-sx_2.^2-sx_2+1).*(sy_2.^3-3*sy_2+2))';
                L12(l2(1)*3,:)=L12(l2(1)*3,:)+1/16* ( c22.*(sx_2.^3-3*sx_2+2).*(sy_2.^3-sy_2.^2-sy_2+1 ))';
                L12(l2(2)*3-2,:)= L12(l2(2)*3-2,:)-1/16* ( c22.*(sx_2.^3-3*sx_2-2).*(sy_2.^3-3*sy_2+2) )';
                L12(l2(2)*3-1,:)=L12(l2(2)*3-1,:)+1/16*(c22.*(sx_2.^3-3*sx_2.^2-sx_1-1).*(sy_2.^3-3*sy_2+2))';
                L12(l2(2)*3,:)=L12(l2(2)*3,:)-1/16*(c22.*((sx_2.^3-3*sx_2-2).*(sy_2.^3-sy_2.^2-sy_2+1 )))';
                L12(l2(3)*3-2,:)=L12(l2(3)*3-2,:)+1/16*(c22.*(sx_2.^3-3*sx_2-2).*(sy_2.^3-3*sy_2-2))';
                L12(l2(3)*3-1,:)=L12(l2(3)*3-1,:)-1/16*(c22.*(sx_2.^3-3*sx_2.^2-sx_2-1).*(sy_2.^3-3*sy_2-2))';
               L12(l2(3)*3,:)=L12(l2(3)*3,:)-1/16*(c22.*(sx_2.^3-3*sx_2-2).*(sy_2.^3-3*sy_2.^2-sy_2-1))';
                L12(l2(4)*3-2,:)=L12(l2(4)*3-2,:)-1/16*(c22.*(sx_2.^3-3*sx_2+2).*(sy_2.^3-3*sy_2-2))';
                L12(l2(4)*3-1,:)=L12(l2(4)*3-1,:)-1/16*(c22.*(sx_2.^3-sx_2.^2-sx_2+1 ).*(sy_2.^3-3*sy_2-2))';
                L12(l2(4)*3,:)=L12(l2(4)*3,:)+1/16*(c22.*(sx_2.^3-3*sx_2+2).*(sy_2.^3-3*sy_2.^2-sy_2-1))';
                %L3
                 L13(l2(1)*3-2,:)= L13(l2(1)*3-2,:)+1/16*( c23.*(sx_3.^3-3*sx_3+2) .*(sy_3.^3-3*sy_3+2))';
               L13(l2(1)*3-1,:)=L13(l2(1)*3-1,:)+1/16* ( c23.*( sx_3.^3-sx_3.^2-sx_3+1).*(sy_3.^3-3*sy_3+2))';
                L13(l2(1)*3,:)=L13(l2(1)*3,:)+1/16* ( c23.*(sx_3.^3-3*sx_3+2).*(sy_3.^3-sy_3.^2-sy_3+1 ))';
                L13(l2(2)*3-2,:)= L13(l2(2)*3-2,:)-1/16* ( c23.*(sx_3.^3-3*sx_3-2).*(sy_3.^3-3*sy_3+2) )';
                L13(l2(2)*3-1,:)=L13(l2(2)*3-1,:)+1/16*(c23.*(sx_3.^3-3*sx_3.^2-sx_3-1).*(sy_3.^3-3*sy_3+2))';
                L13(l2(2)*3,:)=L13(l2(2)*3,:)-1/16*(c23.*((sx_3.^3-3*sx_3-2).*(sy_3.^3-sy_3.^2-sy_3+1 )))';
                L13(l2(3)*3-2,:)=L13(l2(3)*3-2,:)+1/16*(c23.*(sx_3.^3-3*sx_3-2).*(sy_3.^3-3*sy_3-2))';
                L13(l2(3)*3-1,:)=L13(l2(3)*3-1,:)-1/16*(c23.*(sx_3.^3-3*sx_3.^2-sx_3-1).*(sy_3.^3-3*sy_3-2))';
                L13(l2(3)*3,:)=L13(l2(3)*3,:)-1/16*(c23.*(sx_3.^3-3*sx_3-2).*(sy_3.^3-3*sy_3.^2-sy_3-1))';
                L13(l2(4)*3-2,:)=L13(l2(4)*3-2,:)-1/16*(c23.*(sx_3.^3-3*sx_3+2).*(sy_3.^3-3*sy_3-2))';
                L13(l2(4)*3-1,:)=L13(l2(4)*3-1,:)-1/16*(c23.*(sx_3.^3-sx_3.^2-sx_3+1 ).*(sy_3.^3-3*sy_3-2))';
                L13(l2(4)*3,:)= L13(l2(4)*3,:)+1/16*(c23.*(sx_3.^3-3*sx_3+2).*(sy_3.^3-3*sy_3.^2-sy_3-1))';
                %L4
                L14(l2(1)*3-2,:)= L14(l2(1)*3-2,:)+1/16*( c24.*(sx_4.^3-3*sx_4+2) .*(sy_4.^3-3*sy_4+2))';
                L14(l2(1)*3-1,:)=L14(l2(1)*3-1,:)+1/16* ( c24.*( sx_4.^3-sx_4.^2-sx_4+1).*(sy_4.^3-3*sy_4+2))';
                L14(l2(1)*3,:)=L14(l2(1)*3,:)+1/16* ( c24.*(sx_4.^3-3*sx_4+2).*(sy_4.^3-sy_4.^2-sy_4+1 ))';
                L14(l2(2)*3-2,:)= L14(l2(2)*3-2,:)-1/16* ( c24.*(sx_4.^3-3*sx_4-2).*(sy_4.^3-3*sy_4+2) )';
                L14(l2(2)*3-1,:)=L14(l2(2)*3-1,:)+1/16*(c24.*(sx_4.^3-3*sx_4.^2-sx_4-1).*(sy_4.^3-3*sy_4+2))';
                L14(l2(2)*3,:)=L14(l2(2)*3,:)-1/16*(c24.*((sx_4.^3-3*sx_4-2).*(sy_4.^3-sy_4.^2-sy_4+1 )))';
                L14(l2(3)*3-2,:)=L14(l2(3)*3-2,:)+1/16*(c24.*(sx_4.^3-3*sx_4-2).*(sy_4.^3-3*sy_4-2))';
                L14(l2(3)*3-1,:)=L14(l2(3)*3-1,:)-1/16*(c24.*(sx_4.^3-3*sx_4.^2-sx_4-1).*(sy_4.^3-3*sy_4-2))';
               L14(l2(3)*3,:)=L14(l2(3)*3,:)-1/16*(c24.*(sx_4.^3-3*sx_4-2).*(sy_4.^3-3*sy_4.^2-sy_4-1))';
                L14(l2(4)*3-2,:)=L14(l2(4)*3-2,:)-1/16*(c24.*(sx_4.^3-3*sx_4+2).*(sy_4.^3-3*sy_4-2))';
                L14(l2(4)*3-1,:)=L14(l2(4)*3-1,:)-1/16*(c24.*(sx_4.^3-sx_4.^2-sx_1+1 ).*(sy_4.^3-3*sy_4-2))';
                L14(l2(4)*3,:)= L14(l2(4)*3,:)+1/16*(c24.*(sx_4.^3-3*sx_4+2).*(sy_4.^3-3*sy_4.^2-sy_4-1))';
            end
        end
         yt(:,tt)  = [L11 L12 L13 L14]'*y(:,tt);
     end
     
 
    u = rx+yt;
  Fv= [zeros(4,n); [ku1 0 0 0; 0 ku2 0 0; 0 0 ku3 0;0 0 0 ku4]*rx];
    [ ~, z, dz, ddz ] = newmark_beta( Mv,Cv,Kv,Fv,1/1000,zeros(8,3), Gama, Beta );
    %-- bridge system
    P = Mp*(repmat(g,8,n)-ddz);
    Fb= zeros(DoF,n);
    for tt2=1:length(t2)
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
%                 B_c2= repmat(B_c2,4,1);
                 sx_1=(x_x(1,tt2)-da*(i2-1))/da;
                 sx_2=(x_x(2,tt2)-da*(i2-1))/da;
                 sx_3=(x_x(3,tt2)-da*(i2-1))/da;
                 sx_4=(x_x(4,tt2)-da*(i2-1))/da;
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
                L21(l2(1)*3-2,:)= L21(l2(1)*3-2,:)+1/16*( c21.*(sx_1.^3-3*sx_1+2) .*(sy_1.^3-3*sy_1+2))';
                L21(l2(1)*3-1,:)=L21(l2(1)*3-1,:)+1/16* ( c21.*( sx_1.^3-sx_1.^2-sx_1+1).*(sy_1.^3-3*sy_1+2))';
                L21(l2(1)*3,:)=L21(l2(1)*3,:)+1/16* ( c21.*(sx_1.^3-3*sx_1+2).*(sy_1.^3-sy_1.^2-sy_1+1 ))';
                L21(l2(2)*3-2,:)= L21(l2(2)*3-2,:)-1/16* ( c21.*(sx_1.^3-3*sx_1-2).*(sy_1.^3-3*sy_1+2) )';
                L21(l2(2)*3-1,:)=L21(l2(2)*3-1,:)+1/16*(c21.*(sx_1.^3-3*sx_1.^2-sx_1-1).*(sy_1.^3-3*sy_1+2))';
                L21(l2(2)*3,:)=L21(l2(2)*3,:)-1/16*(c21.*((sx_1.^3-3*sx_1-2).*(sy_1.^3-sy_1.^2-sy_1+1 )))';
                L21(l2(3)*3-2,:)=L21(l2(3)*3-2,:)+1/16*(c21.*(sx_1.^3-3*sx_1-2).*(sy_1.^3-3*sy_1-2))';
                L21(l2(3)*3-1,:)=L21(l2(3)*3-1,:)-1/16*(c21.*(sx_1.^3-3*sx_1.^2-sx_1-1).*(sy_1.^3-3*sy_1-2))';
                L21(l2(3)*3,:)=L21(l2(3)*3,:)-1/16*(c21.*(sx_1.^3-3*sx_1-2).*(sy_1.^3-3*sy_1.^2-sy_1-1))';
                L21(l2(4)*3-2,:)=L21(l2(4)*3-2,:)-1/16*(c21.*(sx_1.^3-3*sx_1+2).*(sy_1.^3-3*sy_1-2))';
                L21(l2(4)*3-1,:)=L21(l2(4)*3-1,:)-1/16*(c21.*(sx_1.^3-sx_1.^2-sx_1+1 ).*(sy_1.^3-3*sy_1-2))';
                L21(l2(4)*3,:)= L21(l2(4)*3,:)+1/16*(c21.*(sx_1.^3-3*sx_1+2).*(sy_1.^3-3*sy_1.^2-sy_1-1))';
                %L2
                L22(l2(1)*3-2,:)= L22(l2(1)*3-2,:)+1/16*( c22.*(sx_2.^3-3*sx_2+2) .*(sy_2.^3-3*sy_2+2))';
                L22(l2(1)*3-1,:)=L22(l2(1)*3-1,:)+1/16* ( c22.*( sx_2.^3-sx_2.^2-sx_2+1).*(sy_2.^3-3*sy_2+2))';
                L22(l2(1)*3,:)=L22(l2(1)*3,:)+1/16* ( c22.*(sx_2.^3-3*sx_2+2).*(sy_2.^3-sy_2.^2-sy_2+1 ))';
                L22(l2(2)*3-2,:)= L22(l2(2)*3-2,:)-1/16* ( c22.*(sx_2.^3-3*sx_2-2).*(sy_2.^3-3*sy_2+2) )';
                L22(l2(2)*3-1,:)=L22(l2(2)*3-1,:)+1/16*(c22.*(sx_2.^3-3*sx_2.^2-sx_1-1).*(sy_2.^3-3*sy_2+2))';
                L22(l2(2)*3,:)=L22(l2(2)*3,:)-1/16*(c22.*((sx_2.^3-3*sx_2-2).*(sy_2.^3-sy_2.^2-sy_2+1 )))';
                L22(l2(3)*3-2,:)=L22(l2(3)*3-2,:)+1/16*(c22.*(sx_2.^3-3*sx_2-2).*(sy_2.^3-3*sy_2-2))';
                L22(l2(3)*3-1,:)=L22(l2(3)*3-1,:)-1/16*(c22.*(sx_2.^3-3*sx_2.^2-sx_2-1).*(sy_2.^3-3*sy_2-2))';
                L22(l2(3)*3,:)=L22(l2(3)*3,:)-1/16*(c22.*(sx_2.^3-3*sx_2-2).*(sy_2.^3-3*sy_2.^2-sy_2-1))';
                L22(l2(4)*3-2,:)=L22(l2(4)*3-2,:)-1/16*(c22.*(sx_2.^3-3*sx_2+2).*(sy_2.^3-3*sy_2-2))';
                L22(l2(4)*3-1,:)=L22(l2(4)*3-1,:)-1/16*(c22.*(sx_2.^3-sx_2.^2-sx_2+1 ).*(sy_2.^3-3*sy_2-2))';
                L22(l2(4)*3,:)= L22(l2(4)*3,:)+1/16*(c22.*(sx_2.^3-3*sx_2+2).*(sy_2.^3-3*sy_2.^2-sy_2-1))';
                %L3
                 L23(l2(1)*3-2,:)= L23(l2(1)*3-2,:)+1/16*( c23.*(sx_3.^3-3*sx_3+2) .*(sy_3.^3-3*sy_3+2))';
                L23(l2(1)*3-1,:)=L23(l2(1)*3-1,:)+1/16* ( c23.*( sx_3.^3-sx_3.^2-sx_3+1).*(sy_3.^3-3*sy_3+2))';
                L23(l2(1)*3,:)=L23(l2(1)*3,:)+1/16* ( c23.*(sx_3.^3-3*sx_3+2).*(sy_3.^3-sy_3.^2-sy_3+1 ))';
                L23(l2(2)*3-2,:)= L23(l2(2)*3-2,:)-1/16* ( c23.*(sx_3.^3-3*sx_3-2).*(sy_3.^3-3*sy_3+2) )';
                L23(l2(2)*3-1,:)=L23(l2(2)*3-1,:)+1/16*(c23.*(sx_3.^3-3*sx_3.^2-sx_3-1).*(sy_3.^3-3*sy_3+2))';
                L23(l2(2)*3,:)=L23(l2(2)*3,:)-1/16*(c23.*((sx_3.^3-3*sx_3-2).*(sy_3.^3-sy_3.^2-sy_3+1 )))';
                L23(l2(3)*3-2,:)=L23(l2(3)*3-2,:)+1/16*(c23.*(sx_3.^3-3*sx_3-2).*(sy_3.^3-3*sy_3-2))';
                L23(l2(3)*3-1,:)=L23(l2(3)*3-1,:)-1/16*(c23.*(sx_3.^3-3*sx_3.^2-sx_3-1).*(sy_3.^3-3*sy_3-2))';
                L23(l2(3)*3,:)=L23(l2(3)*3,:)-1/16*(c23.*(sx_3.^3-3*sx_3-2).*(sy_3.^3-3*sy_3.^2-sy_3-1))';
                L23(l2(4)*3-2,:)=L23(l2(4)*3-2,:)-1/16*(c23.*(sx_3.^3-3*sx_3+2).*(sy_3.^3-3*sy_3-2))';
                L23(l2(4)*3-1,:)=L23(l2(4)*3-1,:)-1/16*(c23.*(sx_3.^3-sx_3.^2-sx_3+1 ).*(sy_3.^3-3*sy_3-2))';
                L23(l2(4)*3,:)= L23(l2(4)*3,:)+1/16*(c23.*(sx_3.^3-3*sx_3+2).*(sy_3.^3-3*sy_3.^2-sy_3-1))';
                %L4
                L24(l2(1)*3-2,:)= L24(l2(1)*3-2,:)+1/16*( c24.*(sx_4.^3-3*sx_4+2) .*(sy_4.^3-3*sy_4+2))';
                L24(l2(1)*3-1,:)=L24(l2(1)*3-1,:)+1/16* ( c24.*( sx_4.^3-sx_4.^2-sx_4+1).*(sy_4.^3-3*sy_4+2))';
                L24(l2(1)*3,:)=L24(l2(1)*3,:)+1/16* ( c24.*(sx_4.^3-3*sx_4+2).*(sy_4.^3-sy_4.^2-sy_4+1 ))';
                L24(l2(2)*3-2,:)= L24(l2(2)*3-2,:)-1/16* ( c24.*(sx_4.^3-3*sx_4-2).*(sy_4.^3-3*sy_4+2) )';
                L24(l2(2)*3-1,:)=L24(l2(2)*3-1,:)+1/16*(c24.*(sx_4.^3-3*sx_4.^2-sx_4-1).*(sy_4.^3-3*sy_4+2))';
                L24(l2(2)*3,:)=L24(l2(2)*3,:)-1/16*(c24.*((sx_4.^3-3*sx_4-2).*(sy_4.^3-sy_4.^2-sy_4+1 )))';
                L24(l2(3)*3-2,:)=L24(l2(3)*3-2,:)+1/16*(c24.*(sx_4.^3-3*sx_4-2).*(sy_4.^3-3*sy_4-2))';
                L24(l2(3)*3-1,:)=L24(l2(3)*3-1,:)-1/16*(c24.*(sx_4.^3-3*sx_4.^2-sx_4-1).*(sy_4.^3-3*sy_4-2))';
                L24(l2(3)*3,:)=L24(l2(3)*3,:)-1/16*(c24.*(sx_4.^3-3*sx_4-2).*(sy_4.^3-3*sy_4.^2-sy_4-1))';
                L24(l2(4)*3-2,:)=L24(l2(4)*3-2,:)-1/16*(c24.*(sx_4.^3-3*sx_4+2).*(sy_4.^3-3*sy_4-2))';
                L24(l2(4)*3-1,:)=L24(l2(4)*3-1,:)-1/16*(c24.*(sx_4.^3-sx_4.^2-sx_1+1 ).*(sy_4.^3-3*sy_4-2))';
                L24(l2(4)*3,:)= L24(l2(4)*3,:)+1/16*(c24.*(sx_4.^3-3*sx_4+2).*(sy_4.^3-3*sy_4.^2-sy_4-1))';
            end
        end
        Fb(:,tt2)  = (L21*P(1,tt2)+L22*P(2,tt2)+L23*P(3,tt2)+L24*P(4,tt2));
    end
    
    
    Fb(:,1)=0;
    Fb(:,end)=0;
     % Boundary conditions
%      y0=d2/2;
    
%      x_pos= x0+v*t2;                        % position of force wrt to x-axis
%      y_pos= y0*(ones(1,length(t2)));  % position of force wrt to y-axis
% Set Boundary conditions
bc1= find(NL(:,1)==0); % fixed free
bc2=find(NL(:,1)==d1);

% bc1= find(NL(:,1)==0|NL(:,1)==d1|NL(:,2)==0|NL(:,2)==d2); % fixed fixed on all edges
bc_Dof=[];
bc_dof=[];
nn1= length(bc1); % how many nodes have x=0 or d1 as entry

for i= 1:nn1
    bc_Dof=[bc_Dof; 3*bc1(i)-2; 3*bc1(i)-1;3*bc1(i)];
    bc_dof=[bc_dof;3*bc2(i)-2];
end
bb1=union(bc_Dof,bc_dof);
bb2= union(bc1,bc2);
bc_00=unique(bc1);
% bc_Dof_1= unique(bc_Dof);
bc_Dof_1= unique(bb1);
% LL= size(bc_00,1);
L_bc= length(bc_Dof_1);

%%%%%%%%%%%%
BC_Value= zeros(1,L_bc);
%%%%%%%%%%%%%%%%%inpose constraints%%%%%%%
iddx= 1:DoF;
bc_node_1= 3*iddx-2;

iddx(bc_Dof_1)=[];
nn=length(bc_Dof_1);
bca = repmat(zeros(nn,1),1,n);
bcv = repmat(zeros(nn,1),1,n);
bcy = repmat(zeros(nn,1),1,n);
% y(bc_Dof_1,:)=BC_Value'*ones(1,n); 

p0=Fb(iddx,:)-Mb(iddx,bc_Dof_1)*bca-Cb(iddx,bc_Dof_1)*bcv-Kb(iddx,bc_Dof_1)*bcy;
Mb0= Mb(iddx,iddx);
Cb0= Cb(iddx,iddx);
Kb0=Kb(iddx,iddx);
[ t, y0, dy0, ddy0 ] = newmark_beta( Mb0,Cb0,Kb0,p0,1/1000,zeros(DoF-nn,3), Gama, Beta );
y(iddx,:)=y0;
dy(iddx,:)=dy0;
ddy(iddx,:)=ddy0;
 end
     y(:,end)=0;
     dy(:,end)=0;
     ddy(:,end)=0;
     %-- visualization
% fprintf('outputting the simulation results...\n')
% figure(1)
% subplot(3,1,1); plot(t,u(1:4,:),'-r');hold on
% subplot(3,1,1); plot(t,rx(1:4,:),'-b');hold off
% subplot(3,1,2); plot(t,z(1:4,:),'-r')
% subplot(3,1,3); plot(t,z(5:8,:),'-r')
     
     
     
     
%     
%      x_pos= x0+v*t2;                        % position of force wrt to x-axis
%      y_pos= y0*(ones(1,length(t2)));  % position of force wrt to y-axis
% % Set Boundary conditions
% %bc1= find(NL(:,1)==0|NL(:,1)==d1); % fixed free
% 
% bc1= find(NL(:,1)==0|NL(:,1)==d1|NL(:,2)==0|NL(:,2)==d2); % fixed fixed on all edges
% bc_Dof=[];
% nn1= length(bc1); % how many nodes have x=0 or d1 as entry
% 
% for i= 1:nn1
%     bc_Dof=[bc_Dof; 3*bc1(i)-2; 3*bc1(i)-1;3*bc1(i)];
%     
% end
% bc_00=unique(bc1);
% bc_Dof_1= unique(bc_Dof);
% LL= size(bc_00,1);
% L_bc= length(bc_Dof_1);
% 
% %%%%%%%%%%%%
% BC_Value= zeros(1,L_bc);
% %%%%%%%%%%%%%%%%%inpose constraints%%%%%%%
% iddx= 1:DoF;
% bc_node_1= 3*iddx-2;
% 
% iddx(bc_Dof_1)=[];
% y(bc_Dof_1,:)=BC_Value'*ones(1,n); 
% %-- load
% xp= find(NL(:,1)==x_pos);
% % yp= find(NL(:,1)==0.5&NL(:,2)==0.5);
% yp= find(NL(:,1)==0.5);
% n_yp= length(yp);
% idd=1:DoF;
% n_yp_dof=(3*yp(1):3:3*yp(end))-2 ;
% 
% %##########ASSIGN FORCE###########
% 
% ff(n_yp_dof,:)=-40;      % Uniformly distributed load along y=0.5 line
% 
% y0= y(iddx,:);
% dy0= y(iddx,:);
% ddy0= y(iddx,:);
% p0=ff(iddx,:)-M(iddx,bc_Dof_1)*ddy(bc_Dof_1,:)-Cb(iddx,bc_Dof_1)*dy(bc_Dof_1,:)-K(iddx,bc_Dof_1)*y(bc_Dof_1,:);
% M0= M(iddx,iddx);
% C0= Cb(iddx,iddx);
% K0=K(iddx,iddx);

% Newmarkbeta
% G=1/2; B= 1/4;
% [y0,dy0,ddy0]=newmark_beta(Mb0,Cb0,Kb0,dt,G,B,y0,dy0,ddy0,p0);
% y(iddx,:)=y0;
% dy(iddx,:)=dy0;
% ddy(iddx,:)=ddy0;
% 
%visualization
disp('Visualization');
id=1:3:DoF;
deflect= y(1:3:DoF,:);
Theta_x= y(2:3:DoF,:);
Theta_y= y(3:3:DoF,:);
acc= ddy(id,:);
velo= dy(id,:);
% nm=find(NL(:,1)==d1/2&NL(:,2)==d2/2);
nm=round(NoN/2);
mag=100;
% z=zeros(1,121);
 figure(1)
plot(t2,velo(nm,:),'-r'); legend('Velocity')
figure(2)
plot(t2,acc(nm,:),'-r'); legend('acceleration')
figure(3)
plot(t2,deflect(nm,:),'-r'); legend('Dispacement')
figure(4);set(gcf, 'Unit', 'pixel', 'Position',[120, 280, 800, 300]);
writerObj = VideoWriter('Plate_1.avi','Motion JPEG AVI'); %-- setting movie class
writerObj.FrameRate = 20; %-- flame rate (slide/sec)
open(writerObj); %-- initialization of movie object
for tt=1:10:n
%     plot3(NL(:,1),NL(:,2),zeros(length(NL),1),'k.')
%     hold on
    
    for ii=1:m+1
            xlines=(ii-1)*(p+1)+(1:p+1);
            for j=1:p+1
            ylines=j:p+1:NoN;
            
            plot3(NL(xlines,1),NL(xlines,2),mag*deflect(xlines,tt),'b-')
            hold on
            plot3(NL(ylines,1),NL(ylines,2),mag*deflect(ylines,tt),'b-')
            end
    end
    plot3(NL(:,1),NL(:,2),mag*deflect(:,tt),'r.')
    text(0.2, 0.5, 6.0, ['TIME=' num2str(t2(tt)) ' [s]'])
    hold off
    set(gcf,'Position',[120,141,733,613]);
%     set(gca,'CameraPosition',[-140 -57 32])
%     set(gca,'CameraTarget',[50 4 -0.5])
    set(gca,'CameraPosition',[-19 -24 27])
    set(gca,'CameraTarget',[2.5 2 -0.5])
    set(gca,'zlim',[-5 4])
    frame = getframe(figure(4)); %-- get the figure as a frame of the animation
    writeVideo(writerObj,frame); %-- add the gotten frame to the movie object
end
close(writerObj); %-- end of animation
toc




