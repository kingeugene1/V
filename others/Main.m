tic
clear;  close all;    clc;
d1= 20;                                                          % x-length [m]
d2= 4;                                                          % y- length [m]
p= 10;                                                          % number of elements in x-direction
m=10;                                                          % number of elements in y-direction
da= d1/p;                                                     % increment size in x-direction
db=d2/m;
NoE= p*m;                                                   % Total number of elements
NoN= (p+1)*(m+1);                                     % Total number of nodes
NPE= 4;                                                       % Number of elements per node
PD= 2;                                                         % Problem dimension
N_Dof = 3;                                                   % number of degree of freedom per node
DoF = N_Dof*NoN;                                      % Total degree of freedom of the system
Dof_E= N_Dof*NPE;                                    % degree of freedom per element

E=5e10;       
rhA=3000;
thick= d1/20;                                                        % Plate thickness
rhA_data  = rhA* ones(NoE, 1);
E_data   = E * ones(NoE, 1); % set elements around the target point to have a reduced stiffness for damaged case
% E_DATA([145 146 147 148 149 150 151 152 153  178 179 180 181 182 183 184 185 186 211 212 213 214 215 216 217 218 219])= E*0.2;

nu=0.33;


g=9.81;                                           %[m/s^2]
%Rayleigh Damping coefficients************************************

a=0.7021; b=0.0052;
TRH= 1e-6; %-- Maximum error
%%  vehicle Parameters
fprintf('setting vehicle parameters...\n')
v   = 10;                                                               %-- speed (m/s)
L= 2.4; 
l=L/2;
D   = 4.400;      %-- [m] distance between the front-rear axles
dd1  = 1.215;      %-- [m] distance to the front axle from G
dd2  = D-dd1;       %-- [m] distance to the rear axle from G                                                    %[m]
cs1 = 29e3;  %-- sprung-damping (kg/s)
cs2 = 24.2e3;  %-- sprung-damping (kg/s)
cs3 = 29e3;  %-- sprung-damping (kg/s)
cs4 = 24.2e3;  %-- sprung-damping (kg/s)

ks1= 456e3;  %-- sprung-stifness (N/m)=(kg/s/s)
ks2= 41e4;  %-- sprung-stifness (N/m)=(kg/s/s)
ks3= 456e3;  %-- sprung-stifness (N/m)=(kg/s/s)
ks4= 41e4;  %-- sprung-stifness (N/m)=(kg/s/s)

ku1 = 431e4; %-- unsprung-stifness (N/m)=(kg/s/s)
ku2 = 431e4; %-- unsprung-stifness (N/m)=(kg/s/s)
ku3 = 431e4; %-- unsprung-stifness (N/m)=(kg/s/s)
ku4 = 431e4; %-- unsprung-stifness (N/m)=(kg/s/s)
ms  = 13560;       %-- [kg] mass of the vehicle
mu1      =      751;                                                 %-- unsprung-mass (kg)
mu2      =      469;                                                 %-- unsprung-mass (kg)
mu3      =      751;
mu4      =      469;
ms1=(l*dd2)/((dd1+dd2 )*(2*l))*ms; ms2=(l*dd1)/((dd1+dd2 )*(2*l))*ms ; ms3= (l*dd2)/((dd1+dd2 )*(2*l))*ms ; ms4=(l*dd1)/((dd1+dd2 )*(2*l))*ms;
dt  = 0.001;                                                             %-- time increment
x0  =10;                                                                  %-- start position of G(m)
xb0= 20;    % Entrance position
t1=5;

%-- calculate input
fprintf('calculating vehicle input...\n')
% [k,~]=mk_road_profile_01();
R = load('profile-1.txt');  %-- load: R(x)
% R=k(:,2);
x_1 = (0:0.10:(length(R)-1)*0.10)-xb0;
t2 = 0:dt:t1;
n= length(t2);

y_y= [ 0.68*ones(1,n); L+0.68*ones(1,n); 0.68*ones(1,n);L+0.68*ones(1,n)]; 

x_x= [x0+v*t2+dd1 ;x0+v*t2+dd1; x0+v*t2-dd2 ;x0+v*t2-dd2]-xb0;    %-- vehicle position: x(t)

rx = spline(x_1,R',x_x);                                                     %-- spline interpolation: R(x) to r(t)=R(x(t))
% rx = pchip(x_1,R',x_x);
%-- simulate output
fprintf('vehicle vibration simulation: newmark-beta method...\n')
Gama = 1/2;                                                                        %-- gammma coefficient for newmark-beta method
Beta = 1/4;                                                                         %-- beta coefficient for newmark-beta method                  

[Mv,Cv,Kv,Mp]= Veh_sys(L,dd1,dd2);
%%%%%*******************Bridge System*************************************************
fprintf('bridge matrices...\n')
[NL,EL]=uniform_mesh(d1,d2,p,m);             % this function return the node list as coordinates (x,y) and element list
[D]=Dmatrix(E_data,NoE,thick,nu);                                   % this function return D matrix of isoparametric material
[Kb,Mb]= M_K_mat(NoE,NPE,PD,NL,EL,D,rhA_data);
 Cb= a*Mb+ b*Kb;    % Rayleigh damping
 
 %%
 ddy=zeros(DoF,n);
 dy=zeros(DoF,n);
 y=zeros(DoF,n);
 mem= zeros(DoF,n); %-- memory
 
 fprintf('VBI-system simulation: newmark-beta method...\n')
 for count=1:100
     
     yt = zeros(4,length(t2));
     for tt=1:length(t2)
          [L]= L_distribution(DoF,p,m,EL,NL,da,db,tt,x_x,y_y);
          yt(:,tt)  = L'*y(:,tt);
     end
   
     u = rx+yt;
     Fv= [zeros(4,n); [ku1 0 0 0; 0 ku2 0 0; 0 0 ku3 0;0 0 0 ku4]*u];
     [ t2, z, dz, ddz ] = newmark_beta( Mv,Cv,Kv,Fv,1/1000,zeros(8,3), Gama, Beta );
     %-- bridge system
%      P = Mp*(repmat(g,8,n)-ddz);
     P1 = -ms1*(g+ddz(1,:)) - mu1*(g+ddz(5,:)); %-- [N] at the front axle
    P2 = -ms2*(g+ddz(2,:)) - mu2*(g+ddz(6,:)); %-- [N] at the rear axle
    P3 = -ms3*(g+ddz(3,:)) - mu3*(g+ddz(7,:)); %-- [N] at the front axle
    P4 = -ms4*(g+ddz(4,:)) - mu4*(g+ddz(8,:)); %-- [N] at the rear axle
    P  = [P1;P2;P3;P4];
     Fb= zeros(DoF,n);
     
     for tt2=1:length(t2)
         [L2]= L_distribution(DoF,p,m,EL,NL,da,db,tt2,x_x,y_y);
          Pt=[P(1,tt2);P(2,tt2);P(3,tt2);P(4,tt2)];
           Fb(:,tt2)  = L2*Pt;
    
     end
% 
     Fb(:,1)=0;
     Fb(:,end)=0;
     
     % Set Boundary conditions
     bc1= find(NL(:,1)==0); % fixed fixed
     bc2=find(NL(:,1)==d1);
    nn1= length(bc1); % how many nodes have x=0 or d1 as entry
     bc_Dof=zeros(3*nn1,1);
     bc_dof=zeros(nn1,1);
     for i= 1:nn1
         bc_Dof(3*i-2:3*i,:)=[3*bc1(i)-2; 3*bc1(i)-1;3*bc1(i)];
         bc_dof(3*i-2:3*i,:)=[3*bc2(i)-2;3*bc2(i)-1;3*bc2(i)];
     end
     bb1=union(bc_Dof,bc_dof);
     bc_Dof_1= unique(bb1);
     iddx= 1:DoF;
     iddx(bc_Dof_1)=[];
     nn=length(bc_Dof_1);
     bca = (zeros(nn,1));
     bcv = (zeros(nn,1));
     bcy = (zeros(nn,1));
     p0=Fb(iddx,:)-Mb(iddx,bc_Dof_1)*bca-Cb(iddx,bc_Dof_1)*bcv-Kb(iddx,bc_Dof_1)*bcy;
     Mb0= Mb(iddx,iddx);
     Cb0= Cb(iddx,iddx);
     Kb0=Kb(iddx,iddx);
     [ t2, y0, dy0, ddy0 ] = newmark_beta( Mb0,Cb0,Kb0,p0,1/1000,zeros(DoF-nn,3), Gama, Beta );
     y(iddx,:)=y0;
     dy(iddx,:)=dy0;
     ddy(iddx,:)=ddy0;
     ddy_mem= norm(ddy-mem);
     norm_mem=norm(mem);
     % End itteration condition
     ERR=norm(ddy-mem)/norm(mem);
     if ERR<TRH
         break;
     end
     fprintf([num2str(count,'%0.3d') ': ERR=' num2str(ERR) '\n'])
%      fprintf([num2str(count,'%0.3d') ': norm_ddy-mem=' num2str(ddy_mem) '\n'])
%      fprintf([num2str(count,'%0.3d') ': norm-mem=' num2str(norm_mem) '\n'])
     mem = ddy;
 end
  %%
  % Visualization
[disp,vel,ac]=ploting(t2,u,rx,y,z,ddy,dy,NoN,NL,x_x,m,p);
% Create a table
CreateTable(disp,vel,ac);
toc




