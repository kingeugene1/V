tic;
clear;  close;    clc;
d1= 20;                                                          % x-length [m]
d2= 4;                                                          % y- length [m]
p= 10;                                                          % number of elements in x-direction
m=10;                                                          % number of elements in y-direction
dx= d1/p;                                                     % increment size in x-direction
% dy= d2/m;                                                 % increment size in y-direction
NoE= p*m;                                                   % Total number of elements
NoN= (p+1)*(m+1);                                     % Total number of nodes
NPE= 4;                                                       % Number of elements per node
PD= 2;                                                         % Problem dimension
N_Dof = 3;                                                   % number of degree of freedom per node
DoF = N_Dof*NoN;                                      % Total degree of freedom of the system
Dof_E= N_Dof*NPE;                                    % degree of freedom per element
% E= 10920;                                                   % Plate young's modulus
E= 5e10;
T= 1;                                                        % Plate thickness
nu=0.3;
D0= E*T^3/(12*(1-nu^2));                           % Flexural rigidity of plate
rhA=3000;
%variables for newmark************************************
a=0.7300; b= 0.005200;
t1=4;
dt= 0.001;
t2 = 0:dt:t1;
n= length(t2);
%%%%%**************************************************************************
[NL,EL]=uniform_mesh(d1,d2,p,m);             % this function return the node list as coordinates (x,y) and element list
% [D]=Dmatrix(E,T,nu);                                   % this function return D matrix of isoparametric material
D=E*T^3/(12*(1-nu^2));
K= zeros(DoF,DoF);                                     % initialization of stiffness, Force and Mass matrices
Force= zeros(DoF,n);
M= zeros(DoF,DoF);

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
   f= zeros(Dof_E,1);
   p0=-0.01;        %[N]
% p0=-3.7019e8;
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
          fe= Force_vector(NPE,N,p0);
      
          [A,B]= AB_Matrix(xi,eta,xx,yy);
          ke=ke+det_J*wx*wy*A'*B'*D*B*A;
          f=f+fe*wx*wy*det_J;
          me=me+ rhA*T*m1*wx*wy*det_J;
      end
      
  end
  
  % ASSEMBLE STIFFNESS AND FORCE MATRIX
     idx= elementdof(index,NPE,N_Dof);                            % extract dofs for a given element as a list
     
     [K,M]= Assemble(K,M,ke,idx,me);
     
 end
 Cb= a*M+ b*K;    % Rayleigh damping
 
 %%
 
 ddy=zeros(DoF,n);
 dy=zeros(DoF,n);
 y=zeros(DoF,n);
 ff=zeros(DoF,n);
 x0=0;
     v= 10; %[m/s]
%      y0=d2/2;
    
     x_pos= x0+v*t2;                        % position of force wrt to x-axis
%      y_pos= y0*(ones(1,length(t2)));  % position of force wrt to y-axis
% Set Boundary conditions
%bc1= find(NL(:,1)==0|NL(:,1)==d1); % fixed free

bc1= find(NL(:,1)==0|NL(:,1)==d1|NL(:,2)==0|NL(:,2)==d2); % fixed fixed on all edges
bc_Dof=[];
nn1= length(bc1); % how many nodes have x=0 or d1 as entry

for i= 1:nn1
    bc_Dof=[bc_Dof; 3*bc1(i)-2; 3*bc1(i)-1;3*bc1(i)];
    
end
bc_00=unique(bc1);
bc_Dof_1= unique(bc_Dof);
LL= size(bc_00,1);
L_bc= length(bc_Dof_1);

%%%%%%%%%%%%
BC_Value= zeros(1,L_bc);
%%%%%%%%%%%%%%%%%inpose constraints%%%%%%%
iddx= 1:DoF;
bc_node_1= 3*iddx-2;

iddx(bc_Dof_1)=[];
y(bc_Dof_1,:)=BC_Value'*ones(1,n); 
%-- load
xp= find(NL(:,1)==x_pos);
yp= find(NL(:,1)==d1/2&NL(:,2)==d2/2);
% yp= find(NL(:,1)==0.5);
% n_yp= length(yp);
idd=1:DoF;
n_yp_dof=(3*yp(1):3:3*yp(end))-2 ;
% p0=-40;
 ff(n_yp_dof,:)=p0;      % Uniformly distributed load along x=0.5 line
 
% y0= y(iddx,:);
% dy0= y(iddx,:);
% ddy0= y(iddx,:);
f0=ff(iddx,:)-M(iddx,bc_Dof_1)*ddy(bc_Dof_1,:)-Cb(iddx,bc_Dof_1)*dy(bc_Dof_1,:)-K(iddx,bc_Dof_1)*y(bc_Dof_1,:);

M0= M(iddx,iddx);
C0= Cb(iddx,iddx);
K0=K(iddx,iddx);

% Newmarkbeta
G=1/2; B= 1/4;
[~,y0,dy0,ddy0]=newmark_beta(M0,C0,K0,f0,dt,zeros(DoF-L_bc,3),G,B);
y(iddx,:)=y0;
dy(iddx,:)=dy0;
ddy(iddx,:)=ddy0;

%visualization
id=1:3:DoF;
deflect= y(1:3:DoF,:);
Theta_x= y(2:3:DoF,:);
Theta_y= y(3:3:DoF,:);
% z=zeros(1,121);

xxx= deflect(61,:);
ts=table(xxx,'VariableNames',{'xxx'});
writetable(ts, 'static.txt')
figure(1)
 plot(t2,xxx) %deflection at the middle of plate
 title('Deflection of middle point');
 xlabel('Time [sec]','FontSize',12,'Interpreter','Latex');
 ylabel('deflection [mm]','FontSize',12,'Interpreter','Latex');
figure(2);set(gcf, 'Unit', 'pixel', 'Position',[120, 280, 800, 300]);
writerObj = VideoWriter('Plate_1.avi','Motion JPEG AVI'); %-- setting movie class
writerObj.FrameRate = 20; %-- flame rate (slide/sec)
open(writerObj); %-- initialization of movie object
for tt=1:10:n
%     plot3(NL(:,1),NL(:,2),zeros(length(NL),1),'k.')
%     hold on
    
    for ii=1:11
            xlines=(ii-1)*11+(1:11);
            ylines=ii:11:121;
            plot3(NL(xlines,1),NL(xlines,2),deflect(xlines,tt),'b-')
            hold on
            plot3(NL(ylines,1),NL(ylines,2),deflect(ylines,tt),'b-')
    end
    plot3(NL(:,1),NL(:,2),deflect(:,tt),'r.')
    text(6.0, 0.5, 4, ['TIME=' num2str(t2(tt)) ' [s]'])
    hold off
    set(gcf,'Position',[120,141,733,613]);
    set(gca,'CameraPosition',[-41 -26 22])
    set(gca,'CameraTarget',[10 2 -3])
    set(gca,'zlim',[-6 0])
    frame = getframe(figure(2)); %-- get the figure as a frame of the animation
    writeVideo(writerObj,frame); %-- add the gotten frame to the movie object
end
close(writerObj); %-- end of animation
toc;





