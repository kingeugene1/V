 tic
clear;  close;    clc;
d1= 5;                                                          % x-length [m]
d2= 1;                                                          % y- length [m]
p= 10;                                                          % number of elements in x-direction
m=10;                                                          % number of elements in y-direction
da= d1/p;                                                     % increment size in x-direction
db= d2/m;                                                 % increment size in y-direction
NoE= p*m;                                                   % Total number of elements
NoN= (p+1)*(m+1);                                     % Total number of nodes
NPE= 4;                                                       % Number of elements per node
PD= 2;                                                         % Problem dimension
N_Dof = 3;                                                   % number of degree of freedom per node
DoF = N_Dof*NoN;                                      % Total degree of freedom of the system
Dof_E= N_Dof*NPE;                                    % degree of freedom per element
E= 2.5e11;                                                   % Plate young's modulus
% E= 16.5e11;
T= 0.1;                                                        % Plate thickness
nu=0.3;
D0= E*T^3/(12*(1-nu^2));                           % Flexural rigidity of plate
rhA=3000;
Gama=1/2; Beta=1/4;
a=0.2300; b= 0.00500;
t1=4;
dt= 0.001;
t2 = 0:dt:t1;
n= length(t2);
velo_p= 10; %[m/s]
R = load('profile-1.txt'); 
x0=5; xb0=10;
x_1 = (0:0.10:(length(R)-1)*0.10)-xb0;
x_t= (x0+velo_p*t2)-xb0;  %[m]
y_t=0.55*ones(1,n);
rx= interp1(x_1,R,x_t,'spline');
%%%%%**************************************************************************
[NL,EL]=uniform_mesh(d1,d2,p,m);             % this function return the node list as coordinates (x,y) and element list

[D]=Dmatrix(E,T,nu);                                   % this function return D matrix of isoparametric material
% D=[D0 0 0;0 D0 0;0 0 D0];
Kb= zeros(DoF,DoF);                                     % initialization of stiffness, Force and Mass matrices
Force= zeros(DoF,n);
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
   f= zeros(Dof_E,1);
   pp0=-1;
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
          fe= Force_vector(NPE,N,pp0);
      
%           [A,B]= AB_Matrix(xi,eta,xx,yy);
          [A,G,B]= ke_Matrix(xi,eta,xx,yy);
          ke=ke+det_J*wx*wy*(G.')*(A.')*D*B*G;
          f=f+fe*wx*wy*det_J;
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
 ff=zeros(DoF,n);

% moving load parameters

% x_value= 0:da:d1;
% y_value= 0:db:d2;
% % y_pos= 0.5;
% y_t= (linspace(0,d2,n)); % rand values [0,db]
p1=[-4;0;0];

for counter=1:6

yt = zeros(1,length(t2));
     for tt=1:length(t2)
         L11 = zeros( DoF, 1);
         for jl=1:p
             for il=1:m
                 %%%%
                 l2=EL(((il-1)*p)+jl,:);
                 x_c_2=NL(l2,1)/da;
                 y_c_2=NL(l2,2)/db;
                 z_c_2=[x_c_2 y_c_2];
                 b12=z_c_2(2,:)-z_c_2(1,:);
                 b22=z_c_2(3,:)-z_c_2(2,:);
                 b32=z_c_2(4,:)-z_c_2(3,:);
                 b42=z_c_2(1,:)-z_c_2(4,:);
                 B_c2=[b12;b22;b32;b42];
                 %                 B_c2= repmat(B_c2,4,1);
                 sx_1=(x_t(1,tt))/da;
                 c1 = x_t(1,tt) > da*(jl-1) & x_t(1,tt) <= da*jl;
                 sy_1=y_t(1,tt)/db;
                 
                 k_c1=[sx_1 sy_1];
            
                 
                 %first point
                 a1=k_c1-z_c_2(1,:);
                 a2=k_c1-z_c_2(2,:);
                 a3=k_c1-z_c_2(3,:);
                 a4=k_c1-z_c_2(4,:);
                 
                 A=[a1;a2;a3;a4];
                 
                 B=B_c2.*A;
                 
                 s_21=B(B~=0);
                 
                 c21=all(s_21>0);
                
                 %L1
                 L11(l2(1)*3-2,:)= L11(l2(1)*3-2,:)+1/16*(c21.*((sx_1-1).^2*(sx_1+2)) .*((sy_1-1).^2*(sy_1+2)))';
                 L11(l2(1)*3-1,:)=L11(l2(1)*3-1,:)+1/16* (c21.*( (3*sx_1^2-2*sx_1+3)).*((sy_1-1).^2*(sy_1+2)))';
                 L11(l2(1)*3,:)=L11(l2(1)*3,:)+1/16* ( c21.*((sx_1-1).^2*(sx_1+2)).*((3*sy_1^2-2*sy_1+3)))';
                 
                 L11(l2(2)*3-2,:)= L11(l2(2)*3-2,:)+1/16* ( c21.*((sx_1+1).^2*(sx_1-2)).*((sy_1-1).^2*(sy_1+2)) )';
                 L11(l2(2)*3-1,:)=L11(l2(2)*3-1,:)+1/16*(c21.*((3*sx_1^2+2*sx_1-1)).*((sy_1-1).^2*(sy_1+2)))';
                 L11(l2(2)*3,:)=L11(l2(2)*3,:)+1/16*(c21.*((sx_1+1).^2*(sx_1-2)).*((3*sy_1^2-2*sy_1+3)))';
                
                 L11(l2(3)*3-2,:)=L11(l2(3)*3-2,:)+1/16*(c21.*((sx_1+1).^2*(sx_1-2)).*((sy_1+1).^2*(sy_1-2)))';
                 L11(l2(3)*3-1,:)=L11(l2(3)*3-1,:)+1/16*(c21.*((3*sx_1^2+2*sx_1-1)).*((sy_1+1).^2*(sy_1-2)))';
                 L11(l2(3)*3,:)=L11(l2(3)*3,:)+1/16*(c21.*((sx_1+1).^2*(sx_1-2)).*((3*sy_1^2+2*sy_1-1)))';
          
                 L11(l2(4)*3-2,:)=L11(l2(4)*3-2,:)+1/16*(c21.*((sx_1-1).^2*(sx_1+2)).*((sy_1+1).^2*(sy_1-2)))';
                 L11(l2(4)*3-1,:)=L11(l2(4)*3-1,:)+1/16*(c21.*((3*sx_1^2-2*sx_1+3) ).*((sy_1+1).^2*(sy_1-2)))';
                 L11(l2(4)*3,:)= L11(l2(4)*3,:)+1/16*(c21.*((sx_1-1).^2*(sx_1+2)).*((3*sy_1^2+2*sy_1-1)))';
                 
             end
         end
         yt(:,tt)  = L11'*y(:,tt);
     end
     
 
    u = rx+yt;
NN= zeros(DoF,3);
    Fb= zeros(DoF,n);
    Fb1= zeros(DoF,n);
    for tt2=1:length(t2)
         L21 = zeros( DoF, 1);
         for j2=1:p
             for i2=1:m
                 l2=EL(((i2-1)*p)+j2,:);
                 x_c_2=NL(l2,1)/da;
                 y_c_2=NL(l2,2)/db;
                 z_c_2=[x_c_2 y_c_2];
                 b12=z_c_2(2,:)-z_c_2(1,:);
                 
                 B_c2=[b12;b22;b32;b42];
                 %                 B_c2= repmat(B_c2,4,1);
                 sx_1=(x_t(1,tt2))/da;
                 c2 = x_t(1,tt2) > da*(j2-1) & x_t(1,tt2) <= da*j2;
                 sy_1=y_t(1,tt2)/db;
                 
                 k_c1=[sx_1 sy_1];
                 
                 %first point
                 a12_1=k_c1-z_c_2(1,:);
                 a22_1=k_c1-z_c_2(2,:);
                 a32_1=k_c1-z_c_2(3,:);
                 a42_1=k_c1-z_c_2(4,:);
                 
                 
                 A_c_21=[a12_1;a22_1;a32_1;a42_1];
                 
                 BA_21=B_c2.*A_c_21;
                 
                 s_21=BA_21(BA_21~=0);
                 
                 c21=all(s_21>0);
                 
                 %L1
                 L21(l2(1)*3-2,:)= L21(l2(1)*3-2,:)+1/16*( c21.*((sx_1-1).^2*(sx_1+2)) .*((sy_1-1).^2*(sy_1+2)))';
                 L21(l2(1)*3-1,:)=L21(l2(1)*3-1,:)+1/16* ( c21.*( (3*sx_1^2-2*sx_1+3)).*((sy_1-1).^2*(sy_1+2)))';
                 L21(l2(1)*3,:)=L21(l2(1)*3,:)+1/16* ( c21.*((sx_1-1).^2*(sx_1+2)).*((3*sy_1^2-2*sy_1+3)))';
NN(l2(1)*3-2,1)=L21(l2(1)*3-2);
NN(l2(1)*3-1,2)=L21(l2(1)*3-1);
NN(l2(1)*3,3)=L21(l2(1)*3);
                 L21(l2(2)*3-2,:)= L21(l2(2)*3-2,:)+1/16* ( c21.*((sx_1+1).^2*(sx_1-2)).*((sy_1-1).^2*(sy_1+2)) )';
                 L21(l2(2)*3-1,:)=L21(l2(2)*3-1,:)+1/16*(c21.*((3*sx_1^2+2*sx_1-1)).*((sy_1-1).^2*(sy_1+2)))';
                 L21(l2(2)*3,:)=L21(l2(2)*3,:)+1/16*(c21.*((sx_1+1).^2*(sx_1-2)).*((3*sy_1^2-2*sy_1+3)))';
 NN(l2(2)*3-2,1)=L21(l2(2)*3-2);
NN(l2(2)*3,2)=L21(l2(2)*3-1);
NN(l2(2)*3,3)=L21(l2(2)*3);
                 L21(l2(3)*3-2,:)=L21(l2(3)*3-2,:)+1/16*(c21.*((sx_1+1).^2*(sx_1-2)).*((sy_1+1).^2*(sy_1-2)))';
                 L21(l2(3)*3-1,:)=L21(l2(3)*3-1,:)+1/16*(c21.*((3*sx_1^2+2*sx_1-1)).*((sy_1+1).^2*(sy_1-2)))';
                 L21(l2(3)*3,:)=L21(l2(3)*3,:)+1/16*(c21.*((sx_1+1).^2*(sx_1-2)).*((3*sy_1^2+2*sy_1-1)))';
NN(l2(3)*3-2,1)=L21(l2(3)*3-2);
NN(l2(3)*3,2)=L21(l2(3)*3-1);
NN(l2(3)*3,3)=L21(l2(3)*3);
                 L21(l2(4)*3-2,:)=L21(l2(4)*3-2,:)+1/16*(c21.*((sx_1-1).^2*(sx_1+2)).*((sy_1+1).^2*(sy_1-2)))';
                 L21(l2(4)*3-1,:)=L21(l2(4)*3-1,:)+1/16*(c21.*((3*sx_1^2-2*sx_1+3) ).*((sy_1+1).^2*(sy_1-2)))';
                 L21(l2(4)*3,:)= L21(l2(4)*3,:)+1/16*(c21.*((sx_1-1).^2*(sx_1+2)).*((3*sy_1^2+2*sy_1-1)))';
NN(l2(4)*3-2,1)=L21(l2(4)*3-2);
NN(l2(4)*3,2)=L21(l2(4)*3-1);
NN(l2(4)*3,3)=L21(l2(4)*3);
            end
        end
        Fb(:,tt2)  = L21*pp0;
        Fb1(:,tt2)= NN*p1;
    end
    
    
    Fb(:,1)=0;
    Fb(:,end)=0;
   
% Set Boundary conditions
bc1= find(NL(:,1)==0); % fixed free
bc2=find(NL(:,1)==d1);

%  bc1= find(NL(:,1)==0|NL(:,1)==d1|NL(:,2)==0|NL(:,2)==d2); % fixed fixed on all edges
bc_Dof=[];
bc_dof=[];
nn1= length(bc1); % how many nodes have x=0 or d1 as entry

for i= 1:nn1
    bc_Dof=[bc_Dof; 3*bc1(i)-2; 3*bc1(i)-1;3*bc1(i)];
    bc_dof=[bc_dof;3*bc2(i)-2];
end
bb1=union(bc_Dof,bc_dof);
% bb1=union(bc_Dof);
bb2= union(bc1,bc2);
bc_00=unique(bc1);
% bc_Dof_1= unique(bc_Dof);
bc_Dof_1= unique(bb1);
% LL= size(bc_00,1);
L_bc= length(bc_Dof_1);

iddx= 1:DoF;
bc_node_1= 3*iddx-2;

iddx(bc_Dof_1)=[];
nn=length(bc_Dof_1);
bca = repmat(zeros(nn,1),1,n);
bcv = repmat(zeros(nn,1),1,n);
bcy = repmat(zeros(nn,1),1,n);

% f0=Fb(iddx,:)-Mb(iddx,bc_Dof_1)*bca-Cb(iddx,bc_Dof_1)*bcv-Kb(iddx,bc_Dof_1)*bcy;
f0=Fb1(iddx,:)-Mb(iddx,bc_Dof_1)*bca-Cb(iddx,bc_Dof_1)*bcv-Kb(iddx,bc_Dof_1)*bcy;
Mb0= Mb(iddx,iddx);
Cb0= Cb(iddx,iddx);
Kb0=Kb(iddx,iddx);
[ t, y0, dy0, ddy0 ] = newmark_beta( Mb0,Cb0,Kb0,f0,1/1000,zeros(DoF-nn,3), Gama, Beta );
y(iddx,:)=y0;
dy(iddx,:)=dy0;
ddy(iddx,:)=ddy0;
%%
end

%visualization
mag=1000;
id=1:3:DoF;
deflect= y(1:3:DoF,:);
Theta_x= y(2:3:DoF,:);
Theta_y= y(3:3:DoF,:);

xxx= deflect(round(NoN/2),:);
figure(1)
 plot(t2,xxx) %deflection at the middle of plate
%  ylim([-0.01 0.01])
 title('Deflection of middle point')
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
            plot3(NL(xlines,1),NL(xlines,2),mag*deflect(xlines,tt),'b-')
            hold on
            plot3(NL(ylines,1),NL(ylines,2),mag*deflect(ylines,tt),'b-')
    end
    plot3(NL(:,1),NL(:,2),mag*deflect(:,tt),'r.')
%     line([x_t(1,tt) x_t(1,tt)],[y_t(1,tt) y_t(1,tt)],[0 3],'Color','black','LineWidth',3)
%      line(x_t(1,tt),y_t(1,tt),'Color','k','marker','V','MarkerSize',3,'linewidth',3)
    text(6.0, 0.5, 3.5, ['TIME=' num2str(t2(tt)) ' [s]'])
    hold off
    set(gcf,'Position',[120,141,733,613]);
    set(gca,'CameraPosition',[-26 -7 25])
    set(gca,'CameraTarget',[0.5 0.5 -3])
    set(gca,'zlim',[-6 5])
%     set(gca,'xlim',[0 15])
    frame = getframe(figure(2)); %-- get the figure as a frame of the animation
    writeVideo(writerObj,frame); %-- add the gotten frame to the movie object
end
close(writerObj); %-- end of animation
toc




