function [ t, x, dx, ddx ] = newmark_beta_correct( m,c,k,f,dt,X, G, B)
%% *NEWMARK BETA METHOD* 
% returns the results of dynamic simulation
% 
%% output
% x: the displacement vibration response
% dx: the velocity vibration response
% ddx: the acceleration vibration response
% 
%% input
% m: mass matrix R(N,N)
% c: damping matrix R(N,N)
% k: stiffness matrix R(N,N)
% f: external force R(N,T)
% dt: time interval
% X: initial values of x(1), dx(1) and ddx(1) R(N,3)
% G: gamma of newmark beta method (usually 1/2)
% B: beta of newmark beta method (usually 1/6)
% 
% N: degree of freedom
% T: length of data
% 
%% information
% This function is coded by Kyosuke Yamamoto, used for the course of 
% "Advanced Reliability Engineering", GSSIE, Univ. of Tsukuba, 2018.
% 
% coded on MAY 2, 2018
% 

%% initialization
%-- memory space
x=zeros(size(f));
dx=zeros(size(f));
ddx=zeros(size(f));

%-- input initial values
x(:,1)=X(:,1);
dx(:,1)=X(:,2);
ddx(:,1)=X(:,3);

%-- time vector
T=length(f);
t=(0:T-1)*dt;

%% solving

A = m + dt.*G.*c + dt^2.*B.*k; %-- global matrix
M = inv(A); %-- inverse matrix of A

for tt=2:T
    %-- solving preparation
    b1 = -c * ( dx(:,tt-1) + (1-G)*dt*ddx(:, tt-1) );
    b2 = -k * ( x(:,tt-1) + dt*dx(:,tt-1) + (1/2-B)*dt^2*ddx(:,tt-1) );
    b = f(:,tt) + b1 + b2; %-- right hand side
    
    %-- calculation of acceleration at t=t+1
    ddx(:,tt) = M*b;
    
    %-- calculation of velocity and displacement responses
    dx(:,tt) = dx(:,tt-1) + dt*(1-G)*ddx(:,tt-1) + dt*G*ddx(:,tt);
    x(:,tt) = x(:,tt-1) + dt*dx(:,tt-1) + dt^2*(1/2-B)*ddx(:,tt-1) + dt^2*B*ddx(:,tt);
    
end


end