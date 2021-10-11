function [Mv,Cv,Kv,Mp]= Veh_sys(L,dd1,dd2)
ms  = 13560;       %-- [kg] mass of the vehicle
mu1      =      751;                                                 %-- unsprung-mass (kg)
mu2      =      469;                                                 %-- unsprung-mass (kg)
mu3      =      751;
mu4      =      469;
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

l=L/2;

Iy =( dd1*dd2*ms)/(dd1+dd2);
Ix =0.5* l*ms;

Mv = [(l*dd2)/((dd1+dd2 )*(2*l))*ms           (l*dd1)/((dd1+dd2 )*(2*l))*ms       (l*dd2)/((dd1+dd2 )*(2*l))*ms       (l*dd1)/((dd1+dd2 )*(2*l))*ms      0      0    0     0;
       -dd2*Ix/(2*l*(dd1+dd2))     -dd1*Ix/(2*l*(dd1+dd2))    dd2*Ix/(2*l*(dd1+dd2))    dd1*Ix/(2*l*(dd1+dd2))  0      0    0     0;
         0.5*Iy/(dd1+dd2)     -0.5*Iy/(dd1+dd2)    0.5*Iy/(dd1+dd2)    -0.5*Iy/(dd1+dd2)    0      0    0     0;
         1              -1             -1            1           0      0    0     0;
       0              0             0             0            mu1    0    0     0;
         0              0             0             0            0      mu2  0     0; 
         0              0             0             0            0      0    mu3   0;
         0              0             0             0            0      0    0     mu4];       


 Cv = [ cs1       cs2      cs3        cs4      -cs1      -cs2     -cs3    -cs4;
     -l*cs1    -l*cs2    l*cs3      l*cs4     l*cs1    l*cs2    -l*cs3    -l*cs4;
     dd1* cs1   -dd2* cs2   dd1* cs3    -dd2*cs4    -dd1*cs1    dd2*cs2   -dd1*cs3    dd2* cs4;
     1              -1             -1            1           0      0    0     0;
    -cs1       0       0         0        cs1      0        0        0;
      0       -cs2      0         0        0       cs2       0        0;
      0        0      -cs3        0        0       0        cs3       0;
      0        0       0        -cs4       0       0        0        cs4]; 


 Kv =[ks1       ks2     ks3       ks4     -ks1     -ks2     -ks3     -ks4;
      -l*ks1    -l*ks2   l*ks3     l*ks4    l*ks1    l*ks2   -l*ks3   -l*ks4;
      dd1*ks1    -dd2*ks2   dd1*ks3    -dd2*ks4   -dd1*ks1    dd2*ks2   -dd1*ks3    dd2*ks4;
      1       -1     -1        1       0       0       0       0;
   -ks1       0      0        0       ks1+ku1   0       0       0;
      0       -ks2     0        0       0       ks2+ku2   0       0;
      0        0     -ks3       0       0       0       ks3+ku3   0;
      0        0      0       -ks4      0       0       0       ks4+ku4];
  

 Mp = [(l*dd2)/((dd1+dd2 )*(2*l))*ms          0                    0                  0                   mu1         0       0       0;
              0                      (l*dd1)/((dd1+dd2 )*(2*l))*ms      0                  0                   0        mu2       0       0;
              0                      0                    (l*dd2)/((dd1+dd2 )*(2*l))*ms        0                0         0        mu3      0;
              0                      0                    0                   (l*dd1)/((dd1+dd2 )*(2*l))*ms      0         0        0     mu4];                    

end