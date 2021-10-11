clear 
close all
clc
f=load('Displacement-intact.txt');
t=linspace(-20,100,length(f));
dt=t(2)-t(1);
n= length(f);
fhat=fft(f,n);
psd=fhat.*conj(fhat)/n;
freq=1/(dt*n)*(0:n);
l=1:floor(n/2);
plot(freq(l),psd(l))
indices=(psd>8e-7);
psdclean= psd.*indices;
fhat=indices.*fhat;
fftilt= ifft(fhat);
plot(t,fftilt)
% xlim([-20 50])
% ylim([-0.01 0.01])
y=table(fftilt,'VariableNames',{'fftilt'});
writetable(y, 'accelerationIntact_refined.txt')