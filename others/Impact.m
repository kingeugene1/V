clear; close all; clc;
dyna= load('dynamic.txt');
stat= load('static.txt');
t= 0:0.001:4;
plot(t,dyna);
hold on
plot(t,stat);
hold off
dn=dyna(dyna~=0);
sta= stat(stat~=0);
I=(max(abs(dn))-max(abs(sta)))/max(abs(sta));
fprintf([' Impact factor is=' num2str(I) '\n'])