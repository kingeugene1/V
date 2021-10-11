function CreateTable(d,v,a)
t1=table(d,'VariableNames',{'d'});
writetable(t1, 'Displacement-E(20%).txt');
t2=table(v,'VariableNames',{'v'});
writetable(t2, 'Velocity-E(20%).txt');
t3=table(a,'VariableNames',{'a'});
writetable(t3, 'Acceleration-E(20%).txt')
end