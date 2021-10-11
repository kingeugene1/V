t=0:0.001:5;

d1= load('Displacement-E(0%).txt');
d2= load('Displacement-E(20%).txt');

v1= load('Velocity-E(0%).txt');
v2= load('Velocity-E(20%).txt');

ac1= load('Acceleration-E(0%).txt');
ac2= load('Acceleration-E(20%).txt');

fig=figure (1);
subplot(3,1,1);plot(t,d1,'b','LineWidth',0.5);
hold on 
subplot(3,1,1);plot(t,d2,'r','LineWidth',0.5);
% subplot(3,1,1);plot(t,d3,'k');
% plot(t,d4);
hold off
xlabel('time [sec]','Interpreter','tex','FontSize',8)
ylabel('displacement [m]','Interpreter','tex','FontSize',8)
ylim([-0.0005 9e-5])
legend('E-0%','E-20%');
set(groot, 'DefaultLegendInterpreter', 'tex')
set(groot, 'DefaultLegendFontSize', 8)
 subplot(3,1,2);plot(t,v1,'b','LineWidth',0.5);
hold on 
subplot(3,1,2);plot(t,v2,'r','LineWidth',0.5);

hold off
xlabel('time [sec]','Interpreter','tex','FontSize',8)
ylabel('velocity [m/s]','Interpreter','tex','FontSize',8)
ylim([-0.004 0.004])
legend('E-0%','E-20%');
set(groot, 'DefaultLegendInterpreter', 'tex')
set(groot, 'DefaultLegendFontSize', 8)

 subplot(3,1,3);plot(t,ac1,'b','LineWidth',0.5);
hold on 
subplot(3,1,3);plot(t,ac2,'r','LineWidth',0.5);

hold off
xlabel('time [sec]','Interpreter','tex','FontSize',8)
ylabel('Acceleration m/s^2','Interpreter','tex','FontSize',8)
ylim([-0.3 0.3])
legend('E-0%','E-20%');
set(groot, 'DefaultLegendInterpreter', 'tex')
set(groot, 'DefaultLegendFontSize', 8)
set(fig.Children,'FontName','Times','FontSize',10,'FontWeight','bold');
fig.PaperPositionMode='auto';

savefig(gcf,'images\fig_01.fig')
print('images\fig_01','-dsvg','-r1200')