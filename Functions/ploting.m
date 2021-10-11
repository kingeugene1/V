function [disp,vel,ac]=ploting(t2,u,rx,y,z,ddy,dy,NoN,NL,x_x,m,p)
%visualization
n=size(t2,2);
DoF=3*size(NL,1);
fprintf('outputting the simulation results...\n')
fig1=figure(1);
subplot(3,1,1); plot(t2,u(1:4,:),'-r','LineWidth',0.3);  hold on 
% legend('Input profile');
 plot(t2,rx(1:4,:),'-b','LineWidth',0.3); %legend('Road profile');

hold off 
ylim([-0.005 0.005]);
subplot(3,1,2); plot(t2,z(1:4,:),'-r','LineWidth',0.3);legend('Body vibration');
ylim([-0.005 0.005]);
subplot(3,1,3); plot(t2,z(5:8,:),'-r','LineWidth',0.3);legend('Unsprung mass vibration');
ylim([-0.005 0.005]);
id=1:3:DoF;
deflect= y(1:3:DoF,:);
acc= ddy(id,:);
velo= dy(id,:);
w=zeros(length(t2),1);
nm=round(NoN/2);
mag=1000;
ac= acc(nm,:);
vel=velo(nm,:);
disp= deflect(nm,:);
fig2=figure(2);
subplot(3,1,1);plot(t2,vel,'-k','LineWidth',0.3);
xlabel('time [sec]')
ylabel('Velocity')
legend('Velocity')

subplot(3,1,2);plot(t2,ac,'-k','LineWidth',0.3);
xlabel('time [sec]')
ylabel('acceleration')
legend('acceleration')

subplot(3,1,3);plot(t2,disp,'-k','LineWidth',0.3); 
xlabel('time [sec]')
ylabel('Displacement')
legend('Displacement')
% ylim([-0.005 0.0015])
fig1.Units='centimeters';
% fig1.Position(3)=8;
% fig1.Position(4)=6;
set(fig1.Children,'FontName','Times','FontSize',8,'FontWeight','bold');

fig2.Units='centimeters';
% fig2.Position(3)=8;
% fig2.Position(4)=6;
set(fig2.Children,'FontName','Times','FontSize',8,'FontWeight','bold');
set(gca,'LooseInset',max(get(gca,'TightInset'),0.02));
fig1.PaperPositionMode='auto';
print('images/my_figure1','-dpng','-r600');
fig2.PaperPositionMode='auto';
print('images/my_figure2','-dpng','-r600');
figure(3);set(gcf, 'Unit', 'pixel', 'Position',[120, 280, 800, 300]);
CT=datestr(now,'yyyy_mm_dd');
filename=sprintf('VBI_%s',CT);
writerObj = VideoWriter(filename,'MPEG-4'); %-- setting movie class
writerObj.FrameRate = 20; %-- flame rate (slide/sec)
open(writerObj); %-- initialization of movie object
for tt=1:10:n
    for ii=1:m+1
            xlines=(ii-1)*(p+1)+(1:p+1);
            for j=1:p+1
            ylines=j:p+1:NoN;
            
            plot3(NL(xlines,1),NL(xlines,2),1000*deflect(xlines,tt),'b-')
            hold on
            plot3(NL(ylines,1),NL(ylines,2),1000*deflect(ylines,tt),'b-')
            end
    end
    plot3(NL(:,1),NL(:,2),1000*deflect(:,tt),'r.')
    plot3(x_x(1:2,1:tt)',w(1:tt).*ones(1,2)+0.68,mag*z(5:6,1:tt)'+2,'r-');
    plot3(x_x(1:2,tt)',w(tt).*ones(1,2)+0.68,mag*z(5:6,tt)'+2','ko');
    
    plot3(x_x(3:4,1:tt)',w(1:tt).*ones(1,2)+3.08,mag*z(7:8,1:tt)'+2,'r-');
    plot3(x_x(3:4,tt)',w(tt).*ones(1,2)+3.08,mag*z(7:8,tt)'+2','ko');
    text(10, 5, -2.0, ['TIME=' num2str(t2(tt)) ' [s]'])
    hold off
    set(gcf,'Position',[120,141,733,613]);
    axis tight;
    set(gca,'CameraPosition',[-152 -18  136],'CameraTarget',[10     2    -4],'zlim',[-15 7],'xlim',[-20 40])
    frame = getframe(figure(3)); %-- get the figure as a frame of the animation
    writeVideo(writerObj,frame); %-- add the gotten frame to the movie object
    
end
end