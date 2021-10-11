function [disp,vel,ac]= ploting(t,u,r,y,z,ddy,dy,NoN,NL,x,m,p)
%visualization
n=size(t,2);
DoF=3*size(NL,1);
fprintf('outputting the simulation results...\n')
figure(1)
subplot(3,1,1); plot(t,u(1:4,:),'-r'); legend('Input profile'); hold on
subplot(3,1,1); plot(t,r(1:4,:),'-b'); legend('Road profile');hold off 

subplot(3,1,2); plot(t,z(1:4,:),'-r');legend('Body vibration');
subplot(3,1,3); plot(t,z(5:8,:),'-r');legend('Unsprung mass vibration');
   ylim([-0.01 0.1])
id=1:3:DoF;
deflect= y(1:3:DoF,:);
acc= ddy(id,:);
velo= dy(id,:);
w=zeros(length(t),1);
nm=round(NoN/2);
mag=1000;
ac= acc(nm,:);
vel=velo(nm,:);
disp= deflect(nm,:);
figure(2)
subplot(3,1,1);plot(t,vel,'-r');
xlabel('time [sec]','Interpreter','tex','FontSize',9)
ylabel('Velocity','Interpreter','tex','FontSize',9)
legend('Velocity')

subplot(3,1,2);plot(t,ac,'-r');
xlabel('time [sec]','Interpreter','tex','FontSize',9)
ylabel('acceleration','Interpreter','tex','FontSize',9)
legend('acceleration')

subplot(3,1,3);plot(t,disp,'-r'); 
xlabel('time [sec]','Interpreter','tex','FontSize',9)
ylabel('Displacement','Interpreter','tex','FontSize',9)
legend('Displacement')
% ylim([-0.005 0.0015])
figure(3);set(gcf, 'Unit', 'pixel', 'Position',[120, 280, 800, 300]);
writerObj = VideoWriter('VBI.MPEG-4','MPEG-4'); %-- setting movie class
writerObj.FrameRate = 20; %-- flame rate (slide/sec)
open(writerObj); %-- initialization of movie object
for tt=1:10:n
    for ii=1:m+1
            xlines=(ii-1)*(p+1)+(1:p+1);
            for j=1:p+1
            ylines=j:p+1:NoN;
            
            plot3(NL(xlines,1),NL(xlines,2),6000*deflect(xlines,tt),'b-')
            hold on
            plot3(NL(ylines,1),NL(ylines,2),6000*deflect(ylines,tt),'b-')
            end
    end
    plot3(NL(:,1),NL(:,2),6000*deflect(:,tt),'r.')
    plot3(x([1 3],1:tt)',w(1:tt).*ones(1,2)+0.68,mag*z(5:6,1:tt)'+2,'m-');
    plot3(x([1 3],tt)',w(tt).*ones(1,2)+0.68,mag*z(5:6,tt)'+2','mo');
    
    plot3(x([2 4],1:tt)',w(1:tt).*ones(1,2)+2.68,mag*z(7:8,1:tt)'+2,'m-');
    plot3(x([2 4],tt)',w(tt).*ones(1,2)+2.68,mag*z(7:8,tt)'+2','ko');
    text(10, 5, -2.0, ['TIME=' num2str(t(tt)) ' [s]'])
    hold off
    set(gcf,'Position',[120,141,733,613]);
    set(gca,'CameraPosition',[-194 -16  112])
    set(gca,'CameraTarget',[15 2 -1.5])
    set(gca,'zlim',[-10 7])
    set(gca,'xlim',[-20  50])
    frame = getframe(figure(3)); %-- get the figure as a frame of the animation
    writeVideo(writerObj,frame); %-- add the gotten frame to the movie object
end
close(writerObj); %-- end of animation
end
