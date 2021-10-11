function[rx,ry,time]= moving_load(Lapp,v,w,A,Lb,tk,inv)
%Lapp= approach distance
%v= vehicle speed
%w= wheel base matrix [[wbx1 wbxy1];[wbx2 wby2];......] (m)
%A= cell of axle spacing
%Lb= bridge length
%tk= time step
%inv y coordinates of vehicle path entry
%rx x coordinates of force on plate
%ry y coordinates of force on plate
% time computed for plate solver
% Time
t1 = min((Lapp)./abs(v));
t2 = max((Lapp+Lb(1)+w(:,1)')./abs(v));
time = (0:tk:t2-t1);
% R Location on the plate (Rx and Ry)
rx = []; 
ry = [];
for j1 = 1:length(A)
    % X position of R
    aux3 = [];
    if whbs(j1,2) > 0
        for i = 1:length(A{j1})/2
            aux1 = (time+(t1-Lapp/abs(v)))*abs(v(j1))-A{j1}(i);
            if velj(j1) < 0
                aux1 = -aux1+Lb(1);
            end
            aux2 = find(aux1<=0 | aux1>=Lb(1)); 
            aux1(aux2)=0;
            aux3(i,:) = aux1;
        end
        Rx = [ Rx; aux3; aux3 ];
    elseif whbs(j1,2) == 0
%         aux1 = (time_pla+(t1-L_app(j1)/abs(velj(j1))))*abs(velj(j1))-Ax_sp{j1}(1);
%         if velj(j1) < 0; aux1 = -aux1+Lb(1); end
%         aux2 = find(aux1<=0 | aux1>=Lb(1)); aux1(aux2)=0;
%         aux3(1,:) = aux1;
%         Rx = [ Rx; aux3 ];
        for i = 1:length(Ax_sp{j1})
            aux1 = (time_pla+(t1-L_app(j1)/abs(velj(j1))))*abs(velj(j1))-Ax_sp{j1}(i);
            if velj(j1) < 0; aux1 = -aux1+Lb(1); end
            aux2 = find(aux1<=0 | aux1>=Lb(1)); aux1(aux2)=0;
            aux3(i,:) = aux1;
        end
        Rx = [ Rx; aux3 ];
    end % if whbs(j1,2)
        
    % Y position of R
    aux1 = InVehj(j1);
    if whbs(j1,2) > 0
        aux1 = [aux1*ones(1,length(Ax_sp{j1})/2) (aux1+whbs(j1,2))*ones(1,length(Ax_sp{j1})/2)];
    elseif whbs(j1,2) == 0
        aux1 = [aux1*ones(1,length(Ax_sp{j1}))];
    end % whbs(j1,2) == 0
    if velj(j1) < 0; aux1 = fliplr(aux1); end
    Ry = [ Ry, aux1 ];
end
end