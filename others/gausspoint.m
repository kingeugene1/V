function [x,weight]=gausspoint()

% this function return the gaussian quadrature values
% gp shows at which gaussian point your at , it extract xi,eta and weight
% info.
%
% Second order gaussian
    x=[-1/sqrt(3) -1/sqrt(3);
        -1/sqrt(3)  1/sqrt(3);
        1/sqrt(3) -1/sqrt(3);
        1/sqrt(3) 1/sqrt(3)];
    weight=[1 1;
                  1 1;
                  1 1;
                  1 1];
% x=[-sqrt(3/7-2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5));
%       sqrt(3/7-2/7*sqrt(6/5))  sqrt(3/7+2/7*sqrt(6/5));
%       -sqrt(3/7+2/7*sqrt(6/5))  -sqrt(3/7+2/7*sqrt(6/5));
%       sqrt(3/7+2/7*sqrt(6/5))  sqrt(3/7+2/7*sqrt(6/5))];
%     weight=[(18+sqrt(30))/36  (18-sqrt(30))/36;
%                   (18+sqrt(30))/36  (18-sqrt(30))/36;
%                   (18+sqrt(30))/36  (18-sqrt(30))/36;
%                   (18+sqrt(30))/36  (18-sqrt(30))/36];

end