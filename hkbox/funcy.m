function xpy=funcy(x,b,c,d,lam)

%
% This is a vector field for a 3-D system.
% The ODE comes from Hale-Kocak, p. 515. 
% Langford is the earliest reference, though.
%

xpy=zeros(3,1);

b=lam-b;
d = x(3)+d*(1.0-x(3)^2);

xpy(1)=b*x(1)-c*x(2)+x(1)*d;
xpy(2)=c*x(1)+b*x(2)+x(2)*d;
xpy(3)=lam*x(3) - x(1)^2 - x(2)^2 - x(3)^2;
