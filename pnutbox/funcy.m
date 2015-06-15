function xpy=funcy(x,lam,aa)

%
% This is a vector field for a 2-D system.
% This vector field has a non-convex (peanut) limit cycle.
%

xpy=zeros(2,1);
xpy(1)=x(1)*(lam-x(1)^2-x(2)^2)+aa*x(1)^3/((x(1)^2+x(2)^2)^(3/2))-x(2);
xpy(2)=x(2)*(lam-x(1)^2-x(2)^2)+aa*x(2)*x(1)^2/((x(1)^2+x(2)^2)^(3/2))+x(1);
