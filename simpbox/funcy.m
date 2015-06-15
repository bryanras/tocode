function xpy=funcy(x,lam)

%
% This is a vector field for a 2-D system.
%


%
% This vector field has a simple circular limit cycle of radius lambda^2.
%
%xpy = [lam*x(1) -     x(2) - x(1)*(x(1)^2+x(2)^2); ...
%       x(1) + lam*x(2) - x(2)*(x(1)^2+x(2)^2)];

%
% This is the van der Pol vector field.
%
xpy = [x(2); -x(1) + lam*x(2)*(1-x(1)^2)];

