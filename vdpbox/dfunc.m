function Df = dfunc(xx,lam)

%
% Usage: Df = dfunc(xx,lam)
%
% This is a function for calculating the linearized vector field.
%

Df = [0, 1; -1-2*lam*xx(2)*xx(1), lam*(1-xx(1)^2)];

% If you want to reverse the vector field, ...
%Df=-Df;
