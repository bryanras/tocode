function xpy=funcy(xx,lam)

% Usage: xpy=funcy(xx,lam)
%
% This is a vector field for the Lorenz system.
%
% See notes on 12/4/03. Non-continuation parameter values are
%
%

global sig b;

xpy=zeros(3,1);

xpy(1)=sig*(xx(2)-xx(1));

xpy(2)=lam*xx(1)-xx(2)-xx(1)*xx(3);

xpy(3)=xx(1)*xx(2)-b*xx(3);
