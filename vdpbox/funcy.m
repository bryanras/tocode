function xpy=funcy(xx,lam)

% Usage: xpy=funcy(xx,lam)
%
% This is a vector field for the van der Pol oscillator in 2D.
%

xpy=zeros(2,1);
xpy(1)=xx(2);
xpy(2)=-xx(1)+lam*xx(2)*(1-xx(1)^2);

% If you want to reverse the vector field, ...
%xpy=-xpy;
