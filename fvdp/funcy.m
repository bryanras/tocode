function xpy=funcy(x,lam,pars)

%
% This is a vector field for a 3-D system.
% The ODE comes from a forced van der Pol oscillator
%
% pars = [ww,alp,ofst]
%

% Hard-code for now.
% Can fix later.
ww = pars(1);	
alp = pars(2);
ofst = pars(3);

xpy=zeros(3,1);

ugh=x(1)^2+x(2)^2;
sugh=sqrt(ugh);

xpy(1)=lam*(x(1)^2)/ugh-x(1)*x(3)/sugh-ww*x(2);
xpy(2)=lam*x(1)*x(2)/ugh-x(2)*x(3)/sugh+ww*x(1);
xpy(3)=sugh-ofst+alp*x(3)*(1-(x(3)^2)/3);
