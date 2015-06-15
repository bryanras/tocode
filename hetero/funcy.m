function xpy=funcy(x,lam)

%
% This is a vector field for a 2-D system.
% This vector field has a heteroclinic cycle.
%

xpy=zeros(2,1);
ugh2=x(1)^2+x(2)^2;
ugh=sqrt(ugh2);

xpy(1)=x(1)*(lam-ugh2) - x(2)^2/ugh;
xpy(2)=x(2)*(lam-ugh2+x(1)/ugh);
