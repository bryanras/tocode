function xpy=funcy(xx,alp,lam)

%
% This is a vector field for an 8-D system.
% This vector field comes from the Watanabe and Othmer paper.
%
xpy=zeros(4,1);

% Don't compute more than once.
cp = lam*(xx(1)+xx(2)-xx(3)-xx(4));
sq1=(1- xx(1)^2 - xx(2)^2);
sq2=(1- xx(3)^2 - xx(4)^2);

% The limit cycles.
xpy(1) =  xx(1)*sq1 + xx(2)*alp - cp;
xpy(2) = -xx(1)*alp + xx(2)*sq1 - cp;
xpy(3) =  xx(3)*sq2 + xx(4)*alp + cp;
xpy(4) = -xx(3)*alp + xx(4)*sq2 + cp;

