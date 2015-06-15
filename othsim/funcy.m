function xpy=funcy(xx,lam,alp,eps0)

%
% This is a vector field for an 8-D system.
% This vector field comes from the Watanabe and Othmer paper.
%
xpy=zeros(8,1);

% The limit cycles.
xpy(1) = alp(1)*xx(1)-       xx(2)-xx(1)*(xx(1)^2+xx(2)^2)+lam*(xx(7)-xx(1));
xpy(2) =        xx(1)+alp(1)*xx(2)-xx(2)*(xx(1)^2+xx(2)^2)+lam*(xx(8)-xx(2));
xpy(3) = alp(2)*xx(3)-       xx(4)-xx(3)*(xx(3)^2+xx(4)^2)+lam*(xx(7)-xx(3));
xpy(4) =        xx(3)+alp(2)*xx(4)-xx(4)*(xx(3)^2+xx(4)^2)+lam*(xx(8)-xx(4));
xpy(5) = alp(3)*xx(5)-       xx(6)-xx(5)*(xx(5)^2+xx(6)^2)+lam*(xx(7)-xx(5));
xpy(6) =        xx(5)+alp(3)*xx(6)-xx(6)*(xx(5)^2+xx(6)^2)+lam*(xx(8)-xx(6));

% The end.
xpy(7) = eps0*((xx(1)+xx(3)+xx(5))/3 - xx(7));
xpy(8) = eps0*((xx(2)+xx(4)+xx(6))/3 - xx(8));

