function DP=dPhi(xx,lam)

% Simple oscillator.
% DP=[ 	lam-3*xx(1)^2 -xx(2)^2,		-1-2*xx(1)*xx(2);     ...
%     	1-2*xx(1)*xx(2), 				lam-xx(1)^2-3*xx(2)^2 ];

% van der Pol  oscillator.
 DP=[ 	0,				1; 			...
     	-1-2*lam*xx(1)*xx(2), 		lam*(1-xx(1)^2) ];


