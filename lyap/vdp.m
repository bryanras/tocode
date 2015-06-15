% function f=vdp(x)%%%%
% Van der Pol vector field 
function f=vdp(x)

% Set the parameter.
%para=1; 			% Default.
para=0.6; 	


f=[x(2); -x(1)+para*x(2)*(1-x(1)^2)]; % van der Pol oscillator
