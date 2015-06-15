% function Df=Dvdp(x)%%%%
% Derivative of Van der Pol vector field
function Df=Dvdp(x)

% Set the parameter.
%para=1; % Default.
para=0.6;

% Jacobian of van der Pol oscillator
Df=[0 1; -1-2*para*x(1)*x(2) para*(1-x(1)^2)]; 
