function G=glocfunc(xxl,lam,typ)

%
% Usage: G=glocfunc(xxl,lam)
%
% This is a local function using the box scheme.
% 
% The rows of xxl are xx+ and xx respectively.
%
% Set typ=1 for Cartesian, anything else for transformed.



% Get a bunch of stuff that we need in any case.
xtan=xxl(1,:)-xxl(2,:);
phi = funcy(0.5*(xxl(1,:)+xxl(2,:)),lam);
hi = norm(xtan)/norm(phi);

% If Cartesian,
if typ==1

	G=hi*phi-xtan';

% Othersise, go for transformed.
else

	fi=[xtan(2), -xtan(1)]*phi/norm(xtan);

	G= hi*[ fi; xtan*phi/norm(xtan) - norm(phi)];

end
