function G=gdfunc(xx,lam,typ)

%
% Usage: G=gdfunc(xx,lam)
%
% This is the discrete function that we must set to 
% zero in order to find the invariant torus.
%
% Set typ = 1 for Cartesian, anything else for transformed.


pts=size(xx,1);
G=zeros(pts,1);

for ii=1:pts

	ip=mod(ii,pts)+1;

	% Get the function value.
	G(2*ii-1:2*ii) = glocfunc(xx([ip,ii],:),lam,typ);

end

