function G=dfunc(xx,lam,aa)

% This is the discrete function that we must set to 
% zero in order to find the invariant torus.

pts=length(xx);
G=zeros(pts(1),1);

for ii=1:pts(1)

	ip=mod(ii,pts(1))+1;

	% Get the function value.
	G(ii) = locfunc(xx([ip,ii],:),lam,aa);

end

