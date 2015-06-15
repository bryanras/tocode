function jac=njacob(xx,nb,lam)

%
% Usage: jac=njacob(xx,nb,lam)
%
% This is a numerical computation of the Jacobian of the discrete function.

pts=length(nb);
jac=spalloc(pts,pts,2*pts);

% Get the amount of jiggle needed.
jig=1.0e-9;

% Get the baseline function.
fd=disfunc(xx,lam);

for ii=1:pts

	% Get the indices.
	im = mod(ii-2,pts)+1;
	ip = mod(ii,pts)+1;

	% Jiggle things.
	xx(ii,:)=xx(ii,:)+jig*nb(ii,:);

	% Fill in (pieces of) the matrix.
	jac(im,ii)=(locfunc(xx([ii,im],:),lam)-fd(im))/jig;
	jac(ii,ii)=(locfunc(xx([ip,ii],:),lam)-fd(ii))/jig;

	% Un-jiggle things.
	xx(ii,:)=xx(ii,:)-jig*nb(ii,:);

end


