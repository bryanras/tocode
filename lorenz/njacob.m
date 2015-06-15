function JJ=njacob(xx,nb,lam)

%
% Usage: JJ=njacob(xx,nb,lam)
%
% This is a numerical computation of the Jacobian of the discrete function.
%
% The Jacobian will be 2*N^2.

N=size(xx,1);
JJ=spalloc(2*N,2*N,4*N);

% Get the amount of jiggle needed.
jig=1.0e-9;

% Get the baseline function.
fd=disfunc(xx,nb,lam);

% Skip through the Jacobian.
dim = 3; cdm = 2;
for ii=1:N

	% Get the indices.
	im = mod(ii-2,N)+1;
	ip = mod(ii,N)+1;

	% Jiggle things twice.
	for jj=1:cdm

		xx(ii,:)=xx(ii,:)+jig*nb(ii,(jj-1)*dim+1:jj*dim);

		% Fill in (pieces of) the matrix.
		JJ((im-1)*cdm+1:im*cdm, (ii-1)*cdm+jj)= ...
			(lcfunc(xx([im,ii],:),nb([im,ii],:),lam)- ...
			fd((im-1)*cdm+1:im*cdm))/jig;

		JJ((ii-1)*cdm+1:ii*cdm, (ii-1)*cdm+jj)= ...
			(lcfunc(xx([ii,ip],:),nb([ii,ip],:),lam)- ...
			fd((ii-1)*cdm+1:ii*cdm))/jig;

		% Un-jiggle things.
		xx(ii,:)=xx(ii,:)-jig*nb(ii,(jj-1)*dim+1:jj*dim);

	end

end


