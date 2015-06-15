function G=disfunc(xx,nb,lam,pars)

%
% Usage: G=disfunc(xx,nb,lam,pars)
%
% This is the discrete function that we must set to 
% zero in order to find the invariant torus.
%
% pars = [ww,alp,ofst]
%
% G will be a 2*N - length column vector.

N=size(xx,1);
G=zeros(2*N,1);

for ii=1:N

	ip=mod(ii,N)+1;

	% Get the function value.
	G(2*ii-1:2*ii) = lcfunc(xx([ii,ip],:),nb([ii,ip],:),lam,pars);

end

