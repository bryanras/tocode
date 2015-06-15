function rho = testit(xx,nb,lam)

%
% Usage: rho = testit(xx,lam)
%
%
% This is a script for generating the rhos and product of the modified 
% system as suggested in Hale, p. 220. The final rho is the product.
%
% See the notes on 10/1/03 for details.
%

[pts, dm]=size(xx);


% Run through everything.
rho(1) = 1;
for ii=1:pts

	ip = mod(ii,pts)+1;

	% Get the half-point and tangent.
	% These take the place of the exact derivatives and whatnot.
	xtan=xx(ip,:)-xx(ii,:);
	xh = (xx(ip,:)+xx(ii,:))/2;
	nh = (nb(ip,:)+nb(ii,:))/2;
	dphi = norm(xtan);

	% Get the stuff associated with the vector field.
	dPhih = dfunc(xh,lam);
	pn=norm(funcy(xh,lam));

	% This ain't exactly stable, but it will do for now.
	rho(ii+1) = ( (2*pn+dphi*nh*dPhih*nh')/(2*pn-dphi*nh*dPhih*nh') )*rho(ii);


end


