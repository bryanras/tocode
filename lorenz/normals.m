function nn=normals(xx,nb)

%
% Usage: nn=normals(xx,nn)
%
% This routine calculates unit normal vectors at each point
% just like in disfunc.m.
%
% nb = N x 6 matrix column-split into two normals for each ii.
%
% The calculation method is QR decomposition with old normals and center-
% difference tangents.
%

N=size(xx,1);
nn=zeros(N,2*size(xx,2));

for ii=1:N

	im=mod(ii-2,N)+1;		% Get the ii-1 index.
	ip=mod(ii,N)+1;			% Get the ii+1 index.

	% Tangent "bar".
	tn=xx(ip,:)-xx(im,:);
	[Q,R]   = qr( [tn; nb(ii,1:3); nb(ii,4:6)]' );

	nn(ii,1:3) = sign(R(2,2))*Q(:,2)';
	nn(ii,4:6) = sign(R(3,3))*Q(:,3)';

end

