function [nb, Ubs, Uhs]=normals(xx)

%
% Usage: [nb, Ubs, Uhs]=normals(xx)
%
% This routine calculates unit normal vectors at each point
% just like in disfunc.m.
%
% nb = nbars.
% Ubs = Ubars.
% Uhs = U updated.
%
% See the notes on or about Sept. 14 for details about what these mean.
%
%


pts=length(xx);
nb=zeros(pts,2);

for ii=1:pts

	im=mod(ii-2,pts)+1;			% Get the ii-1 index.
	ip=mod(ii,pts)+1;			% Get the ii+1 index.

	% Tangent "bar".
	xtan=xx(ip,:)-xx(im,:);  xtan=xtan/norm(xtan);
	nb(ii,:) = [ xtan(2) , -xtan(1) ];

	% U "bar".
	Ubs(:,:,ii) = [ xtan; nb(ii,:)]';

	% U "new" (i.e., at half points).
	 xtan=xx(ip,:)-xx(ii,:); xtan=xtan/norm(xtan);
	Uhs(:,:,ii) = [xtan; xtan(2), -xtan(1)]';

end

