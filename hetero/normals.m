function nms=normals(x)

% This routine calculates unit normal vectors at each point, for
% use when moving to a new lambda value.

[pts, tmp]=size(x);
nms=zeros(pts,2);

for ii=0:pts-1

	im=mod(ii-1,pts);	% Get the ii-1 index.
	ip=mod(ii+1,pts);	% Get the ii+1 index.

	xtan=x(ip+1,:)-x(im+1,:);		% Get the tangent and normalize.
	xtan=xtan/norm(xtan);

	nms(ii+1,:)= [ xtan(2) , -xtan(1) ];

end

