function nms=normals(x)

% This routine calculates unit normal vectors at each point
% using a quadratic derivative, just like in dfunc.m.

% fittype=0 => linear tangent
% fittype=1 => quadratic tangent
fittype=0;

pts=length(x);
nms=zeros(pts,2);

for ii=0:pts-1

	im=mod(ii-1,pts);	% Get the ii-1 index.
	ip=mod(ii+1,pts);	% Get the ii+1 index.
	xi=x(ii+1,:);
	xp=x(ip+1,:);
	xm=x(im+1,:);

	% Quadratic fit
	if fittype==1
		s=norm(xi-xm)/(norm(xp-xi)+norm(xi-xm));
		aa=(xi-xm-s*(xp-xm))/(s^2-s);
		bb=xp-xm-aa;
		xtan=2*s*aa+bb;

	% Linear fit
	elseif fittype==0
		xtan=xp-xm;
	end

	xtan=xtan/norm(xtan); % Normalize

	nms(ii+1,:)= [ xtan(2) , -xtan(1) ];

end

