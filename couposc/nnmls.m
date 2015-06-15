function enn2=nnmls(xx,enns,pts)

% This is a function for updating the normal vectors at each new
% lambda step.

% Get the dimensions.
[aps,ndim]=size(enns);
[aps,dim]=size(xx);
cdm=round(ndim/dim);

% Allocate.
enn2=zeros(aps,ndim);

% Use a center difference at each point.
for ii=1:pts(1)

	iisc = (ii-1)*pts(2);
	ipsc = mod(ii,pts(1))*pts(2);
	imsc = mod(ii-2,pts(1))*pts(2);

	for jj=1:pts(2)

		jp = mod(jj,pts(2))+1;
		jm = mod(jj-2,pts(2))+1;

		% Get the tangent directions.
		tt1 = xx(ipsc+jj,:)-xx(imsc+jj,:);
		tt2 = xx(iisc+jp,:)-xx(iisc+jm,:);

		% Get the old normals.
		olnt = enns(iisc+jj,:);

		oln = reshape(olnt,dim,cdm)';
	
		% Orthonormalize.
		[Q,R]=qr([tt1;tt2;oln]');
		oln(1,:) = sign(R(3,3))*Q(:,3)';
		oln(2,:) = sign(R(4,4))*Q(:,4)';
			
		% Store it back in the original normal vector slot.
		enn2(iisc+jj,:)=reshape(oln',1,ndim);

	end
end

