function enn2=nnmls(xx,N,pts)

% This is a function for updating the normal vectors at each new
% lambda step.

% Get the dimensions.
[aps,ndim]=size(N);
[aps,dim]=size(xx);
cdm=round(ndim/dim);

% Allocate.
enn2=zeros(aps,ndim);

% Use a center difference at each point.
isk=pts(2)*pts(3);
for ii=1:pts(1)

	iisc = (ii-1)*isk;
	ipsc = mod(ii,pts(1))*isk;
	imsc = mod(ii-2,pts(1))*isk;

	for jj=1:pts(2)

		jjsc = (jj-1)*pts(3);
		jpsc = mod(jj,pts(2))*pts(3);
		jmsc = mod(jj-2,pts(2))*pts(3);

		for kk=1:pts(3)

			km = mod(kk-2,pts(3))+1;
			kp = mod(kk,pts(3))+1;

			% Get the tangent directions.
			tti = xx(ipsc+jjsc+kk,:)-xx(imsc+jjsc+kk,:);
			ttj = xx(iisc+jpsc+kk,:)-xx(iisc+jmsc+kk,:);
			ttk = xx(iisc+jjsc+kp,:)-xx(iisc+jjsc+km,:);

			% Get the old normals.
%			olnt = mean(N([ 	iisc+jjsc+kk,ipsc+jjsc+kk,imsc+jjsc+kk, ...
%								iisc+jpsc+kk,iisc+jmsc+kk,iisc+jjsc+kp,...
%								iisc+jjsc+km ],:),1);

			olnt = N(iisc+jjsc+kk,:);

			oln = reshape(olnt,dim,cdm)';
	
			% Orthonormalize.
			[Q,R]=qr([tti;ttj;ttk;oln]');
			oln(1,:) = sign(R(4,4))*Q(:,4)';
			oln(2,:) = sign(R(5,5))*Q(:,5)';
			oln(3,:) = sign(R(6,6))*Q(:,6)';
			oln(4,:) = sign(R(7,7))*Q(:,7)';
			oln(5,:) = sign(R(8,8))*Q(:,8)';
			
			% Store it back in the original normal vector slot.
			enn2(iisc+jjsc+kk,:)=reshape(oln',1,ndim);

		end
	end
end

