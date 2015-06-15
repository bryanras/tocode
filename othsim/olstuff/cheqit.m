function [maxi, maxj, maxk]=cheqit(pts)

% This is a function for checking smoothness of the updated normal vectors.

global xx enns;

% Get the dimensions.
[aps,ndim]=size(enns);
[aps,dim]=size(xx);
cdm=round(ndim/dim);

maxi=0; maxj=0; maxk=0;
% Use a center difference at each point.
isk=pts(2)*pts(3);
for ii=1:pts(1)

	iisc = (ii-1)*isk;
	ipsc = mod(ii,pts(1))*isk;

	for jj=1:pts(2)

		jjsc = (jj-1)*pts(3);
		jpsc = mod(jj,pts(1))*pts(3);

		for kk=1:pts(3)

			kp = mod(kk,pts(3))+1;
			pos=iisc+jjsc+kk;

			% Check each of the norms.
			tmp=norm(enns(ipsc+jjsc+kk,:)-enns(pos,:));
			if tmp > maxi maxi=tmp; end

			tmp=norm(enns(iisc+jpsc+kk,:)-enns(pos,:));
			if tmp > maxj maxj=tmp; end

			tmp=norm(enns(iisc+jjsc+kp,:)-enns(pos,:));
			if tmp > maxk maxk=tmp; end
	
		end
	end
end

