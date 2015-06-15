function G=disfunc(xx,enns,avenns,pts,alp,lam)

%
% Usage: G = disfunc(xx,enns,avenns,pts,alp,lam)
%
% This is the discrete function that we must set to 
% zero in order to find the invariant torus.
%

[tp,dim]=size(xx); [ug,lee]=size(enns); 		% Get the sizes.
cdm = round(lee/dim);							% Co-dimension.

G=zeros(cdm*tp,1);								% Allocate some memory.

% Loop through.
pos=0;
for ii=1:pts(1)
	
	% Get the sections.
	isc = pts(2)*(ii-1); ipsc = mod(ii,pts(1))*pts(2);

	for jj=1:pts(2)
	
		jp = mod(jj,pts(2))+1; 

		% Remember the ordering in lcfunc.	
		pos=pos+1;
		xxl = xx([isc+jj,ipsc+jj,isc+jp,ipsc+jp],:);
						
		% Get the normals in rows.
		oln = reshape(avenns(pos,:),dim,cdm)';

		% Get the function value.
		G(cdm*(pos-1)+1:cdm*pos) = lcfunc(xxl,oln,alp,lam);

	end
end

