function G=disfunc(xx,N,avN,lam,pts,alp,eps0)

% This is the discrete function that we must set to 
% zero in order to find the invariant torus.

[tp,dim]=size(xx); [ug,lee]=size(N); 		% Get the sizes.
cdm = round(lee/dim);						% Co-dimension.

G=zeros(cdm*tp,1);			% Allocate some memory.

% This calculation should only occur once:
pjpk=pts(2)*pts(3);

pos=0;
for ii=1:pts(1)
	
	% Get the sections.
	isc = pjpk*(ii-1); ipsc = mod(ii,pts(1))*pjpk;

	for jj=1:pts(2)
	
		% Ditto.
		jsc = pts(3)*(jj-1); jpsc = mod(jj,pts(2))*pts(3); 

		for kk=1:pts(3)

			kp = mod(kk,pts(3))+1;

			% Remember the ordering in lcfunc.	
			pos=pos+1;
			crnr(1:4) = [pos,ipsc+jsc+kk,ipsc+jpsc+kk,isc+jpsc+kk];
			crnr(5:8) = crnr(1:4)-kk+kp;

			xxl = xx(crnr,:);
						
			% Get the normals in rows.
			oln = reshape(avN(pos,:),dim,cdm)';

			% Get the function value.
			G(cdm*(pos-1)+1:cdm*pos) = lcfunc(xxl,oln,lam,alp,eps0);

		end
	end
end

