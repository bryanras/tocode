function jac=njacob(xx,nb,ptsi,ptsj,lam,pars)

%
% Usage:  jac=njacob(xx,nb,ptsi,ptsj,lam,pars)
%
% pars = [omega,alpha,ofst]
%
% This is a numerical computation of the Jacobian of the discrete function.

% Determine if we are using forward or backwards differencing.
% updn=1 => fwd; updn=0 => bkwd;
updn=1;

% Allocate some memory.
jac=spalloc(ptsi*ptsj,ptsi*ptsj,5*ptsi*ptsj);
fut=zeros(4,1);

% Get the amount of jiggle needed.
jig=updn*1.0e-8;

% Get the baseline function.
fd=dfunc(xx,ptsi,ptsj,lam,pars);

for ii=1:ptsi

	% Get the indices in the plus/minus i-direction.
	im = mod(ii-2,ptsi);
	ip = mod(ii,ptsi);

	% Figure out what ptsj section we are in.
	sim=im*ptsj;
	sii=(ii-1)*ptsj;
	sip=ip*ptsj;

	for jj=1:ptsj

		% Get the indices in the plus/minus j-direction.
		jm = mod(jj-2,ptsj)+1;
		jp = mod(jj,ptsj)+1;

		% Jiggle the point.
		xx(sii+jj,:)=xx(sii+jj,:)+jig*nb(sii+jj,:);

		% Work on a piece of the torus each time.
		% jumpy is how we jump around the rectangles. 
		% Each row is a different square.
		jumpy=[	sim+jm,		sim+jj,		sii+jm,		sii+jj; ...
				sim+jj,		sim+jp,		sii+jj,		sii+jp; ...
				sii+jm,		sii+jj,		sip+jm,		sip+jj; ...
				sii+jj,		sii+jp,		sip+jj,		sip+jp ];

		% Note that we only have four nonzero entries in each column.
		for kk=1:4

			% Give everything a name to make debugging easier.
			xll=xx(jumpy(kk,1),:); 	xlr=xx(jumpy(kk,2),:); 	
			xul=xx(jumpy(kk,3),:); 	xur=xx(jumpy(kk,4),:); 	
			xave=0.25*(xll+xlr+xul+xur);

			% Get the normal.
			v=locn(xll,xlr,xul,xur);
			
			% Calculate the vector field and function.
    		fut(kk)=v*funcy(xave,lam,pars);

		end

		% Reset the original torus.
		xx(sii+jj,:)=xx(sii+jj,:)-jig*nb(sii+jj,:);

		% Fill in (pieces of) the matrix.
		for kk=1:4
			jac(jumpy(kk,1),sii+jj) = (fut(kk)-fd(jumpy(kk,1)))/jig;
		end

	end

end

