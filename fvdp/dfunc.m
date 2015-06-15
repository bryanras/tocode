function G=dfunc(xx,ptsi,ptsj,lam,pars)

% This is the discrete function that we must set to 
% zero in order to find the invariant torus.
G=zeros(ptsi*ptsj,1);

% Normal method
for ii=1:ptsi

	sii=(ii-1)*ptsj;			% Get the ii section.
	sip=mod(ii,ptsi)*ptsj;		% Get the ii+1 section.

	for jj=1:ptsj

		jp=mod(jj,ptsj)+1;      % Get the jj+1 index.
			
		% Get the box of points.
		xll = xx(sii+jj,:);
		xur = xx(sip+jp,:);
		xlr = xx(sii+jp,:);
		xul = xx(sip+jj,:);

		% Get the vector field and function value.
		ff=funcy(0.25*(xll+xur+xlr+xul),lam,pars);
		G(sii+jj) = locn(xll,xlr,xul,xur)*ff;

	end

end

