function N=nnmls(xx,ptsi,ptsj)

% This routine calculates unit normal vectors at each point.

% Loop through the i's.
for ii=1:ptsi
	sii=(ii-1)*ptsj;				% Get the ii section.
	sim=mod(ii-2,ptsi)*ptsj;		% Get the ii-1 section.
	sip=mod(ii,ptsi)*ptsj;			% Get the ii+1 section.

	% Loop through the j's.
	for jj=1:ptsj
		jm=mod(jj-2,ptsj)+1;	% Get the jj-1 index.
		jp=mod(jj,ptsj)+1;		% Get the jj+1 index.

		% Assign a normal vector.
		xll = xx(sim+jj,:); xlr = xx(sii+jp,:);
		xul = xx(sii+jm,:); xur = xx(sip+jj,:);

		N(sii+jj,:) = locn(xll,xlr,xul,xur); 

	end
end
