function nms=normals(xx,ptsi,ptsj,typ)


% This routine calculates unit normal vectors at each point.

% Loop through the i's.
for ii=0:ptsi-1
	sii=ii*ptsj;					% Get the ii section.
	sim=mod(ii-1,ptsi)*ptsj;		% Get the ii-1 section.
	sip=mod(ii+1,ptsi)*ptsj;		% Get the ii+1 section.

	% Loop through the j's.
	for jj=1:ptsj
		jm=mod(jj-2,ptsj)+1;	% Get the jj-1 index.
		jp=mod(jj,ptsj)+1;		% Get the jj+1 index.

		% Define the nearby points.

		% Normals to the center of the box.
		% Used for calculating a function value.
		xll = xx(sii+jj,:);
		xlr = xx(sii+jp,:);
		xul = xx(sip+jj,:);
		xur = xx(sip+jp,:);
			
		% Assign a normal vector.
		nms(sii+jj,:) = locn(xll,xlr,xul,xur); 

	end
end

% If we are getting the NEW nbars, then do an average.
if (typ==0)
	for ii=0:ptsi-1
		sii=ii*ptsj;					% Get the ii section.
		sim=mod(ii-1,ptsi)*ptsj;		% Get the ii-1 section.
		sip=mod(ii+1,ptsi)*ptsj;		% Get the ii+1 section.

		% Loop through the j's.
		for jj=1:ptsj
			jm=mod(jj-2,ptsj)+1;	% Get the jj-1 index.
			jp=mod(jj,ptsj)+1;		% Get the jj+1 index.
	
			% Assign a normal vector.
			nmt(sii+jj,:) = mean(nms([sii+jj,sii+jm,sim+jj,sim+jm],:),1);

		end
	end

	nms=nmt;
end

