function G=dfunc(xx,N,ptsi,ptsj,b,c,d,lam)

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
		xll = xx(sii+jj,:); xur = xx(sip+jp,:);
		xlr = xx(sii+jp,:); xul = xx(sip+jj,:);

		% Get the vector field and function value.
		ff=funcy(0.25*(xll+xur+xlr+xul), b,c,d,lam);
		v = locn(xll,xlr,xul,xur);
		G(sii+jj) = v*ff;

%		G(sii+jj) = lcfunc(xll,xlr,xul,xur,b,c,d,lam);

		% Added stuff.  See notes on 6/3/03.
%		Nbo=mean(N([sii+jj,sip+jj,sii+jp,sip+jp],:),1)';
%		Nbo=Nbo/norm(Nbo);
%		G(sii+jj) = G(sii+jj)*(v*Nbo); 

	end

end

