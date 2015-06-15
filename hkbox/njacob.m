function jac=njacob(xx,N,ptsi,ptsj,b,c,d,lam)

% This is a numerical computation of the Jacobian of the discrete function.

% Determine if we are using forward or backwards differencing.
% updn=1 => fwd; updn=0 => bkwd;
updn=1;

% NUMERICAL EXPERIMENT.
%angle=1;
%diffex=0;
%proj=0;

% Allocate some memory.
jac=spalloc(ptsi*ptsj,ptsi*ptsj,5*ptsi*ptsj);

% Get the amount of jiggle needed.
jig=updn*1.0e-8;

% Get the baseline function.
fd=dfunc(xx,N,ptsi,ptsj,b,c,d,lam);

flipchk=0;

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
		xx(sii+jj,:)=xx(sii+jj,:)+jig*N(sii+jj,:);

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

			% NUMERICAL EXPERIMENT.
%			if kk==4
%
%				xll=xll-jig*N(jumpy(kk,1),:);
%
%				% Numerical difference of dn/dr's.
%%				tempyi=(locn(xll,xlr,xul,xur+jig*N(jumpy(kk,4),:)) - ...
%					locn(xll,xlr,xul,xur-jig*N(jumpy(kk,4),:)) )/(2*jig);
%				tempyii=(locn(xll+jig*N(jumpy(kk,1),:),xlr,xul,xur) - ...
%					locn(xll-jig*N(jumpy(kk,1),:),xlr,xul,xur) )/(2*jig);
%				ugly=max([norm(tempyi+tempyii),ugly]);
%
%				% Analytical difference of dn/dr's.
%				tv = cross( xur-xll, xul-xlr );
%				ntv = norm(tv);
%				tempit=cross(N(jumpy(kk,4),:)-N(jumpy(kk,1),:),xul-xlr)'/ntv;
%				ntmp=norm(tempit);
%
%				proj=max([ 1-abs(tv*tempit)/(ntv*ntmp), proj ]);
%				tempit = (eye(3)-(tv'*tv)/ntv^2)*tempit;
%				ntmp=norm(tempit);
%
%				if norm(tempit)>diffex
%					diffex=norm(tempit);
%					indx=jj;
%				end
%
%
%				% Angles.
%				angle=min([ angle, norm(tv)/(norm(xur-xll)*norm(xul-xlr)) ]);
%
%				xll=xll+jig*N(jumpy(kk,1),:);
%
%			end

			% If using Matlab code ...
			v=locn(xll,xlr,xul,xur);
			
			% Check to make sure we are not flipping the normals.
			if (flipchk==0) & (v*N(jumpy(kk,1),:)'<0)
				disp(' Warning: A normal vector is being flipped!')
				flipchk=1;
			end
	  		fut(kk)=v*funcy(xave,b,c,d,lam);

			% Otherwise...
%			fut(kk)=lcfunc(xll,xlr,xul,xur,b,c,d,lam);

			% Added stuff. See notes on 6/3/03.
%			Nbo=mean(N(jumpy(kk,:),:),1)'; Nbo=Nbo/norm(Nbo);
%			fut(kk)=fut(kk)*(locn(xll,xlr,xul,xur)*Nbo); 

		end

		% Reset the original torus.
		xx(sii+jj,:)=xx(sii+jj,:)-jig*N(sii+jj,:);

		% Fill in (pieces of) the matrix.
		for kk=1:4
			jac(jumpy(kk,1),sii+jj) = (fut(kk)-fd(jumpy(kk,1)))/jig;
		end

	end

end

% NUMERICAL EXPERIMENT.
%disp(sprintf('\b Angle: %g \t Diffex: %g \t Proj: %g \t Index: %d', ...
%			angle,diffex,proj,indx));
