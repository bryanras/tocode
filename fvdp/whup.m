function xpy=whup(ptsi,ptsj,ofst)

% Get the radii of the torus.
qh=2;		% Radius of arm cross-section.

% Build the torus.
dphi1=2*pi/ptsi;	dphi2=2*pi/ptsj;

pos=0;
for ii=1:ptsi
		phi1=(ii-1)*dphi1;				% Get phi1.
	for jj=1:ptsj
		phi2=(jj-1)*dphi2;					% Get phi2.
		pos=pos+1;
		xpy(pos,:)=[cos(phi1)*(ofst+qh*cos(phi2)), ...
						-sin(phi1)*(ofst+qh*cos(phi2)), ...
						qh*sin(phi2)];
	end	
end


