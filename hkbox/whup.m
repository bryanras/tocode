function xpy=whup(ptsi,ptsj)

% Get the radii of the torus.
rh=0.95;	% Height of whole thing.
rw=0.95;	% Width of whole thing.
qh=0.35;	% Height of arm cross-section.
qw=0.2;		% Width of arm cross-section.

% Offsets
xoff=0.0; yoff=0.0; zoff=1.0;

% Build the torus.
dphi1=2*pi/ptsi;	dphi2=2*pi/ptsj;

pos=0;
for ii=1:ptsi
		phi1=(ii-1)*dphi1;				% Get phi1.
	for jj=1:ptsj
		phi2=(jj-1)*dphi2;				% Get phi2.
		pos=pos+1;
		xpy(pos,:)=[cos(phi1)*(rw+qw*cos(phi2))+xoff, ...
						-(sin(phi1)*(rh+qw*cos(phi2))+yoff), ...
						qh*sin(phi2)+zoff ];
	end	
end


