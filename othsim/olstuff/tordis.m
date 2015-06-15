function newtor=tordis(tor,pts)


% This function needs work.

% This function re-distributes a torus by applying arc-lenngth
% re-parameterization repeatedly.  It is by no means perfect.

% Re-distribute along cross-sections.
for ii=1:pts(1) 
	for jj=1:pts(3)
	newtor((ii-1)*pts(2)+1:ii*pts(2),:)=alen(tor((ii-1)*pts(2)+1:ii*pts(2),:));
end


function newtor=alen(tor)

% This function re-distributes a periodic curve according to arclength.
% tor should be an nxm matrix (n-length vector of R^m points).

% Get the size of the curve.
sz=size(tor);
pts=sz(1);
dim=sz(2);

% Get a list of distances.  dis(i) = jsumtoi(| x(j+1) - x(j) |).
dis(1)=norm(tor(2,:)-tor(1,:));
for kk=2:pts-1
	dis(kk)=dis(kk-1)+norm(tor(kk+1,:)-tor(kk,:));
end
	dis(pts)=dis(pts-1)+norm(tor(pts,:)-tor(1,:));

% Get the distance increment between points.
inc=dis(pts)/pts;
dis(pts)=dis(pts)*1.00000001;			% Prevent numerical error.

% Leave x(1,:) in the same place and re-distribute everything else.
nn=1;							% Counter along dis.
loc=0;							% Location along arc.
newtor(1,:)=tor(1,:);			% Fix the first point.

% Step through the interpolation loop.
for kk=2:pts

	% Figure out where we are along the road.
	loc=loc+inc;
	while (dis(nn)<loc)
		nn=nn+1;
	end

	% Interpolate.  Check to see if we are in the first region.
	if (nn==1)
		t = loc/dis(1);
	else
		t = (loc-dis(nn-1))/(dis(nn)-dis(nn-1));
	end

	% Here is the big line.
	newtor(kk,:)=tor(nn,:)+t*(tor(mod(nn,pts)+1,:)-tor(nn,:));

end


