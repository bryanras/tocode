function newtor=alen(ttor)

% This function re-distributes a periodic curve according to arclength.

% ttor should be an nx2 matrix (vector of R2 points).
% Get the size of the curve.
pts=length(ttor);

% Get a list of distances.  dis(i) = jsumtoi(| x(j+1) - x(j) |).
dis=zeros(pts,1);
dis(1)=norm(ttor(2,:)-ttor(1,:));
for kk=2:pts-1
	dis(kk)=dis(kk-1)+norm(ttor(kk+1,:)-ttor(kk,:));
end
	dis(pts)=dis(pts-1)+norm(ttor(pts,:)-ttor(1,:));

% Get the distance increment between points.
inc=dis(pts)/pts;
dis(pts)=dis(pts)*1.000001;			% Prevent numerical error.

% Leave x(1,:) in the same place and re-distribute everything else.
nn=1;							% Counter along dis.
loc=0;							% Location along arc.
newtor(1,:)=ttor(1,:);			% Fix the first point.

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
	newtor(kk,:)=ttor(nn,:)+t*(ttor(mod(nn,pts)+1,:)-ttor(nn,:));

end


