function newtor=tordis(xx,pts)


% This function needs work.

% This function re-distributes a torus by applying arc-lenngth
% re-parameterization repeatedly.  It is by no means perfect.

newtor=zeros(size(xx));	% Memory.

ibk = pts(2)*pts(3); jbk=pts(3);

% Re-distribute sections of constant i and j.
for ii=1:pts(1) 
	ioff = (ii-1)*ibk;

	for jj=1:pts(2)
		joff = (jj-1)*jbk;
		scn=ioff+joff+1:ioff+joff+pts(3); 	
		newtor(scn,:)=alen(xx(scn,:));
	end

end

% Re-distribute sections of constant i and k.
for ii=1:pts(1) 
	ioff = (ii-1)*ibk;

	for kk=1:pts(3)
		koff = (kk-1);
		scn=ioff+koff+1:jbk:ioff+koff+jbk*pts(2); 	
		newtor(scn,:)=alen(xx(scn,:));
	end

end


% Re-distribute sections of constant j and k.
for jj=1:pts(2) 
	joff = (jj-1)*jbk;

	for kk=1:pts(3)
		koff = (kk-1);
		scn=joff+koff+1:ibk:joff+koff+ibk*pts(1); 	
		newtor(scn,:)=alen(xx(scn,:));
	end

end

function newtor=alen(tor)

% This function re-distributes a periodic curve according to arclength.
% tor should be an nxm matrix (n-length vector of R^m points).

% Get the size of the curve.
tp=size(tor,1);

% Get a list of distances.  dis(i) = jsumtoi(| x(j+1) - x(j) |).
dis(1)=norm(tor(2,:)-tor(1,:));
for kk=2:tp-1
	dis(kk)=dis(kk-1)+norm(tor(kk+1,:)-tor(kk,:));
end
	dis(tp)=dis(tp-1)+norm(tor(tp,:)-tor(1,:));

% Get the distance increment between points.
inc=dis(tp)/tp;
dis(tp)=dis(tp)*1.00000001;			% Prevent numerical error.

% Leave x(1,:) in the same place and re-distribute everything else.
nn=1;							% Counter along dis.
loc=0;							% Location along arc.
newtor(1,:)=tor(1,:);			% Fix the first point.

% Step through the interpolation loop.
for kk=2:tp

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
	newtor(kk,:)=tor(nn,:)+t*(tor(mod(nn,tp)+1,:)-tor(nn,:));

end


