function T=tanny(xx,pts)

% This is a fucntion for generating tangent vectors.
%
% Usage: T=tanny(xx,pts)
%
%

% Do we want to orthonormalize?  If so, give on an nonzero value.
on = 1;

% Get the overall dimensions.
[tp, dim]=size(xx);

% Initialize.
T=zeros(tp,dim*3);

% In this case, we happen to know the normal directions.
pos=0;
for ii = 1:pts(1)
		
	% Keep some of the thing that do not change.
	phi1 = (ii-1)*2*pi/pts(1);							% Angle.
	c1=cos(phi1); s1=sin(phi1);							% Trig stuff.

	ivec = [-s1, c1, 0, 0, 0, 0, (c1-s1)/6, (c1+s1)/6];	% Vector.

	% Normalize.
	if on ~= 0 
		ivec=ivec/norm(ivec);
	end

	for jj=1:pts(2)

		% Ditto last comment.
		phi2 = (jj-1)*2*pi/pts(2);							% Angle.
		c2=cos(phi2); s2=sin(phi2);							% Trig stuff.

		jvec = [0, 0, -s2, c2, 0, 0, (c2-s2)/6, (s2+c2)/6];	% Vector.

		% If you want to orthonormalize.
		if on ~= 0 
			jvec = jvec-ivec*(jvec*ivec');
			jvec=jvec/norm(jvec);
		end

		for kk=1:pts(3);
			phi3 = (kk-1)*2*pi/pts(3);						% Angle.
			c3=cos(phi3); s3=sin(phi3);						% Trig stuff.

			kvec = [0, 0, 0, 0, -s3, c3, (c3-s3)/6, (s3+c3)/6];	% Vector.

			% Normalize.
			if on ~= 0 
				kvec = kvec - jvec*(kvec*jvec') - ivec*(kvec*ivec');
				kvec=kvec/norm(kvec);
			end


			% Update the position.
			pos=pos+1;

			% Set the normal vectors.
			T(pos,:) = [ivec,jvec,kvec]; 

		end
	end
end
