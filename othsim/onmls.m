function N=onmls(xx,pts)

% This is a fucntion for generating the original normal vectors.

% Get the overall dimensions.
[tp, dim]=size(xx);

% Initialize the solution. Carry the extra zeros.
N=zeros(tp,dim*5);

% Set the cut-offs.
n1=1:dim; n2=dim+1:2*dim; n3=2*dim+1:3*dim; 
n4=3*dim+1:4*dim; n5=4*dim+1:5*dim;

% In this case, we happen to know the normal directions.
pos=0;
for ii = 1:pts(1)
		
	% Keep some of the thing that do not change.
	phi1 = (ii-1)*2*pi/pts(1);							% Angle.
	c1=cos(phi1); s1=sin(phi1);							% Trig stuff.

	% Various and sundry.
	mu1=(c1-s1)/6; nu1 = (c1+s1)/6; mu12=mu1^2; mn1=mu1*nu1;

	ivec = [c1, s1, 0, 0, 0, 0, 0, 0];					% Vector.

	for jj=1:pts(2)

		% Ditto last comment.
		phi2 = (jj-1)*2*pi/pts(2);							% Angle.
		c2=cos(phi2); s2=sin(phi2);							% Trig stuff.

		% Various and sundry.
		mu2=(c2-s2)/6; nu2=(c2+s2)/6; mu22=mu2^2; mn2=mu2*nu2;

		jvec = [0, 0, c2, s2, 0, 0, 0, 0];					% Vector.

		for kk=1:pts(3);
			phi3 = (kk-1)*2*pi/pts(3);						% Angle.
			c3=cos(phi3); s3=sin(phi3);						% Trig stuff.

			% Various & sundry.
			mu3=(c3-s3)/6; nu3=(c3+s3)/6; mu32=mu3^2; mn3=mu3*nu3;
			xi = (mn1+mn2+mn3)/(1+mu12+mu22+mu32);

			kvec = [0, 0, 0, 0, c3, s3, 0, 0];				% Vector.

			% Update the position.
			pos=pos+1;

			% Set the normal vectors.
			N(pos,n1) = ivec; N(pos,n2) = jvec; N(pos,n3) = kvec;

			% These do involve some extra calculations.
			% Maybe I can repair that later.
			N(pos,n4) = [ s1*mu1, -c1*mu1, s2*mu2, -c2*mu2, ...
							s3*mu3, -c3*mu3, 1, 0];

			N(pos,n5) = [	(nu1-mu1*xi)*s1, (mu1*xi-nu1)*c1, ...
							(nu2-mu2*xi)*s2, (mu2*xi-nu2)*c2, ...
							(nu3-mu3*xi)*s3, (mu3*xi-nu3)*c3, -xi, 1];

			% Normalize the last two vectors.
			N(pos,n4) = N(pos,n4)/norm(N(pos,n4));
			N(pos,n5) = N(pos,n5)/norm(N(pos,n5));

		end
	end
end
