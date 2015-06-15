function JJ=rnjacob(xx,enns,avenns,fd,pts,alp,lam)

% 
%  Usage: JJ = rnjacob(xx,enns,avenns,pts,fd,alp,lam)
%
% This is a numerical computation of the Jacobian of the discrete function.
% It goes basically block row by block row.
%
%

% Allocate memory.  
[tp tmp] = size(fd);
cdm = round(tp/prod(pts)); 		% Co-dimension.
[tmp dim] = size(xx);

% Here is the memory for the Jacobian.
JJ=spalloc(tp,tp,4*prod(pts)*cdm^2);

% Tell the user what's up.
disp(sprintf( '\b Entries: ii=    '));

% Set the amount of jiggle.
jig=1e-8; 

% Step through. Don't indent because that takes up too much space.
scn=(-cdm+1):0;
pos=0;
for ii=1:pts(1)

	% Get the next levels up and down.
	isc = (ii-1)*pts(2);
	ipsc = mod(ii,pts(1))*pts(2);

	for jj=1:pts(2)

		jp = mod(jj,pts(2))+1;

		% Let the user know where we are
		disp(sprintf('\b\b\b\b\b%3d ', ii));

		% Get the point and section. This saves a little work.
		pos = pos+1; scn = scn+cdm;

		% Set out the ordering of the corners of the box.
		crnr = [ pos, ipsc+jj, isc+jp, ipsc+jp ];

		% The local normal vectors.
		oln = reshape(avenns(pos,:),dim,cdm)';

		% Jiggle each corner.
		for pp=1:4

		jpt = crnr(pp); 		% Jiggled point.
		xxt = xx(jpt,:); 		% Temporary point.
		col = cdm*(jpt-1); 		% The column.

		for nn=1:cdm

			% Jiggle it.
			xx(jpt,:) = xx(jpt,:) + jig*enns(jpt,(nn-1)*dim+1:nn*dim);

			% Fill in the transpose of the Jacobian. (See bottom.)
			ugly=((lcfunc(xx(crnr,:),oln,alp,lam)-fd(scn))/jig)';

			JJ(col+nn,scn)=ugly;

			% Un-jiggle it.
			xx(jpt,:) = xxt;


		end
	end

end
end

% We have to fill in the transpose of JJ because of the way the Matlab 
% uses contiguous memory.
JJ=JJ';

