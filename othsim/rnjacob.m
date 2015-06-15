function JJ=rnjacob(lam,pts,fd,alp,eps0,xx,N,avN)

% This is a numerical computation of the Jacobian of the discrete function.
% It goes basically row-by-row.

% Allocate memory.  
[tp tmp] = size(fd);
cdm = round(tp/prod(pts)); 		% Co-dimension.
[tmp dim] = size(xx);

% Here is the memory for the Jacobian.
JJ=spalloc(tp,tp,8*prod(pts)*cdm^2);

% Tell the user what's up.
cool =sprintf( ' Generating Jacobian Entries: ii=     jj=    '); disp(cool)

% Set the amount of jiggle.
jig=1e-8; 

% Size of blocks.
ibk=pts(2)*pts(3);

% Step through. Don't indent because that takes up too much space.
scn=(-cdm+1):0;
pos=0;
for ii=1:pts(1)

% Get the next levels up and down.
isc = (ii-1)*ibk;
ipsc = mod(ii,pts(1))*ibk;

for jj=1:pts(2)

jsc = (jj-1)*pts(3);
jpsc = mod(jj,pts(2))*pts(3);

cool=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bii=%3d   jj=%3d', ii,jj);
disp(cool)

% Do this as before, box-by-box.
% It is almost as efficient as doing the single-jiggle.
for kk=1:pts(3)

	kp = mod(kk,pts(3))+1;

	% Get the point and section. This saves a little work.
	pos = pos+1; scn = scn+cdm;

	% Set out the ordering of the corners of the box.
	crnr(1:4) = [ pos, ipsc+jsc+kk, ipsc+jpsc+kk, isc+jpsc+kk ];
	crnr(5:8) = crnr(1:4)-kk+kp;

	% The local normal vectors.
	oln = reshape(avN(pos,:),dim,cdm)';

	% Jiggle each corner.
	for pp=1:8

		jpt = crnr(pp); 		% Jiggled point.
		xxt = xx(jpt,:); 		% Temporary point.
		col = cdm*(jpt-1); 		% The column.

		for nn=1:cdm

			% Jiggle it.
			xx(jpt,:) = xx(jpt,:) + jig*N(jpt,(nn-1)*dim+1:nn*dim);

			% Fill in the transpose of the Jacobian. (See bottom.)
			ugly=((lcfunc(xx(crnr,:),oln,lam,alp,eps0)-fd(scn))/jig)';

			% Experiment.
			%ugly=ugly.*mask(nn,:);
	
			JJ(col+nn,scn)=ugly;

			% This is an experiment.
			
			% Un-jiggle it.
			xx(jpt,:) = xxt;

		end

	end


end
end
end

% We have to fill in the transpose of JJ because of the way the Matlab 
% uses contiguous memory.
JJ=JJ';

% That's all, folks.
cool = sprintf('\b  Done.');
disp(cool)

