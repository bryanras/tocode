function JJ=cnjacob(lam,pts,fd)

% This is a numerical computation of the Jacobian of the discrete function.
% This algorithm works column-wise.

% Get the global stuff.
global xx enns avenns

% Allocate memory.  
[tp, dim]=size(xx); [ugh, ndim]=size(enns);
cdm = round(ndim/dim); 							% Co-dimension.
JJ=spalloc(5*tp,5*tp,25*8*tp);

nt = 1:ndim;
nt=reshape(nt,dim,cdm)';

% Tell the user what's up.
cool =sprintf( '\b  Generating Entries: ii=     jj=    '); disp(cool)

% Set the amount of jiggle needed.
jig=1.0e-8; 

% Size of blocks.
ibk=pts(2)*pts(3);

% Step through. Don't indent because that takes up too much space.
cpt=0; col=0;
for ii=1:pts(1)

% Get the next levels up and down.
imsc = mod(ii-2,pts(1))*ibk;
iisc = (ii-1)*ibk;
ipsc = mod(ii,pts(1))*ibk;

for jj=1:pts(2)

jmsc = mod(jj-2,pts(2))*pts(3);
jjsc = (jj-1)*pts(3);
jpsc = mod(jj,pts(2))*pts(3);

cool=sprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bii=%3d   jj=%3d', ii,jj);
disp(cool)

for kk=1:pts(3)

	km = mod(kk-2,pts(3))+1;
	kp = mod(kk,pts(3))+1;

	cpt = cpt+1; 				% Bump the point.

	% Set out the ordering of the corners of the box.
	% There are 8 boxes. This is a little inefficient, but not too bad.
	crnr(1,1:4)=[imsc+jmsc+km, iisc+jmsc+km, iisc+jjsc+km, imsc+jjsc+km];
	crnr(1,5:8)=crnr(1,1:4)-km+kk;

	crnr(2,:)= [crnr(1,5:8), crnr(1,5:8)-kk+kp] ;

	crnr(3,1:4)= [imsc+jjsc+km, iisc+jjsc+km, iisc+jpsc+km, imsc+jpsc+km];
	crnr(3,5:8)= crnr(3,1:4)-km+kk;

	crnr(4,:)= [crnr(3,5:8), crnr(3,5:8)-kk+kp];

	crnr(5,1:4)= [iisc+jmsc+km, ipsc+jmsc+km, ipsc+jjsc+km, iisc+jjsc+km];
	crnr(5,5:8)= crnr(5,1:4)-km+kk;

	crnr(6,:)= [crnr(5,5:8), crnr(5,5:8)-kk+kp];

	crnr(7,1:4)= [iisc+jjsc+km, ipsc+jjsc+km, ipsc+jpsc+km, iisc+jpsc+km];
	crnr(7,5:8)= crnr(7,1:4)-km+kk;

	crnr(8,:)= [ crnr(7,5:8), crnr(7,5:8)-kk+kp];


	% Jiggle the point in each normal direction.
	xxt=xx(cpt,:);
	for nn=1:cdm

		col=col+1; 				% Bump up the column.

		% Jiggle.
		xx(cpt,:) = xx(cpt,:) + jig*enns(cpt,nt(nn,:)); 

		% Now see how this affects each of the boxes.
		for pp=1:8


			% Get the section and column affected.
			tmp = 5*(crnr(pp,1)-1); scn=tmp+1:tmp+cdm;

			% Set the local normal vectors.
			oln = reshape(avenns(crnr(pp,1),:),dim,cdm)';

			% Fill in the Jacobian.
			JJ(scn,col)=((lcfunc(xx(crnr(pp,:),:),oln,lam)-fd(scn))/jig);
			
		end

		% Un-jiggle.
		xx(cpt,:) = xxt;

	end


end
end
end

% That's all, folks.
cool = sprintf('\b  Done.');
disp(cool)

