function yy=sheeh(ff,JJ,pts)

% This is an experiment with using J^TJ instead of J.
% This uses a modified form of the bordering algorithm.
[tp,tmp]=size(JJ);

% Get the stats.
cdm = tp/prod(pts);
ibk = pts(2)*pts(3)*cdm;

% Declare the memory.
yy=zeros(tp,1);

% Make this global to save memory.
global ctil;
global fbk;

% Start and end sections.
ssec=1:ibk;
esec=(tp-ibk+1):tp;
esecm=esec-ibk;
dsec=1:(tp-2*ibk);
csec=1:(tp-ibk);

% Get the new stuff.   This adds memory requirements, but we only have
% to store half (but not really even that).
ff=JJ'*ff;

% cc and dd are as in the notes from 3/19/03.
disp(sprintf(' Generating Cs and Ds ...'));
[cc,dd]=csnds(JJ,pts);
disp(sprintf('\b Done.'));

% Define b.
bb=sparse(tp-ibk,ibk);
bb(ssec,:)=dd(esec,:)';
bb(esecm,:)=dd(esecm,:);

% For now, we are saving two full blocks.
disp(' Sweeping up ...')
ftil=sweep([ff(csec,:),bb],cc(csec,:),dd(dsec,:),'u');
disp(sprintf('\b Done.'))

% Modify g and calculate one block of bTM^-1b.
ff(esec,:)=ff(esec,:)-dd(esec,:)*( ctil\ftil(:,1) );
fbk=cc(esec,:)-dd(esec,:)*( ctil\ftil(:,2:ibk+1) );

disp(' Sweeping down ...')
ftil=sweep([ff(csec,:),bb],cc(csec,:),dd(dsec,:),'d');
disp(sprintf('\b Done.'))

% Modify g and calculate the other block of bTM^-1b.
ff(esec,:)=ff(esec,:)-dd(esecm,:)'*( ctil\ftil(:,1) );
fbk=fbk-dd(esecm,:)'*( ctil\ftil(:,2:ibk+1) );

yy(esec)=fbk\ff(esec,:);

% We do eventually have to solve for M\ff in its entirety,
disp(' Solving ...')

% Get the rhs and set the tolerance. 
xx=ff(csec)-bb*yy(esec);
tolly=max([5e-4,5e-8/(5*norm(xx))]);

% Solve.
yy(csec)=pcg(@Mmult,xx,tolly,50,@dsub,[],[],cc,dd);



%
% ------------------------------------------------------------------
%

%
% This is a function for doing the full block tri-diagonal substitution.
% This really does take too much time.
%

function xx=msub(ff,cc,dd)


disp(' Sweeping repeatedly ...');

global ctil;

% Get the block sizes.
[tp , ibk]=size(cc);
bks=tp/ibk;
xx=zeros(size(ff));

% Do the last section.
ftil= sweep(ff,cc,dd,'d');

sec=tp-ibk+1:tp;
xx(tp-ibk+1:tp,:)=ctil\ftil;

% We need to sweep bks-1 more times.
for ii=bks-1:-1:1

	% Set the sections.
	secp=sec; sec=sec-ibk; csec=1:ii*ibk; dsec=1:(ii-1)*ibk;

	% Sweep down.
	ftil=sweep(ff(csec,:), cc(csec,:) , dd(dsec,:), 'd' );

	% Get the next section of the solution vector.
	xx(sec,:) = ctil\(ftil - dd(sec,:)*xx(secp,:));
	
end

disp(sprintf('\b Done.'));


%
% ------------------------------------------------------------------
%
	
% This is a function for doing a single up- or downsweep.
function ftil = sweep(ff, cc, dd, swp)

%
% If swp='d', then sweep down and return the last  block element of M^-1ff;
%
% If swp='u', then sweep up and return the first  block element of M^-1ff;
% 
% Remember that cc is the diagonal block, while dd is the off-diagonal.
% There should be a correct number of elements in each.

% Get the sizes.
[tp, ibk]=size(cc);
bks=tp/ibk;

% Define ctil globally to save memory.
global ctil;

if strcmp( swp, 'd' )
	% Declare C~ and f~.
	sec=1:ibk;
	ctil=full(cc(sec,:));	ftil=full(ff(sec,:));
	
	for ii=2:bks

		secm=sec;	sec=sec+ibk; 	% Update the sections.

		% We'd might as well invert.
		ctil=dd(secm,:)'/ctil;

		% Adjust the rhs and diagonal.
		ftil=ff(sec,:)-ctil*ftil;	
		ctil=cc(sec,:)-ctil*dd(secm,:);

	end

% Otherwise, sweep up.
else

	sec=tp-ibk+1:tp;

	% Kick it off.
	ctil=full(cc(sec,:));	ftil=full(ff(sec,:));

	for ii=bks-1:-1:1

		secp=sec;	sec=sec-ibk; 	% Update the sections.

		% We'd might as well invert.
		ctil=dd(sec,:)/ctil;

		% Adjust the rhs and diagonal.
		ftil=ff(sec,:)-ctil*ftil;
		ctil=cc(sec,:)-ctil*dd(sec,:)';

	end

end


%
% ------------------------------------------------------------------
%
	
% This is a function for multiplying Mx, where M is the block tri-diagonal
% part in the bordering algorithm.
function yy = Mmult(xx, cc, dd)

% First, get the number of blocks.  You'd better pass in xx, c's and d's 
% of the correct size. The length of xx determines the overall size.
[tp, ibk]=size(cc);
[tp, tmp]=size(xx);
bks=tp/ibk;					% Should be N_i-1.

% Just loop down. The ends are a little different.
yy=zeros(tp,tmp);
sec = 1:ibk; secp = sec+ibk;
yy(sec,:) = cc(sec,:)*xx(sec,:) + dd(sec,:)*xx(secp,:);

for ii=2:bks-1

	% Update the sections.
	secm=sec; sec=secp; secp=secp+ibk;

	% Multiply.
	yy(sec,:) = dd(secm,:)'*xx(secm,:)+ cc(sec,:)*xx(sec,:) + ...
					dd(sec,:)*xx(secp,:);

end

% Take care of the last block.
yy(secp,:) = dd(sec,:)'*xx(sec,:) + cc(secp,:)*xx(secp,:);


%
% ------------------------------------------------------------------
%
	
% This is a pre-conditioning function for the M-solution.  This preconditioner
% is simply the diagonal.
% The extra inputs are necessary for use in the pcg function.

function xx = dsub(xx, cc, dd)

% Block sizes, as before.
[tp, ibk]=size(cc); [tp, tmp]=size(xx);
bks=tp/ibk;									% Should be N_i-1.

for ii=1:bks

	scn=(ii-1)*ibk+1:ii*ibk;
	xx(scn,:)=cc(scn,:)\xx(scn,:);

end
	
	
