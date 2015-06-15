
%
% This is the main Newton's Method loop.
%


disp(' ')
disp(' ')
disp('Invariant tori are your friends.')
disp(' ')
disp(' ')

% Get our necessary input.
pts=input('How many points in the torus? ');
disp(' ')

lams=input('Enter lambdas (default=[0.1]): ');
disp(' ')
if isempty(lams) lams=0.1; end

% Initialize the solution.
xx=whup(pts);				% The "torus".

% Allocate some memory.
huge=zeros(pts,2*length(lams));	% The answer.
FF=zeros(pts,1);					% Discrete function.
yy=zeros(pts,1);					% Newton adjustment to rr.
nb=zeros(pts,2);					% Normal vectors.

% Set some tolerances.
ytol=0.001;
ftol=0.001;
nomas=50;

% Here is the main loop.
for ii=1:length(lams)

	lam=lams(ii); 		% Update lambda.
	nb=normals(xx);		% Get the normal vectors.
	rr=zeros(pts,1);	% Dist. at each pt. in normal direction.

	% Newton's Method

	jj=0; ny=ytol+1;
	FF=disfunc(xx,lam); 	nf=norm(FF);  % Get the function.
	while (jj<nomas) & (nf>ftol) & (ny>ytol)

		% Get the new stuff. We are actually solving for rr.
		JJ=njacob(xx+[rr,rr].*nb,nb,lam);	
		yy=-(JJ\FF);			% Solve.
		ny=norm(yy);

		% Print out the norms to show us how we're doing.	
		disp(sprintf( ' |FF| =  %g  \tk(JJ) = %g  \t |yy| = %g', ...
			nf, condest(JJ), ny));

		rr=rr+yy;								% Update the parameters.

		FF=disfunc(xx+[rr,rr].*nb,lam); 		% Get the new function.
		nf=norm(FF); ny=norm(yy);				% Update the norms.

		jj=jj+1;
		

	end

	% Update the torus.
	xx=xx+[rr,rr].*nb;

	% Show the norms.
	cool=sprintf( ' |FF| =  %g',  nf);
	disp( cool )

	if (jj >= nomas)
		damn=sprintf('Bad way to go: We maxed out Newton iter. Step #: %d', jj);
		disp( damn )
		break;
	elseif  (nf <= ftol)
		cool=sprintf('Completed Newton iteration due to F-norm. Alpha = %g', lam);
		disp( cool )
	else 
		cool=sprintf( 'Completed Newton iteration due to y-norm. Alpha = %g', lam);
		disp( cool )
	end

	% Save the torus in the huge matrix.
	huge(:,2*ii-1:2*ii)=xx;
	
	% Re-distribute.
	if ii<length(lams)
		disp(' Re-distributing according to arc length...')
		xx=alen(xx);
	end

end

