
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

% Change this for the initial guess.
deflam=0.6;

disp(sprintf('Default Lambda = [%g].', deflam));
lams=input('Lambda values (RET to accept default): ');
disp(' ')
if isempty(lams) lams=deflam; end
steps=max(size(lams));

% Initialize the solution.
xx=whup(pts);				% The "torus".

% Allocate some memory.
huge=zeros(pts,2*steps);	% The answer.
FF=zeros(pts);					% Discrete function.
yy=zeros(pts);					% Newton adjustment to rr.
nb=zeros(pts,2);				% Normal vectors.

% Set some tolerances.
ytol=1e-6;
ftol=1e-6;
nomas=50;

% Here is the main loop.

for ii=1:steps

	lam=lams(ii);	 	% Update lambda.
	nb=normals(xx);		% Get the normal vectors.
	rr=zeros(pts,1);	% Dist. at each pt. in normal direction.

	% Newton's Method


	jj=0; ny=ytol+1;
	FF=dfunc(xx,lam); 	nf=norm(FF);  % Get the function.
	while (jj<nomas) & (nf>ftol) & (ny>ytol)

		% Get the new stuff. We are actually solving for rr.
		JJ=njacob(xx+[rr,rr].*nb,nb,lam);	
		cool=sprintf( 'Condition #: %g', condest(JJ));
		disp(cool);
		yy=-(JJ\FF);			% Solve.
		ny=norm(yy);

		% Print out the norms to show us how we're doing.	
		cool=sprintf( ' |FF| =  %g    \t |yy| = %g', nf, ny);
		disp( cool )

		rr=rr+yy;							% Update the parameters.

		FF=dfunc(xx+[rr,rr].*nb,lam); 		% Get the new function.
		nf=norm(FF); ny=norm(yy);			% Update the norms.

		jj=jj+1;
		

	end

	% Update the torus.
	xx=xx+[rr,rr].*nb;
%	disp(' Re-distributing according to arc length...')
%	xx=alen(xx);

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

end

%nb=normals(xx);
%JJ=njacob(xx,nb,lam);
