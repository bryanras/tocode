function [huge,dvec] = runit(pts,lams,ww)

%
% Usage: 	huge = runit(pts,lams,omega);
%			[huge,dvec] = runit(pts,lams,omega);
%
% This is the main Newton's Method loop.
%
% pts =[ptsi,ptsj] is the number of points.
% lams is the vector of lambdas.
% omega is the main parameter.
%
% Optional output, dvec, is diary vector:
%   [lambda, ||F||, ||y|| , k(Jacobian)]
%
% At the end of each Newton convergence, the diary vector
% tacks on a labmda, and a final ||F||. Everything else
% receives a value of -1. 
%
%

disp(' ')
disp(' ')
disp('Invariant tori are your friends.')
disp(' ')
disp(' ')

% Generic parameters. We only let the user set one.
ofst = 5;
alp = 0.4;
pars = [0,alp,ofst];

% Set the non-continuation parameters.
if nargin<3 || isempty(ww)
	ww = sqrt(0.84);
end
disp(sprintf(...
	'Parameters are omega = %g\t alpha = %g\t offst = %g.\n\n', ...
	ww,alp,ofst))
pars(1) = ww;

% Flag whether the user wants a diary. This is better
% than repeatedly checking nargout.
dflg = nargout>1;

% Initialize the solution.
rr=zeros(pts(1)*pts(2),1);

% Number of continuation steps.
llm = length(lams);

% Initialize the solution.
dvec= [];					% Diary vector: [lam, nf, k(JJ), ny] 
dim =3;
huge=zeros(pts(1)*pts(2),dim*llm);		% The answer.
xx=whup(pts(1),pts(2),ofst);			% The initial torus.
nb = nnmls(xx,pts(1),pts(2));

% Set some tolerances.
ftol=1e-6;
ytol=ftol;;
nomas=15;

% Here is the main loop.
for ii=1:llm

	lam=lams(ii);
	rr(:)=0;

	% Newton's Method parameters
	jj=0; 									% Initialize.
	FF=dfunc(xx,pts(1),pts(2),lam,pars); 	% Get the function.
	nf=norm(FF); 							% Get the norm.
	ny=ytol+1;
	nyo=ny;
	xxt=xx;

	% Here's the Newton loop. It can terminate four ways:
	%	two good, two bad.
	while (nf>ftol) && (ny>ytol) && (jj<nomas)

		% Get the new stuff. We are actually solving for rr.
		% Fix the Jacobian for some iterations, if desired.
		if mod(jj,1)==0
			disp(' Calculating Jacobian ...');
			JJ=njacob(xxt,nb,pts(1),pts(2),lam,pars);	
			disp(sprintf('\b Done.'));
		end

		disp(' Solving system ...');
		yy=-(JJ\FF);
		disp(sprintf('\b Done.'));

		% Get norms.
		ny=norm(yy);

		% Only compute condition numbers if the user is 
		% asking for it.
		if dflg

			% User information.
			kJJ = condest(JJ);
			disp(sprintf( ...
				' |FF| =  %g  \tk(JJ) = %g  \t |yy| = %g',nf, kJJ, ny));

			% Record our progress.
			dvec=[dvec; [lam, nf, ny, kJJ]]; 

		else

			% Same thing, no condition number.
			disp(sprintf(' |FF| =  %g  \t|yy| = %g',nf, ny));
			
		end
		
		
		rr=rr+yy;			% Update the parameters.

		% Get the new function.
		xxt=xx+[rr,rr,rr].*nb;
		FF=dfunc(xxt,pts(1),pts(2),lam,pars); 	
		nf=norm(FF); 			

		jj=jj+1;

	end

	% Show the norms and store everything if the user wants.
	disp(sprintf( ' |FF| =  %g',  nf ));
	if dflg
		dvec=[dvec; [lam, nf, -1, -1]]; 
	end

	% Save the torus in the huge matrix.
	huge(:,dim*ii-dim+1:dim*ii)=xx;

	% Tell the user what is happening.
	% Upon failure, exit ignominiously.
	if (jj >= nomas)
		disp(sprintf( ...
			'Bad way to go: We maxed out Newton iter. Step #: %d', jj));
		huge = huge(:,1:3*ii);
		break;

	elseif isnan(ny) || isnan(nf)
		disp(sprintf( ...
			'Bad way to go: Something weird happened. Step #: %d', jj));
		huge = huge(:,1:3*ii);
		break;
		
	% Otherwise, keep going, or, if we have reached the end, stop.
	elseif  (nf <= ftol)
		disp(sprintf( ...
			'Completed Newton iteration due to F-norm. Lambda = %g', lam));

	else 
		disp(sprintf( ...
			'Completed Newton iteration due to y-norm. Lambda = %g', lam));
	end

	xx=xxt;						% Reset x.

	% Re-distribute if desired.
	if (lam*Inf < 0.35+10*eps) && (ii < llm)
		disp(' Re-distributing torus ...');
		xx=tordis(xx,pts(1),pts(2));
		disp(sprintf('\b Done.'));
	end

end

