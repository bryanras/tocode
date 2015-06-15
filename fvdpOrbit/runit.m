function [huge,dvec] = runit(N,lams,ww)

%
% Usage: 	[huge,dvec] = runit(N,lams,omega);
%
% This is the main Newton's Method loop.
%
% N is the number of points.
% lams is the vector of lambdas.
% omega is the main parameter.
%
% Optional output, dvec, is diary vector:
%   [lambda, ||F||, ||y|| , k(Jacobian), period]
%
% At the end of each Newton convergence, the diary vector
% tacks on a labmda, ||F||, and period. Everything else
% receives a value of -1. 
%
% During the Newton process, the routine does not calculate periods,
% so the "period" column receives a -1.
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
	disp(sprintf(...
		'Parameters are omega = %g\t alpha = %g\t offst = %g.\n\n', ...
		ww,alp,ofst))
end
pars(1) = ww;

% Flag whether the user wants a diary. This is better
% than repeatedly checking nargout.
dflg = nargout>1;
	
% Number of continuation steps.
llm = length(lams);

% Initialize the solution.
% Note that we have to grow the diary vector dynamically
% because  we don't know how many iterations we will need.
dvec= [];				% Diary vector: [lam, nf, k(JJ), ny, per] 
dim=3;
huge=zeros(N,dim*llm);	% The answer.
[xx,nb]=whup(N);				

% Set some tolerances.
ytol=1e-4;
ftol=1e-4;
nomas=20;

% Here is the main loop.
for ii=1:llm

	% Get the normal vectors.
	if  ii>1 
		nb=normals(xx,nb);
	end 	

	lam=lams(ii); 						% Update lambda.

	% Newton's Method
	jj=0; ny=ytol+1;
	FF=disfunc(xx,nb,lam,pars); 	nf=norm(FF);  % Get the function.
	
	% There are several reasons why we might break out of here.
	while (jj<nomas) && (nf>ftol) && (ny>ytol) && ~isnan(nf) && ~isnan(ny)

		% Get the new stuff. We are actually solving for rr.
		JJ=njacob(xx,nb,lam,pars);	
		yy=-(JJ\FF);			% Solve.
		ny=norm(yy);

		% Only compute condition numbers if the user is 
		% asking for it.
		if dflg

			% User information.
			kJJ = condest(JJ);
			disp(sprintf( ...
				' |FF| =  %g  \tk(JJ) = %g  \t |yy| = %g',nf, kJJ, ny));

			% Record our progress.
			dvec=[dvec; [lam, nf, ny, kJJ, -1]]; 

		else

			% Same thing, no condition number.
			disp(sprintf(' |FF| =  %g  \t|yy| = %g',nf, ny));
		end
	

		% Update the solution.
		y1=yy(1:2:2*N-1);   y2=yy(2:2:2*N); 
		xx=xx + [y1,y1,y1].*nb(:,1:3) + [y2,y2,y2].*nb(:,4:6);

		% Get the new function and update the norms.
		FF=disfunc(xx,nb,lam,pars); 
		nf=norm(FF); ny=norm(yy);

		jj=jj+1;
		

	end

	% Estimate the total period.
	per = 0;
	for kk=1:N
		kp = mod(kk,N)+1; 
		Phih = norm(funcy(mean(xx([kk,kp],:)),lam,pars));
		sh = norm(xx(kp,:)-xx(kk,:));
		per=per+sh/Phih;
	end
	
	% Show the norms and store everything if the user wants.
	disp(sprintf( ' |FF| =  %g \t Estimated period: %g',  nf, per));
	if dflg
		dvec=[dvec; [lam, nf, -1, -1, per]]; 
	end

	% Save the torus in the huge matrix.
	huge(:,dim*ii-dim+1:dim*ii)=xx;
	
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
		
	elseif  (nf <= ftol)
		disp(sprintf( ...
			'Completed Newton iteration due to F-norm. Lambda = %g', lam));

	else 
		disp(sprintf( ...
			'Completed Newton iteration due to y-norm. Lambda = %g', lam));
	end

	% Re-distribute (unless we're at the end, of course.
	if ii<Inf
		disp(' Re-distributing according to arc length ...')
		xx=alen(xx);
	end

end

