%
% This is the main Newton's Method loop.
%

disp(' ')
disp(' ')
disp('Invariant tori are your friends.')
disp(' ')
disp(' ')

% Define the variables.
alp=0.55;
cool=sprintf('Using alpha = %g\n', alp);
disp(cool)

% Load the original guess.
pts(1)=input(' Points in the i-direction: ');
pts(2)=input(' Points in the j-direction: ');
disp(' ');
aps=prod(pts);

% Get the initial lambda value
lams=input('Lambda values (defaults to [0]): ');
if isempty(lams) lams=0; end
steps=length(lams);

disp(sprintf('\n Generating the initial torus and normals ...'));
[ xx, enns ] = whup(pts);  % Starting from scratch.

%xx = huge(:,61:64);				% Picking up an old solution.
%enns=nnmls(xx,enns,pts);

avenns = avnmls(pts,enns);
disp(sprintf('\b Done.'));

% Initialize the solution and set some dimensions.
dim=4; cdm=2; tp=aps*cdm; 

% Allocate some memory.
ff=zeros(tp,1);				% Discrete function.
yy=zeros(tp,1);				% Newton adjustment to radials.

% Allocate the memory for one continuation step.
huge=zeros(aps,dim);

% Set some tolerances.
ytol=1e-6;
ftol=1e-6;
nomas=20;

% The main loop iterates lambda.
for ii=1:steps

	% Update lambda;
	lam=lams(ii);

	% Get the function and enter the Newton loop.
	disp(sprintf(' Calculating the initial function ...'));
	ff=disfunc(xx,enns,avenns,pts,alp,lam);
	nf=norm(ff);
	disp(sprintf('\b Done.'));

	% Here is the Newton loop.
	jj=0; ny=ytol+1;
	while (jj<nomas) & (nf>ftol) & (ny>ytol)

		% Calculate the Jacobian.
		if mod(jj,1)==0
			clear JJ;			
			disp(sprintf(' Generating the Jacobian ...'));
			JJ=rnjacob(xx,enns,avenns,ff,pts,alp,lam);
			disp(sprintf('\b Done.'));

		end

		% Here are the solvers with different preconditioners.
		disp(sprintf(' Condition #: %g', condest(JJ)));
		disp(sprintf(' Solving the system ...'));
		yy=-JJ\ff;
		ny=norm(yy);
		disp(sprintf('\b Done.'));

		% Print out the norms to show us how we're doing.	
		disp(sprintf( ' |ff| =  %g    \t |yy| = %g', nf, ny));

		% Update the solution.
		% I can't figure out a slick way to do this.
		rw=1;
		for pp=1:aps
			nscn=1:dim;
			for ss=1:cdm;
				xx(pp,:)=xx(pp,:)+yy(rw)*enns(pp,nscn);
				nscn = nscn+dim;
				rw=rw+1;
			end
		end

		ff=disfunc(xx,enns,avenns,pts,alp,lam);		% Get the new function.
		nf=norm(ff); ny=norm(yy);					% Update the norms.

		jj=jj+1;

	end

	% Show the norms.
	cool=sprintf( ' |ff| =  %g',  nf);
	disp( cool )

	% Tell the user about the exit condition.
	if (jj >= nomas)
		damn=sprintf('Bad way to go: We maxed out Newton iter. Step #: %d', jj);
		disp( damn )
		break;
	elseif  (nf <= ftol)
		cool=sprintf('Completed Newton iteration due to F-norm. Lambda = %g', lam);
		disp( cool )
	else 
		cool=sprintf('Completed Newton iteration due to y-norm. Lambda = %g', lam);
		disp( cool )
	end

	% Save everything in a huge matrix.
	huge(:,(ii-1)*dim+1:ii*dim) = xx;

	% Update the normals.
	if ii<steps
		disp(sprintf('\n\n Updating and averaging the normal vectors ...'));

		if lam < 0.25000001
			xx=tordis(xx(:,[3,4,1,2]),pts(1),pts(2));
			xx=tordis(xx(:,[3,4,1,2]),pts(1),pts(2));
		end
		enns=nnmls(xx,enns,pts);
		avenns=avnmls(pts,enns);
		disp(sprintf('\b Done.\n\n'));
	end

end

clear cool
clear damn
clear relres
clear resvec
clear steps
clear sv
clear tmp
clear tmpii

