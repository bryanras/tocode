%
% This is the main Newton's Method loop.
%

disp(' ')
disp(' ')
disp('Invariant tori are your friends.')
disp(' ')
disp(' ')

% Define the constants.
alp=[1;1;1];
eps0=1;

disp(sprintf('Using alpha1 = %g, alpha2 = %g, alpha3 = %g, and eps0 = %g\n' ...
				, alp(1), alp(2), alp(3), eps0));

% Load the original guess.
pts(1)=input(' Points in the i-direction: ');
pts(2)=input(' Points in the j-direction: ');
pts(3)=input(' Points in the k-direction: ');
disp(' ');
aps=prod(pts);

% Get the initial and final lambda values.
lam=input('Initial lambda (default=0): ');
if isempty(lam) lam=0; end

finlam = input('Final lambda (default=initial): ');
disp(' ')
if isempty(finlam) finlam=lam; end

if finlam>lam
    dellam = input('Initial lambda step (default=0.05): ');
    disp(' ')
    if isempty(dellam) dellam=0.05; end
else
    dellam=10;
end

% Save it to disk?
sv=input('Save the new solution to disk? (y/n default="n") ','s');
if isempty(sv) sv='n'; end

% Generate the initial torus and normals.
olrn=input('(g)enerate a new torus (default) or (l)oad an old one? ','s');

% If loading an old torus, ...
if strcmp(olrn,'l') | strcmp(olrn,'L') 


	stt=input(' Start position (defaults to last-7): ');
	disp(sprintf('\n Loading the initial torus and normals ...')); 

	% Load it up and set xx.
	load huge; load Nb; load lams;
	hgcl = size(huge,2);

	if isempty(stt) stt=hgcl-7; end
	lmcl = size(lams,2);
	xx= huge(:,stt:stt+7);
	avN = avnmls(pts,Nb);

	disp( sprintf('\b Done.') );

% If generating a new torus, ...
else
	disp(sprintf('\n Generating the initial torus and normals ...'));

	clear huge; clear Nb; clear avN;
	xx=whup(pts,alp);
	Nb=onmls(xx,pts);
	avN=avnmls(pts,Nb);
	hgcl=0; lmcl=0;

	disp(sprintf('\b Done.'));
end

% Initialize the solution.
dim=8; cdm=5; tp=aps*cdm; bksz = tp/pts(1);

% Allocate some memory.
ff=zeros(tp,1);				% Discrete function.
yy=zeros(tp,1);				% Newton adjustment to radials.

% Set some tolerances.
ytol=sqrt(tp)*1e-5;
ftol=sqrt(tp)*1e-5;
yinc=2;
nomas=10;
stepmax=200;
deltol=1e-5;
ycap=2*tp;

% The main loop iterates lambda.
newone=1; lam=lam-dellam; ii=0; 
while (lam<finlam-1.5*deltol) & (ii<=stepmax) & (dellam>deltol)

	dellam = min([dellam, finlam-lam]);		% Don't overshoot.

	if newone==1
		lam=lam+dellam;						% Update lambda.
	
		% Get the normal vectors.
		if ii>0
			disp(sprintf('\n Updating and averaging the normal vectors ...'));
			Nb=nnmls(xx,Nb,pts);
			if strcmp(sv,'y') save Nb Nb; end  % Save, if required.
			avN=avnmls(pts,Nb);
			disp(sprintf('\b Done.\n'));
		end

		newone=0;

	end

	% Get the function and enter the Newton loop.
	disp(sprintf(' Calculating the initial function ...'));
	ff=disfunc(xx,Nb,avN,lam,pts,alp,eps0);
	nf=norm(ff);  
	disp(sprintf('\b Done. Initial norm = %g',nf));

	% This tells us which solver to use. Try GMRES at least once per lambda.
	flag = 0;  

	% Here is the Newton loop there are a bunch of exit conditions.
	jj=0; ny=ytol+ycap/2; nyo=ny; death = 0;
	while (nf>ftol) & (ny>ytol) & (max(death)==0)

		nyo=ny;		% This keeps the expansion from getting out of hand.

		% Calculate the Jacobian.
		if mod(jj,1)==0
			clear JJ; pack;
			
			disp(sprintf(' Generating the Jacobian ...'));
			JJ=njacob(lam,pts,ff,alp,eps0,xx,Nb,avN);
			disp(sprintf('\b Done.'));

		end

		% Here are the solvers with different preconditioners.
		disp(sprintf(' Solving the system ...'));
		tolly=max( [ 1e-3, min([ftol/(5*nf),0.5]) ] );

		if flag == 0
			[yy,flag,relres,iter]= ...
				gmres(JJ,ff,25,tolly,5,@upmlsub,[],[],JJ,pts,1);

			% If something goes wrong ...
			if flag ~= 0
				disp(sprintf( ...
				'\b GMRES flag %d. relres = %g Going with direct soln.\n', ...
				flag, relres));
			end

		end
		
		% If GMRES no longer works, ...
		if flag ~=0 
			relres = -1;
			yy=solit(JJ,ff,pts);
			flag = 0;
		end

		% Tell the user what happened.
		yy=-yy; ny=norm(yy);
		disp(sprintf('\b   Residual: %g \t (inf. norm): %g', ...
				 norm(JJ*yy+ff),  norm(JJ*yy+ff,inf)));

		% GMRES/PCG/CGS stuff.
		if relres >= 0
			disp(sprintf( ...
				' Returned outer: %d \tInner: %d \tRelative residual: %g', ...
					iter(1), iter(2), relres));
		end

		% Print out the norms to show us how we're doing.	
		disp(sprintf( ' |ff| =  %g    \t |yy| = %g', nf, ny));

		% Update the solution.
		% I can't figure out a slick way to do this.
		rw=1;
		for pp=1:aps
			nscn=1:dim;
			for ss=1:cdm;
				xx(pp,:)=xx(pp,:)+yy(rw)*Nb(pp,nscn);
				nscn = nscn+dim;
				rw=rw+1;
			end
		end

		ff=disfunc(xx,Nb,avN,lam,pts,alp,eps0);	% Get the new function.
		nf=norm(ff); ny=norm(yy);				% Update the norms.

		jj=jj+1;

		% If things didn't work out, stop the loop.
		death = [ny/nyo>=yinc, jj>=nomas, ny>=ycap, ~isfinite([ny,nf]) ];

	end

	% Show the norms.
	disp(sprintf( ' |ff| =  %g',  nf));

    % Tell the user what is happening.

    % Abject failure.
    if max(death)>0

        disp(sprintf('D''oh! Newton failed at step #: %d Exit: %s',jj, ...
				sprintf('%d', death) ));

        % If it dies on the first try, we're screwed.
        if ii==0; break; end

        % Otherwise, keep going.
        dellam=dellam/2;			% Reset lambda.
        lam=lam-dellam;

		xx=huge(:,hgcl-dim+1:hgcl); 			% Reset xx and nf.
		ff=disfunc(xx,Nb,avN,lam,pts,alp,eps0);	
		nf=norm(ff);

        disp(sprintf('New lambda: %g Step size: %g.\n',lam,dellam));

    % Success!
    else

        % Figure out why we converged.
        if (nf <= ftol) & (ny<=ytol)
            way='F- and y-';
        elseif (nf <= ftol)
            way='F-';
        else
            way='y-';
        end

        % Tell the user. (That's you.)
        disp(sprintf( ...
            'Completed Newton iteration due to %snorm. Lambda = %g.\n'...
            ,way,lam));

		% Save everything and reset the loop.
		ii=ii+1; lmcl=lmcl+1;
        newone=1;
		lams(lmcl)=lam;
		huge(:,hgcl+1:hgcl+dim) = xx; hgcl=hgcl+dim;

		% If we want to save this to disk, go ahead.
		if strcmp(sv,'y') 
			save lams lams alp pts;
			save huge huge;
		end

	end

end

clear death
clear *cl
clear relres
clear resvec
clear steps
clear stt
clear sv

