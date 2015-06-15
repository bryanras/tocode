
%
% This is the main Newton's Method loop.
%


disp(' ')
disp(' ')
disp('Invariant tori are your friends.')
disp(' ')
disp(' ')

% Set the non-continuation parameters.
global sig b;
sig = 10;
b = 8/3;

disp(sprintf('Parameters are sigma = %g and b = %g.\n\n',sig,b))


% Get our necessary input.
N=input('How many points in the orbit? ');
disp(' ')

lams=input('Enter lambdas (default=[24.5]): ');
disp(' ')
if isempty(lams) lams=24.5; end

% Initialize the solution.
dvec= [];  			  			% Diary vector: [lam, nf, k(JJ), ny], 
dim=3; cdm=2;
huge=zeros(N,dim*length(lams));	% The answer.
[xx,nb]=whup(N);				% If using the output from some other source.
%load xxnb; 						% If using previously-calculated values.

% Set some tolerances.
ytol=1e-5;
ftol=1e-5;
nomas=50;

% Here is the main loop.
for ii=1:length(lams)

	if  ii>1 nb=normals(xx,nb);  end 	% Get the normal vectors.

	lam=lams(ii); 						% Update lambda.

	% If we want to save at a particular value.
	if lam<=13.93 
		save xxnb huge xx nb dvec
	end

	% Newton's Method
	jj=0; ny=ytol+1;
	FF=disfunc(xx,nb,lam); 	nf=norm(FF);  % Get the function.
	while (jj<nomas) & (nf>ftol) & (ny>ytol)

		% Get the new stuff. We are actually solving for rr.
		JJ=njacob(xx,nb,lam);	
		yy=-(JJ\FF);			% Solve.
		ny=norm(yy);

		% Print out the norms to show us how we're doing.	
		kJJ = condest(JJ);
		disp(sprintf( ' |FF| =  %g  \tk(JJ) = %g  \t |yy| = %g',nf, kJJ, ny));
	
		% Record our progress.
		dvec=[dvec; [lam, nf, kJJ, ny]]; 

		% Update the solution.
		y1=yy(1:2:2*N-1);   y2=yy(2:2:2*N); 
		xx=xx + [y1,y1,y1].*nb(:,1:3) + [y2,y2,y2].*nb(:,4:6);

		FF=disfunc(xx,nb,lam); 			% Get the new function.
		nf=norm(FF); ny=norm(yy);		% Update the norms.

		jj=jj+1;

	end

	% Estimate the total period.
	per = 0;
	for kk=1:N
		kp = mod(kk,N)+1; 
		Phih = norm(funcy(mean(xx([kk,kp],:)),lam));
		sh = norm(xx(kp,:)-xx(kk,:));
		per=per+sh/Phih;
	end
	
	% Show the norms and store everything.
	disp(sprintf( ' |FF| =  %g \t Estimated period: %g',  nf, per));
	dvec=[dvec; [lam, nf, -1, -1]]; 

	if (jj >= nomas)
		disp(sprintf('Bad way to go: We maxed out Newton iter. Step #: %d', jj));
		break;
	elseif  (nf <= ftol)
		disp(sprintf('Completed Newton iteration due to F-norm. Lambda = %g', lam));
	else 
		disp(sprintf( 'Completed Newton iteration due to y-norm. Lambda = %g', lam));
	end

	% Save the torus in the huge matrix.
	huge(:,dim*ii-dim+1:dim*ii)=xx;

	
	% Re-distribute.
	if ii<length(lams)
		disp(' Re-numbering and distributing according to arc length ...')
		[xx,nb]=reord(xx,nb,lam);
		xx=alen(xx);
	end

end



