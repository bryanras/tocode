
%
% This is the main Newton's Method loop.
%


disp(' ')
disp(' ')
disp('Invariant tori are your friends.')
disp(' ')
disp(' ')

% Set some parameters.
b=3.0;c=0.25;d=0.2;dim=3;
cool=sprintf( 'I am using the following parameters:');
disp( cool )
cool=sprintf( 'b = %g \nc = %g\nd = %g',b,c,d);
disp( cool )

% Get our necessary input.
ptsi=input('Points in the i-direction: ');
disp(' ')
ptsj=input('Points in the j-direction: ');
disp(' ')

% NUMERICAL EXPERIMENT
scn=1:ptsj; scnp=scn+ptsj;

lams=input('Enter lambdas: (default = [2.005]): ');
disp(' ')
if isempty(lams) lams=2.005; end

% Initialize the solution.
pts=ptsi*ptsj;
xx=whup(ptsi,ptsj);			% The initial torus.

% Re-distribute if desired.
disp(' Re-distributing torus ...');
xx=tordis(xx,ptsi,ptsj); 
disp(sprintf('\b Done.'));

% Save every "svinc" steps.
svinc=1;
	
% Set some tolerances.
ftol=1e-6;
ytol=ftol;
nomas=10;

% Here is the main loop.
clear huge; huge=[0,0];

for ii=1:length(lams)

	lam=lams(ii); 				% Update lambda.
	N=nnmls(xx,ptsi,ptsj);		% Get the normal vectors.
	rr=zeros(pts,1);			% Dist. in normal direction.

	% Newton's Method

	jj=0; ny=ytol+1;
	xxt=xx+[rr,rr,rr].*N;				% Torus at a Newton iteration.
	FF=dfunc(xxt,N,ptsi,ptsj,b,c,d,lam); 	% Get the function.
	nf2=norm(FF); nf=norm(FF,inf); 		% Get the norm.

	% Here's the Newton loop.
	while (jj<nomas) & (nf>ftol) & (ny>ytol)

		% Get the new stuff. We are actually solving for rr.
		% Fix the jacobian for some iterations.
		if mod(jj,1)==0
			disp(' Calculating Jacobian ...');
			JJ=njacob(xxt,N,ptsi,ptsj,b,c,d,lam);	
			cool=sprintf('\b Done.'); disp(cool);
break;
		end

		% NUMERICAL EXPERIMENT
% 		AmB=JJ(scn,scn)-JJ(scn,scnp);
% 		holy = 0;
% 		for kk=1:ptsj
% 			kp=mod(kk,ptsj)+1;
% 			holy= max([ abs((AmB(kk,kk)-AmB(kk,kp))*2/ ...
% 						(AmB(kk,kk)+AmB(kk,kp))),holy]);
% 		end
% 		disp(sprintf(' Condition # : %g \t Weird term : %g', ...
% 				condest(JJ), holy ));

		disp(' Solving system ...');
		yy=-(JJ\FF);
		cool=sprintf('\b Done.'); disp(cool);

		% Get norms.
		ny=norm(yy,inf); ny2=norm(yy);

		% Print out the norms to show us how we're doing.	
		cool=sprintf( ' |FF| =  %g \t |yy| = %g \t |FF|2 = %g \t |yy|2 = %g', nf, ny, nf2, ny2);
		disp( cool )

		rr=rr+yy;			% Update the parameters.

		% Get the new function.
		xxt=xx+[rr,rr,rr].*N;
		FF=dfunc(xxt,N,ptsi,ptsj,b,c,d,lam); 	
		nf=norm(FF,inf); nf2=norm(FF); 				% Update the norms.

		jj=jj+1;

	end


	% Update the torus.
	xx=xxt;

	% Show the norms.
	cool=sprintf( ' |FF| =  %g   |FF|2 = %g ',  nf,nf2);
	disp( cool )

	% Tell the user what is happening.
	if (jj >= nomas)
		disp(sprintf('Bad way to go: We maxed out Newton iter. Step #: %d', jj));
		break;
	elseif  (nf <= ftol)
		cool=sprintf( 'Completed Newton iteration due to F-norm. Lambda = %g', lam);
		disp( cool )
	else 
		cool=sprintf( 'Completed Newton iteration due to y-norm. Lambda = %g', lam);
		disp( cool )
	end

	% Save the torus in the huge matrix.
	% Get our position.  If this is the first time, start at 1.
	temp=size(huge);
	if temp(1)==1
		clear huge;
		temp(2)=0;
	end

	% Always save the last one.
	if ii==length(lams)
		huge(:,temp(2)+1:temp(2)+dim)=xx;
	elseif mod(ii,svinc)==0 
		huge(:,temp(2)+1:temp(2)+dim)=xx;
	end

	% Re-distribute, if desired.
	disp(' Re-distributing torus ...');
	xx=tordis(xx,ptsi,ptsj); 
	disp(sprintf('\b Done.\n\n'));

	
end

% Clean up the mess.
%clear FF;
%clear JJ; 
%clear N; 
%clear b;
%clear c;
clear cool;
clear cut;
%clear d;
clear damn;
%clear dim;
%clear ftol;
%clear ii;
%clear jj;
%clear lam;
%clear nf;
%clear nf2;
%clear ny;
%clear ny2;
%clear pts*;
%clear rr;
%clear steps;
clear svinc;
clear temp;
%clear xx;
%clear xxt;
%clear ytol;
%clear yy;

