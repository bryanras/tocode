function x1=nsol(ff,JJ,pts)

% This is the "Newton-x1" solver.  It's almost assuredly a bad idea.
% Check out the notes on 4/14/03 and 4/15/03.

% Block sizes, as usual.
[tp, tmp]=size(JJ);
ibk=tp/pts(1);

% Declare the memory. Note that x1 is actually a matrix full of guesses,
% some of which are perturbed.
global x1t
x1=zeros(ibk,1);
x1t=zeros(ibk,ibk+1);

% Set the tolerance and jiggle amount.
tolly=1e-6; ny=tolly+1; jig=1e-8;
jiggy=jig*speye(ibk);

% Here is the Newton step.
while ny>tolly

	% Calculate the Jacobian. Set up the initial to be processed.
	for ii=1:ibk+1  x1t(:,ii)=x1; end

	% Jiggle it.
	x1t(1:ibk,1:ibk)=x1t(1:ibk,1:ibk)+jiggy;

	% Traipse through the outer Jacobian.
	disp(' Effing sub-J ...')
	effit(ff,JJ,0);
	disp(sprintf('\b  Done.'));

	% Now this is the Jacobian, with the extra column added.
	disp(sprintf('\b  Dividing to get Jacobian ...'))
	for ii=1:ibk
		x1t(:,ii)=(x1t(:,ii)-x1t(:,ibk+1))/jig;
		x1t(ii,ii)=x1t(ii,ii)-1;
	end
	disp(sprintf('\b  Done.'));

	% Get the correction.
	disp(sprintf('\b  Solving ...'))


	% This is just a check.
	multy=eye(ibk,ibk);
	for ii=1:pts(1)
		scna=(ii-1)*ibk+1:ii*ibk;
		scnb=scna+ibk;
		if scnb(1)>tp scnb=scnb-tp; end

		multy=JJ(scna,scnb)\(JJ(scna,scna)*multy);
	end
	
	ugly=x1t(1:ibk,1:ibk)+multy;
	ugly(1:50,1:50)
	norm(multy,inf) 
	norm(ugly,inf)
	norm(x1t(1:ibk,1:ibk),inf)

	yy=x1t(1:ibk,1:ibk)\(x1t(:,ibk+1)-x1); ny=norm(yy);
	disp(sprintf('\b  Done.'));

	% Add it up.
	x1=x1-yy;

end

% That's it.
x1t=x1; effit(ff,JJ,1); x1=x1t;

clear x1t
