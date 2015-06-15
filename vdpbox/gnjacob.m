function JJ=gnjacob(xx,nb,lam,typ)

%
% Usage: JJ=gnjacob(xx,nb,lam,typ)
%
% This is a numerical computation of the Jacobian of the 
% "g-function", as developed on and around 7/28/03. 
%
% Set typ = 1 for Cartesian, anything else for local transformed.

pts=size(nb,1);
JJ=spalloc(2*pts,2*pts,8*pts);

% Get the amount of jiggle needed.
jig=1.0e-9;

% Get the baseline function.
fd=gdfunc(xx,lam,typ);

% Just do it column-wise for now.
for ii=1:pts

	% Get the indices.
	is = 2*ii-1;
	isp = 2*ii;

	% Set which type of jiggle vector we want.
	if typ == 1
		jgvec = [1,0];
	else
		jgvec = nb(ii,:);
	end

	% Jiggle and calculate.
	xx(ii,:)=xx(ii,:)+jig*jgvec;
	JJ(:,is)=(gdfunc(xx,lam,typ)-fd)/jig;
	xx(ii,:)=xx(ii,:)-jig*jgvec;

	xx(ii,:)=xx(ii,:)+jig*[jgvec(2),-jgvec(1)];
	JJ(:,isp)=(gdfunc(xx,lam,typ)-fd)/jig;
	xx(ii,:)=xx(ii,:)-jig*[jgvec(2),-jgvec(1)];

end


