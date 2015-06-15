function yy = solit(JJ,ff,pts)

%
% Usage: yy = solit(JJ,ff,pts)
%
% This is a function for solving a periodic, block bi-diagonal system, JJyy=ff.
%

[tp, tmp] = size(ff);
yy=zeros(tp,tmp);

reccut = 2;						% Recursion cutoff.  Use reccut>3 for none.
ibk=round(tp/pts(1)); 			% ibk size of the blocks.

% First, get y1, the first row block of the solution.
% We have to build a multiplier and keep track of changes in f.
sec=tp-ibk+1:tp; secp=1:ibk;

% Initialize.
multy=-eye(ibk); ft=zeros(ibk,tmp);

% If we are dealing with the bottom recursion level ...
if max(size(pts))>reccut
	disp(sprintf('\b  Dividing block:   '));
end
for ii=1:pts(1)

	if max(size(pts))>reccut
		disp(sprintf('\b\b\b\b\b%3d.',ii));
	end
		
	% Use this only if you want a recursive solution.
	if max(size(pts))>reccut

		ft=ft-multy*solit(JJ(sec,secp),ff(sec,:),pts(2:3));
		multy=-multy*solit(JJ(sec,secp),JJ(sec,sec),pts(2:3));

	else

		multy=-multy/full(JJ(sec,secp));
		ft=ft+multy*ff(sec,:);
		multy=multy*full(JJ(sec,sec));

	end

	% Update the sections.
	secp=sec; sec=sec-ibk;

end


% Get x1.
multy=speye(ibk)+multy;
sec=1:ibk; secp=sec+ibk;
yy(sec,:)=multy\ft;

if max(size(pts))>2
	disp(sprintf('Norm of weird matrix: %g \tCondition #: %g', ...
			norm(multy,inf),condest(multy)));
end

% Substitute forward.
for ii=2:pts(1)

	tempy = ff(sec,:)-full(JJ(sec,sec))*yy(sec,:);
	yy(secp,:) = full(JJ(sec,secp))\tempy;

	sec=secp; secp=secp+ibk;  		% Update the sections.
end

