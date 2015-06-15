function effit(ff,JJ,typ)

% This is my latest crazy idea that probably won't work.
% The idea is to guess an x1t, traipse though the Jacobian, then
% see what that would give as a new x1t.

% Save memory.
global x1t

% Block sizes, as usual.
[tp, tmp]=size(JJ);
[ibk, tmp]=size(x1t);
bks=tp/ibk;

% If we want the whole damn thing returned, ...
if typ==1 yy=zeros(tp,tmp); end

% There isn't much to this.
scn=1:ibk; scnb=scn+ibk;
for ii=1:bks

	% Do it in parts to handle differently-sized right-hand sides.
	x1t=JJ(scn,scn)*x1t;
	effy=ff(scn,:);
	for jj=1:tmp
		x1t(:,jj)=effy-x1t(:,jj);
	end

	x1t=JJ(scn,scnb)\x1t;

	% Save it if desired.
	if typ==1 yy(scnb,:)=x1t; end

	% Update the sections.
	scn=scnb; scnb=scnb+ibk;
	if max(scnb)>tp scnb=scnb-tp; end

end

% Again, save the whole thing if desired.
if typ==1 x1t=yy; end

