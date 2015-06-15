function yy = lastel(JJ,ff,pts)

% This is a function for solving a periodic, block bi-diagonal system, JJx = q.
% See the notes on 2/3/03, p.2.

[tp, jnk] = size(JJ);
yy=zeros(tp,1);

bksz=round(tp/pts(1)); 			% bksz = size of the blocks.
esec=tp-bksz+1:tp; 

% Get the last block of yy.  Go as many terms as necessary.
sn=-1; tolly=tp*1e-5; cor=tolly*2; tms=1;
while (tms<=2) & (norm(cor)>tolly)

	sn=-sn;			% Update sign of correction.

	Nb = tp-(tms+1)*bksz+1:tp-tms*bksz;		% Get the starting block.
	cor = ff(Nb);

	for jj=tms:-1:2
		
		Na=Nb+bksz;
		cor=JJ(Na,Na)*(JJ(Nb,Na)\cor);
		Nb=Na;

	end

	% Add it in.
	yy(esec)=yy(esec)+sn*cor;

	% Put in more terms.
	tms=tms+1;

end

% End it.
yy(esec)=JJ(esec-bksz,esec)\yy(esec);

% Now go with the only possible back substitution that 
% uses all the information.
Na=1:bksz; Nb=Na;
yy(Na) = JJ(esec,Na)\(ff(esec)-JJ(esec,esec)*yy(esec));

for ii=2:pts(1)-2
	Na=Nb; Nb=Nb+bksz;
	yy(Na) = JJ(Na,Nb)\(ff(Na) - JJ(Na,Na)*yy(Na));
end

