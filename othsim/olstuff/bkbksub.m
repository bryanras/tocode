function xx=bkbksub(bb,A,bksz)

% This is a function for solving a block-upper-bi-diagonal matrix
% system, Ax=b.

% Get the sizes of everything.

[tp tmp]=size(bb);

xx=zeros(tp,1);						% Declare some memory.

% Take care of the end.
sec1 = tp-bksz+1:tp;
xx(sec1)=A(sec1,sec1)\bb(sec1);

% Now go block by block.
while sec1(1)>bksz

	sec2=sec1; sec1=sec1-bksz;	% Update the sections.

	% Back substitute.
	xx(sec1) = A(sec1,sec1)\(bb(sec1)-A(sec1,sec2)*xx(sec2));
	
end
