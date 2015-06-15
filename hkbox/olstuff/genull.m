function zz=genull(ptsi,ptsj)

% This is a function for generating a "null" eigenvector.

zz=zeros(ptsi*ptsj,1);

st = -1;
pos=1;
for ii=1:ptsi

	st = -st;
	for jj=1:ptsj
		zz(pos)=st;
		st=-st;
		pos=pos+1;
	end

end
		
