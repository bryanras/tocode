function xx=genx(ptsi,ptsj)

% This function generates a "null" eigenvector.

ar=-1;
pls=0;
for ii=1:ptsi

	ar=-ar;
	for jj=1:ptsj

		pls=pls+1;
		xx(pls,1)=ar;

		ar=-ar;
	end
end
	
xx=xx/norm(xx);
	
