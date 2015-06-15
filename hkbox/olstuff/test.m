
% This is a file for testing the positive-definiteness of 
% some of the blocks.

scn=1:ptsj;
A=JJ(scn,scn);
B=JJ(scn,scn+ptsj);

ugly=[];
for ii=1:40

	xx=rand(ptsj,1)-0.5;
	xx=xx/norm(xx);
	ugly(ii,1:4)=xx'*[A*xx,B*xx,(A+B)*xx,(B\A)*xx];

%	yy=rand(ptsj*ptsi,1)-0.5;
%	yy=yy/norm(yy);
%	ugly(ii,5)=yy'*JJ*yy;

end


ugly



