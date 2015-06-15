function cl=centerline(torus,ptsj)

%
% cl = centerline(torus,ptsj)
%
% This is a function for calculating the centerline of a torus.
% It's not that complicated.

[pts,tmp]=size(torus);

for ii=1:floor(pts/ptsj)

	cl(ii,:)=mean(torus((ii-1)*ptsj+1:ii*ptsj,:),1);

end

cl(ii+1,:)=cl(1,:);

