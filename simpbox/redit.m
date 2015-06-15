
% This is a little ditty for saving every kth column of the huge
% matrix that comes out of our algorithm.

k=5;

blip=size(huge);
n=blip(2);

% Always get the first row.
clear stillbig;
stillbig(:,1:2)=huge(:,1:2);

for ii=2*k:2*k:n
	blip=round(ii/k);
	stillbig(:,blip+1:blip+2)=huge(:,ii+1:ii+2);
end

