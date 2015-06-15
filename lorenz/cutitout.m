
% This is just a script for cutting the output of runit down to size.

nth = 1;		% Pick every nth point.

ii=1; jj=0; ffs = []; kjjs=[]; yys=[]; 
sz = size(dvec,1);

for jj=1:length(lams)

	lammy = lams(jj);		% This is the lambda value.

	while (dvec(ii,1)==lammy)

		if mod(jj,nth)==0 | jj==length(lams)  % If we want this lambda,
		
			% Fill up the vectors.
			ffs=[ffs;dvec(ii,1:2)];
			if dvec(ii,3)>0
				kjjs=[kjjs;dvec(ii,[1,3])];
				yys=[yys;dvec(ii,[1,4])];
			end

		end

		ii=ii+1;		% Go to the next row in dvec.
		if ii>sz
			break;
		end

	end

end
		

% Now do the same thing for huge.
huges=[];lammy=[];
for ii=1:length(lams)
	if mod((ii-1),nth)==0 | ii==length(lams)
		huges=[huges,huge(:,(3*ii-2):3*ii)];
		lammy=[lammy;lams(ii)];
	end
end
huges(size(huge,1)+1,:)=huges(1,:);

