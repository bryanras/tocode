function JJ=fixj(JJ,pts)


% This is a function for repairing the Jacobian.
% It's basically like pre-pivoting.
cdm=5;

% Loop through.
sec=-cdm+1:0;
for ii=1:pts(1)
	ipsc=mod(ii,pts(1))*pts(2)*pts(3);

for jj=1:pts(2)
	jpsc=mod(jj,pts(2))*pts(3);

for kk=1:pts(3)
	kp=mod(kk,pts(3))+1;

	% Get the sections.
	sec=sec+cdm;
	ps=ipsc+jpsc+kp;
	psec=(ps-1)*cdm+1:ps*cdm;
	if (max(sec)*cdm > 5*length(JJ(:,1))) | (ps*cdm > 5*length(JJ(:,1)))
		ii
		jj
		kk
	end

	% Flip it.
	JJ([sec,psec],:)=JJ([psec,sec],:);
	
	
end
end
end
