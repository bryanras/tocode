function xx=upmlsub(xx,JJ,pts,sgn)

%
% Usage: xx=upmlsub(xx,JJ,pts,sgn)
%
% This is just a function for dividing by U sgn L
%
% sgn should be 1 or -1.
%
% See the notes on 2/12/03.

% Get the necessary stuff.
[tp tmp] = size(JJ);
bksz=round(tp/pts(1));

% Save the last thing.
esec = tp-bksz+1:tp; rowscn=esec;
xxt=xx(esec);

% Just substitute along the blocks.
for ii=pts(1):-1:2

	colscn=rowscn;
	rowscn=rowscn-bksz; 			% Set the sections.

	% Go with any solution technique here.
	% This conditional checks the level of recursion.
	if length(pts)>5
		[xx(colscn,:),ugh] = ...
			gmres(JJ(rowscn,colscn),xx(rowscn,:),25,1e-3,40,@upmlsub, ...
				[],[],JJ(rowscn,colscn),pts(2:3),1);

			if ugh>0
				disp('shite')
			end
	else
%		xx(colscn,:) = JJ(rowscn,colscn)\xx(rowscn,:);
		xx(colscn,:) = solit(JJ(rowscn,colscn),xx(rowscn,:),pts(2:3));
	end

end
	
% Get the final block.
if length(pts)>5
	[xx(rowscn,:),ugh] = gmres(JJ(esec,rowscn),xxt,25,1e-3,40,@upmlsub, ...
				[],[],JJ(esec,rowscn),pts(2:3),1);

	xx(rowscn,:)=sgn*xx(rowscn);

else
%	xx(rowscn,:) = sgn*JJ(esec,rowscn)\xxt;
	xx(rowscn,:) = sgn*solit(JJ(esec,rowscn),xxt,pts(2:3));
end


