function avenns=avnmls(pts,enns)

% This is a function for ``averaging'' (actually, summing) the normal vectors.
% This is just so we don't have to do the calculation repeatedly.

% This calculation should only occur once:

for ii=1:pts(1)
   
    % Get the sections.
    isc = pts(2)*(ii-1); ipsc = mod(ii,pts(1))*pts(2);

    for jj=1:pts(2)
   
        jp = mod(jj,pts(2))+1;

        % Average it all.
		avenns(isc+jj,:) = sum( enns([isc+jj,ipsc+jj,isc+jp,ipsc+jp],:), 1 );

    end
end

