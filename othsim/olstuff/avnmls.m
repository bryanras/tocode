function avN=avnmls(pts,N)

% This is a function for averaging the normal vectors.
% This is just so we don't have to do the calculation repeatedly.

% This calculation should only occur once:
pjpk=pts(2)*pts(3);

pos=0;
for ii=1:pts(1)
   
    % Get the sections.
    isc = pjpk*(ii-1); ipsc = mod(ii,pts(1))*pjpk;

    for jj=1:pts(2)
   
        % Ditto.
        jsc = pts(3)*(jj-1); jpsc = mod(jj,pts(2))*pts(3);

        for kk=1:pts(3)

            kp = mod(kk,pts(3))+1;

            % Remember the ordering in lcfunc.
            pos=pos+1;
            crnr(1:4) = [pos,ipsc+jsc+kk,ipsc+jpsc+kk,isc+jpsc+kk];
            crnr(5:8) = crnr(1:4)-kk+kp;

            % Average it all.
            avN(pos,:) = sum(N(crnr,:));

        end
    end
end

