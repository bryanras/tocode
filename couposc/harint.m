function [y1, y2] = harint(xx,pts)

% This is an interpolation function that returns the stable
% orbit from the coupled oscillators.
%
% Specifically, it return the intersection of the torus with the
% x1=x3, x2=x4 plane, plus the intersection with the x1=-x3, x2=-x4.
%
% Usage: [y1,y2] = harint(xx,pts)
%
% yi looks like yi = [ phi1, phi2, x1, x2, x3, x4 ]
%
% y1 is the intersection with x1=x3, x2=x4; y2 is the other intersection.
%
% Both vectors are originally sorted according to phi1.

% This is necessary for dealing with problem points.
tolly=5e-3;

% We have to expand out the torus to get the wraparound.
xxt=[];
for ii=1:pts(1)
    scn=(ii-1)*pts(2)+1:ii*pts(2);
    xxt = [xxt; xx( [scn,scn(1)], : ) ];
end
xxt = [xxt; xxt(1:pts(2)+1,:)];

% Find the points where x1=x3 and x1=-x3.
y1=[]; y2=[];
for ii=1:pts(1)

    isc = (ii-1)*(pts(2)+1);
    ipsc= ii*(pts(2)+1);

    for jj=1:pts(2)

        % These are the differences for the x1=x3, x2=x4 plane.
        u1 = xxt(isc+jj,[1:2])  -xxt(isc+jj,[3,4]);
        u2 = xxt(isc+jj+1,[1:2])-xxt(isc+jj+1,[3:4]);
        u3 = xxt(ipsc+jj,[1:2]) -xxt(ipsc+jj,[3:4]);

        % If we hit it going up the torus ...
        if (u1(1)*u2(1) <= 0)  & (u1(2)*u2(2) <= 0)

            % Interpolate linearly. Use an average.
            if max(abs(u1))+max(abs(u2)) < tolly   % If we are lined up, ...
                s=1/2;
            else
                s = mean( [u1(1)/(u1(1)-u2(1)), u1(2)/(u1(2)-u2(2))] );
            end

            y1=[y1; [   ii, (1-s)*jj+s*(jj+1), ...
                        (1-s)*xxt(isc+jj,:)+s*xxt(isc+jj+1,:) ] ];
        end

        % Otherwise, if we hit it going around the torus ...
        if (u1(1)*u3(1) <= 0)  & (u1(2)*u3(2) <= 0)

            % Interpolate linearly. Use an average.
            if max(abs(u1))+max(abs(u3)) < tolly
                s=1/2;
            else
                s = mean( [u1(1)/(u1(1)-u3(1)), u1(2)/(u1(2)-u3(2))] );
            end

            y1=[y1; [   (1-s)*ii+s*(ii+1), jj, ...
                        (1-s)*xxt(isc+jj,:)+s*xxt(ipsc+jj,:) ] ];
        end


        % These are the differences for the x1=-x3, x2=-x4 plane.
        u1 = xxt(isc+jj,[1:2])  +xxt(isc+jj,[3,4]);
        u2 = xxt(isc+jj+1,[1:2])+xxt(isc+jj+1,[3:4]);
        u3 = xxt(ipsc+jj,[1:2]) +xxt(ipsc+jj,[3:4]);

        % If we hit it going up the torus ...
        if (u1(1)*u2(1) <= 0)  & (u1(2)*u2(2) <= 0)

            % Interpolate linearly. Use an average.
            if max(abs(u1))+max(abs(u2)) < tolly   % If we are lined up, ...
                s=1/2;
            else
                %s = mean( [u1(1)/(u1(1)-u2(1)), u1(2)/(u1(2)-u2(2))] );
                s = u1(1)/(u1(1)-u2(1));
            end

            y2=[y2; [   ii, (1-s)*jj+s*(jj+1), ...
                        (1-s)*xxt(isc+jj,:)+s*xxt(isc+jj+1,:) ] ];
        end
           
        % Otherwise, if we hit it going around the torus ...
        if (u1(1)*u3(1) <= 0)  & (u1(2)*u3(2) <= 0)

            % Interpolate linearly. Use an average.
            if max(abs(u1))+max(abs(u3)) < tolly
                s=1/2;
            else
                %s = mean( [u1(1)/(u1(1)-u3(1)), u1(2)/(u1(2)-u3(2))] );
                s = u1(1)/(u1(1)-u3(1));
            end

            y2=[y2; [   (1-s)*ii+s*(ii+1), jj, ...
                        (1-s)*xxt(isc+jj,:)+s*xxt(ipsc+jj,:) ] ];
        end

	end
end


% Sort it all out.
y1=sortrows(y1,1);
y2=sortrows(y2,1);

