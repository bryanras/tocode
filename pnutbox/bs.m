function yy = bs( ll, u1, u2, pp, yy)

% The function does back-substitution on the output of 
% luit.m
[N N2]=size(u1);

% Solve Lz = Py 
yy=yy(pp);
for ii=2:N
		if pp(ii)~=N
			cc=pp(ii);
		else
			cc=1;
		end
		yy(ii) = yy(ii) - ll(cc:ii-1)'*yy(cc:ii-1);
end

% Solve Ux = z
yy(N) = yy(N)/u1(N);
for ii=N-1:-1:1
	
	% Check for permutation.
	if pp(ii) == ii
		yy(ii) = (yy(ii) - yy(ii+1)*u2(ii))/u1(ii);
	else
		yy(ii) = (yy(ii) - yy(N)*u2(ii))/u1(ii);
	end

end

