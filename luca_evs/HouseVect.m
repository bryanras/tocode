function w=HouseVect(x)

% w is computed suct that H=I-w*w' is a unitary matrix satisfying (H*x)(2:n)=0

a=abs(x(1)); nx=norm(x); w=x;
if a>0 , w(1)=w(1)*(1+nx/a); w=w/sqrt((a+nx)*nx);
   else w(1)=nx; w=w/nx; end