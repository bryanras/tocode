
multyi = eye(ptsj);
multyii = multyi;

for ii=1:ptsi

	scn=(ii-1)*ptsj+1:ii*ptsj;
	scnp=scn+ptsj;
	if scnp(1)>ptsi*ptsj scnp=scnp-ptsi*ptsj; end

	multyi=(JJ(scn,scnp)\JJ(scn,scn))*multyi;
	multyii=(JJ(scn,scn)\JJ(scn,scnp))*multyii;

end
