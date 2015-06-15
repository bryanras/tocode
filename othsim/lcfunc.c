#include "locf.h"

/* This is the "gateway" routine for the local function calculation. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  	double *xx, *oln, *outy, *alp;
	double lam, eps0;
  
	/* I'm not doing ANY error checking on the inputs. Be careful. */

	/* Assign pointers to input and get the scalar. */
	xx = mxGetPr(prhs[0]);
	oln = mxGetPr(prhs[1]);
	lam = mxGetScalar(prhs[2]);
	alp = mxGetPr(prhs[3]);
	eps0 = mxGetScalar(prhs[4]);

	/* Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix(5, 1, mxREAL);
	outy = mxGetPr(plhs[0]);
  
	/* Call the computational subroutine. */
	lcfunc(xx, oln, lam, alp, eps0, outy);
	
}

