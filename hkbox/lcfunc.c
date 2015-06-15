#include "locf.h"
#include <math.h>

/* This is the local function calculation, including the "gateway" routine. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  	double *outy;
	double *xll, *xul, *xur, *xlr;
	double b, c, d, lam;
	double *ff;
	unsigned int ii;
	double xx[3];
  
	/* I'm not doing ANY error checking on the inputs. Be careful. */

	/* Assign pointers to input and get the scalar. */
	xll = mxGetPr(prhs[0]);
	xlr = mxGetPr(prhs[1]);
	xul = mxGetPr(prhs[2]);
	xur = mxGetPr(prhs[3]);
	b = mxGetScalar(prhs[4]);
	c = mxGetScalar(prhs[5]);
	d = mxGetScalar(prhs[6]);
	lam = mxGetScalar(prhs[7]);

	/* Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	outy = mxGetPr(plhs[0]);
  
	/* Subtract and average the vectors. */
	for (ii=0; ii<3; ii++) {
		xx[ii] = 0.25*(xll[ii] + xul[ii] + xlr[ii] + xur[ii]);
		xur[ii]-=xll[ii];
		xul[ii]-=xlr[ii];
	}

	/* Call the computational subroutine. */
	ff=fcy(xx, b, c, d, lam);

	/* Store the crossed product in xll. */
	xll[0] = xur[1]*xul[2] - xur[2]*xul[1];
	xll[1] = xur[2]*xul[0] - xur[0]*xul[2];
	xll[2] = xur[0]*xul[1] - xur[1]*xul[0];
 
	/* Store the norm of n in b. */
	b = sqrt(xll[0]*xll[0] + xll[1]*xll[1] + xll[2]*xll[2]);

	/* Calculate f*v. */
	*outy=0.0;
	for (ii=0; ii<3; ii++)
		*outy+=xll[ii]*ff[ii];

	*outy/=b;					/* Normalization of n. */
}

