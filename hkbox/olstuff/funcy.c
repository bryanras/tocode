#include "mex.h"

/* This is the vector field calculation, including the "gateway" routine. */

/* This function returns the vector field. */
/* It is basically the same thing as funcy.m */
void fcy(double *xx, double b, double c, double d, double lam, double *outy)
{
    double dit;

    /* Don't do this calculation more than once. */
    dit=xx[2]+d*(1-xx[2]*xx[2]);

    /* This is everything. */
    outy[0] = (lam-b)*xx[0]-c*xx[1]+xx[0]*dit;
    outy[1] = c*xx[0]+(lam-b)*xx[1]+xx[1]*dit;
    outy[2] = lam*xx[2]-xx[0]*xx[0]-xx[1]*xx[1]-xx[2]*xx[2];

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  	double *xx, *outy;
	double b, c, d, lam;
  
	/* I'm not doing ANY error checking on the inputs. Be careful. */

	/* Assign pointers to input and get the scalar. */
	xx = mxGetPr(prhs[0]);
	b = mxGetScalar(prhs[1]);
	c = mxGetScalar(prhs[2]);
	d = mxGetScalar(prhs[3]);
	lam = mxGetScalar(prhs[4]);

	/* Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
	outy = mxGetPr(plhs[0]);
  
	/* Call the computational subroutine. */
	fcy(xx, b, c, d, lam, outy);
	
}

