#include "locf.h"

/* Here is the main function for calculating the Jacobian numerically. */
/* Note that this routine requires xx, enns, and avenns (along with */
/* some other stuff). */

void njacob(double lam, double *pts, double *fd, double *alp, double eps0, 
				int tp, int cdm, double *xx, double *enns, double *avenns,
				double *JJ, int *ir, int *jc)
{

	int ii, jj, kk, nn, pp, bb, qq, pose, posx;
	int iisc, ipsc, imsc, jjsc, jpsc, jmsc, kp, km;
	int col, posr;

	int pt0, pt1, pt2;						/* Stupid type forcing. */
	int ibk = (int)(pts[1]*pts[2]);
	int aps = (int)(pts[0]*pts[1]*pts[2]);

	double jig=1e-8;
	int crnr[27];						/* Indices of box corners. */
	int llft[] = {0,9,12,3,1,10,13,4};	/* Lower-left corners of meta-box. */
	double xit[64];						/* Array to pass to lcfunc. */
	double xxt[8]; 						/* Temporary array. */
	double oln[40];						/* Average normals. */
	double fft[5];						/* The function value. */

	/* pts has to be a pointer to a double. */
	pt0=(int)pts[0];
	pt1=(int)pts[1];
	pt2=(int)pts[2];

	jc[0]=0; 			/* This should always be zero, right? */
	posr=0;				/* Row position. */
	col=-1;				/* Column position. */

	/* Give props to the user. 
	printf(" Generating Jacobian Entries: ii=    ");*/

	/* Loop it column-wise. */
	/* Don't indent to save space. */
	for (ii=0; ii<pt0; ii++)	{
		
		imsc=((ii-1+pt0)%pt0)*ibk;
		iisc=ii*ibk;
		ipsc=((ii+1)%pt0)*ibk;

		/* Moving right along.
		printf("\b\b\b%3d",ii+1);*/

	for (jj=0; jj<pt1; jj++)	{

		jmsc=((jj-1+pt1)%pt1)*pt2;
		jjsc=jj*pt2;
		jpsc=((jj+1)%pt1)*pt2;

	for (kk=0; kk<pt2; kk++)	{

		km=(kk-1+pt2)%pt2;
		kp=(kk+1)%pt2;

		/* Now we have the point. */
		/* Build the boxes that the point will affect. */
		crnr[0]=imsc+jmsc+km;	crnr[3]=imsc+jjsc+km;	crnr[6]=imsc+jpsc+km;
		crnr[1]=imsc+jmsc+kk;	crnr[4]=imsc+jjsc+kk;	crnr[7]=imsc+jpsc+kk;
		crnr[2]=imsc+jmsc+kp;	crnr[5]=imsc+jjsc+kp;	crnr[8]=imsc+jpsc+kp;

		crnr[9]=iisc+jmsc+km;	crnr[12]=iisc+jjsc+km;	crnr[15]=iisc+jpsc+km;
		crnr[10]=iisc+jmsc+kk;	crnr[13]=iisc+jjsc+kk;	crnr[16]=iisc+jpsc+kk;
		crnr[11]=iisc+jmsc+kp;	crnr[14]=iisc+jjsc+kp;	crnr[17]=iisc+jpsc+kp;
		
		crnr[18]=ipsc+jmsc+km;	crnr[21]=ipsc+jjsc+km;	crnr[24]=ipsc+jpsc+km;
		crnr[19]=ipsc+jmsc+kk;	crnr[22]=ipsc+jjsc+kk;	crnr[25]=ipsc+jpsc+kk;
		crnr[20]=ipsc+jmsc+kp;	crnr[23]=ipsc+jjsc+kp;	crnr[26]=ipsc+jpsc+kp;

		/* We need to re-sort llft according to row-ordering. */
		/* This does get a little goofy. I'm using my own bubble sort. */
		for (pp=7; pp>0; pp--)
			for (qq=0; qq<pp; qq++)
				if (crnr[llft[qq]] > crnr[llft[qq+1]])	{
					nn=llft[qq];
					llft[qq] = llft[qq+1];
					llft[qq+1]=nn;
				}

		/* Tweak each direction individually. */
		pose=crnr[13];
		for (nn=0; nn<5; nn++){

			/* Each column has 40 elements in it. */
			jc[++col+1]=jc[col]+40;

			/* Jiggle the point in the correct direction. */
			posx=crnr[13];
			for (pp=0; pp<8; pp++) {
				xxt[pp]=xx[posx];
				xx[posx]+=jig*enns[pose];
				pose+=aps;
				posx+=aps;
			}

			/* Generate the box O' points (8 boxes). */
			for (bb=0; bb<8; bb++)	{

				/* This loop is over the eight components. */
				posx=0; qq=0;
				for (pp=0; pp<64; pp+=8)	{

					/* Fill in the points. */
					xit[0+pp]=xx[ posx + crnr[llft[bb]] ]; 
					xit[1+pp]=xx[ posx + crnr[llft[bb]+9] ]; 
					xit[2+pp]=xx[ posx + crnr[llft[bb]+12] ]; 
					xit[3+pp]=xx[ posx + crnr[llft[bb]+3] ]; 
					xit[4+pp]=xx[ posx + crnr[llft[bb]+1] ]; 
					xit[5+pp]=xx[ posx + crnr[llft[bb]+10] ]; 
					xit[6+pp]=xx[ posx + crnr[llft[bb]+13] ]; 
					xit[7+pp]=xx[ posx + crnr[llft[bb]+4] ]; 

					/* Fill in the average normal vectors. */
					oln[0+qq]=avenns[crnr[llft[bb]]+posx];
					oln[1+qq]=avenns[crnr[llft[bb]]+posx+8*aps];
					oln[2+qq]=avenns[crnr[llft[bb]]+posx+16*aps];
					oln[3+qq]=avenns[crnr[llft[bb]]+posx+24*aps];
					oln[4+qq]=avenns[crnr[llft[bb]]+posx+32*aps];

					/* Update the positions. */
					posx+=aps;
					qq+=5;
				}

				/* Now that we have the box, we can calculate the function. */
				lcfunc(xit, oln, lam, alp, eps0, fft);

				/* Divide the difference by jiggle. */
				/* Stuff it into the Jacobian. */
				for (pp=0; pp<5; pp++)	{

					ir[posr] = 5*crnr[llft[bb]]+pp;
					JJ[posr] = ( fft[pp]-fd[ir[posr]] )/jig;
					posr++;
				}

			}

			/* Un-jiggle the point. */
			posx=crnr[13];
			for (pp=0; pp<8; pp++) {
				xx[posx]=xxt[pp];
				posx+=aps;
			}

			
		}
		

	}
	}
	}

}


/* "Gateway" routine. */
/* As above, note that we have to pass in the variables 
	xx, enns, and avenns. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  	double *fd, *alp, *xx, *enns, *avenns, *JJ;
	int *ir, *jc;
	double lam, eps0;
	double *pts;			/* I don't like this. */
	int mm, tp, cdm;
  
	/* I'm not doing ANY error checking on the inputs. Be careful. */

	/* Assign pointers to input and get the scalar. */
	lam = mxGetScalar(prhs[0]);
	pts = mxGetPr(prhs[1]);
	fd = mxGetPr(prhs[2]);
	alp = mxGetPr(prhs[3]);
	eps0 = mxGetScalar(prhs[4]);
	xx = mxGetPr(prhs[5]);
	enns = mxGetPr(prhs[6]);
	avenns = mxGetPr(prhs[7]);

	/* Get the dimensions of the matrix, and  */
	/* get the pointers to the parts. */
	mm = mxGetM(prhs[5]);
	tp = mxGetM(prhs[2]);
	cdm = (int)tp/mm;

	/* Set the output pointer to the output matrix. */
	plhs[0] = mxCreateSparse(tp, tp, 8*mm*cdm*cdm, mxREAL);
	JJ = mxGetPr(plhs[0]);
	ir = mxGetIr(plhs[0]);
	jc = mxGetJc(plhs[0]);

	/* Call the computational subroutine. */
	njacob(lam, pts, fd, alp, eps0, tp, cdm, xx, enns, avenns, 
				JJ, ir, jc);
	
}

