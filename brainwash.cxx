// bug fix (8/1/11): function check_atlas() no longer checks to see if atlasID.nii and atlasID_manual.nii 
// have the same voxel dimensions.

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  brainwash.cxx                                             *
*  Copyright 2011 by Babak A. Ardekani                       *
*  ALL RIGHTS RESERVED.  No part of this program may be      *
*  used, transferred, or modified by any means without       *
*  prior written permission.  Making copies of any part      *
*  of this program for any purpose is a violation of         *
*  copyright laws.                                           *
\ * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <time.h>       
#include <sys/types.h>  
#include <sys/stat.h>  
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include "volume.h"
#include "spm_analyze.h"
#include "babak_lib.h"
#include <pthread.h> // required for multithreading
#include "smooth.h"

#define YES 1
#define NO 0
#define MAXITER 10 
#define NTHREADSMAX 256
   
///////////////////////////////////////////////////
//GLOBAL VARIABLES
char subacpcmodelfile[512]="";
char trgacpcmodelfile[512]="";

float FWHM=1.0;

int Lx,Ly,Lz;
int Wx,Wy,Wz;
int search_win=5;  // default value

int HRnx, HRny, HRnz, HRnp;
float HRdx, HRdy, HRdz;

short *HRtrghead, *HRtrgbrain;
short *HRsubhead;

float *Xw, *Yw, *Zw;

int nthreads=1; // number of threads

// searchradius[0] is for RP
// searchradius[1] is for AC
// searchradius[2] is for PC
double searchradius[3]; // in units of mm
char subOrient[4]="";  // orientation code for the subject image
char trgOrient[4]="";  // orientation code for the target (atlas) image
char Tfile[512]=""; // subject to target affine transformation file
DIM trgdim;
DIM subdim;

float AC_sub[4]={0.0, 0.0, 0.0, 1.0}; 
float PC_sub[4]={0.0, 0.0, 0.0, 1.0}; 
float RP_sub[4]={0.0, 0.0, 0.0, 1.0}; 
///////////////////////////////////////////////////

void print_matrix(const char * title, float *T);

extern float *resizeXYZ(float *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

extern short *resizeXYZ(short *image1, 
int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2);

static float v1,v2,v3,v4;
static float w1,w2;

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-h",0,'h'},
   {"-help",0,'h'},

   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},

   {"-version",0,'V'},
   {"-Version",0,'V'},
   {"-V",0,'V'},

   {"-a", 1, 'r'},     // atlas id
   {"-atlas", 1, 'r'}, // atlas id

   {"-sub", 1, 's'},
   {"-s", 1, 's'},

   {"-suborient",1,'B'},
   {"-trgorient",1,'R'},

   {"-ppm",0,'p'},
   {"-txt",0,'t'},
   {"-I",0,'I'},
   {"-iter", 1, 'i'},

   {"-AC",1,'A'},
   {"-PC",1,'P'},
   {"-VSPS",1,'S'},

   {"-iter8", 1, '8'},
   {"-iter4", 1, '4'},
   {"-iter2", 1, '2'},
   {"-iter1", 1, '1'},

   {"-o", 1, 'o'},
   {"-T", 1, 'T'},
   {"-FWHM", 1, 'F'},
   {"-threads", 1, 'H'},

   {"-rmpj",1,'0'},
   {"-w", 1, '7'},
   {"-sw", 1, '9'},
   {"-rac",1,'a'},
   {"-sd", 1, 'd'},
   {"-ms", 1, 'm'},
   {"-ma", 1, 'M'},

   {"-rpc",1,'c'},
   {0, 0, 0}
};

int opt_I=NO;
int opt_v=NO;
int opt_w=NO;
int opt_sd=NO;
int opt_cubicspline=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_help_and_exit();

short *computeReslicedImage2(short *im1, DIM dim1, DIM dim2, float *Xwarp, float *Ywarp, float *Zwarp);

short *computeReslicedImage2(float *,int,int,int,float,float,float,int,int,int,float,float,float,float *,float *,float *,float *);

void print_help_and_exit()
{
   printf("\nUsage: brainwash [-I -version -verbose -suborient <orientation code> -trgorient <orientation code>]\n"
   "[-ppm -txt -I -iter <n>] [-T <filename>] [-o <filename>] [-AC i j k] [-PC i j k] [-VSPS i j k]\n"
   "[-m <filename>] -atlas <atlasID> -sub <subjectID.nii>\n\n");
 
   printf("-atlas/-a <atlasID>: specifies the atlas to be used.  The program expects two images\n"
   "<atlasID>.nii and <atlasID>_manual.nii to exist in the $ARTHOME/atlas directory.\n\n"

   "-sub/-s <subjectID.nii>: specifies the MRI volume to be skull-stripped.\n\n"

   "-version/-V: prints software version.\n\n"

   "-verbose/-v: enables verbose mode.\n\n"

   "-suborient <orientation code>: overrides the subject image NIFTI header orientation information.\n\n"

   "-trgorient <orientation code>: overrides the atlas image NIFTI header orientation information.\n\n"

   "-ppm: outputs two images <subjectID>_ACPC_axial.ppm  and <subjectID>_ACPC_sagittal.ppm in PPM format\n"
   "that indicate the locations of the detected mid-sagittal plane and AC/PC on the subject volume.\n\n" 

   "-txt: writes the (i,j,k) coordinates of the detected AC/PC and the equation of the detected\n"
   "mid-sagittal plane to text file <subjectID>_ACPC.txt.\n\n" 

   "-I: skips the initial rigid-body alignment between subject and atlas volume based on automatic\n"
   "AC/PC and mid-sagittal plane detection.  E.g., this is required for skull-stripping rat images.\n\n"

   "-iter <n>: number of iterations used for find an affine transformation between the subject and\n"
   "atlas volumes (default=4).\n\n"

   "-T <filename>: saves the subject to atlas affine transformation in this file.\n\n"

   "-o <filename>: saves the output (subject brain mask) in this file (default=<subjectID>_<atlasID>_mask.nii).\n\n"

   "-AC i j k: overrides the automatically detected AC location on the subject volume.\n\n"

   "-PC i j k: overrides the automatically detected PC location on the subject volume.\n\n"

   "-VSPS i j k: overrides the automatically detected VSPS landmark location on the subject volume.\n\n"

   "-ms <filename>: Modelfile to be used for AC/PC detection of the subject image (default=T1acpc.mdl).\n\n"

   "-ma <filename>: Modelfile to be used for AC/PC detection of the atlas image (default=T1acpc.mdl).\n\n"
   );
   exit(0);
}

short *computeReslicedImage2(float *im1, int nx1, int ny1, int nz1, float dx1, float dy1, float dz1,
int nx2, int ny2, int nz2, float dx2, float dy2, float dz2, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
   float  x,y,z;   
   float  xx,yy,zz;   
   int q;
   int np1;
   short *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float *beta, del;
   float *c;

   if(opt_cubicspline)
   {
      beta=computeBeta(&del);
      c = (float *)calloc(nx1*ny1*nz1, sizeof(float));

      if(c==NULL)
      {
         printf("\nMemory allocation error for variable `c' in computeReslicedImage2() ...\n");
         exit(1);
      }

      cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
   }

	np1=nx1*ny1;

   im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));
   if(im2==NULL)
   {
      printf("\nMemory allocation error (im2) ...\n");
      exit(1);
   }

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

	q=0;
   
   float z0, y0;
   for(int k=0;k<nz2;k++) 
   {
      z0 = k*dz2 - zc2;
      for(int j=0;j<ny2;j++) 
      {
         y0 = j*dy2 - yc2;
         for(int i=0;i<nx2;i++) 
         {
            zz = z0 + Zwarp[q];
            yy = y0 + Ywarp[q];
            xx = i*dx2 - xc2 + Xwarp[q];

            x = ( T[0]*xx +T[1]*yy +T[2]*zz  +T[3]   + xc1 )/dx1;
            y = ( T[4]*xx +T[5]*yy +T[6]*zz  +T[7]   + yc1 )/dy1;
            z = ( T[8]*xx +T[9]*yy +T[10]*zz +T[11]  + zc1 )/dz1;

            if(opt_cubicspline)
               im2[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
            else
               im2[q++]=(short)(linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1)+0.5);
         }
      }
   }

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

void print_matrix(const char * title, float *T)
{
	printf("\n%s:",title);
	printf("\n%f\t%f\t%f\t%f",T[0],T[1],T[2],T[3]);
	printf("\n%f\t%f\t%f\t%f",T[4],T[5],T[6],T[7]);
	printf("\n%f\t%f\t%f\t%f",T[8],T[9],T[10],T[11]);
	printf("\n%f\t%f\t%f\t%f",T[12],T[13],T[14],T[15]);
	printf("\n");
}

short *computeReslicedImage2(short *im1, DIM dim1, DIM dim2, float *Xwarp, float *Ywarp, float *Zwarp)
{
   float  x,y,z;   
   int q;
   int np1;
   short *im2;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float *beta, del;
   float *c;

   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   nx1 = dim1.nx; ny1 = dim1.ny; nz1 = dim1.nz;
   nx2 = dim2.nx; ny2 = dim2.ny; nz2 = dim2.nz;

   dx1 = dim1.dx; dy1 = dim1.dy; dz1 = dim1.dz;
   dx2 = dim2.dx; dy2 = dim2.dy; dz2 = dim2.dz;

   if(opt_cubicspline)
   {
      beta=computeBeta(&del);
      c = (float *)calloc(nx1*ny1*nz1, sizeof(float));

      if(c==NULL)
      {
         printf("\nMemory allocation error for variable `c' in computeReslicedImage2() ...\n");
         exit(1);
      }

      cubicSplineAnalysis(im1, c, nx1, ny1, nz1);
   }

	np1=nx1*ny1;

   im2=(short *)calloc(nx2*ny2*nz2,sizeof(short));
   if(im2==NULL)
   {
      printf("\nMemory allocation error (im2) ...\n");
      exit(1);
   }

	xc1=dx1*(nx1-1)/2.0;     /* +---+---+ */
	yc1=dy1*(ny1-1)/2.0;
	zc1=dz1*(nz1-1)/2.0;

   	xc2=dx2*(nx2-1)/2.0;     /* +---+---+ */
	yc2=dy2*(ny2-1)/2.0;
	zc2=dz2*(nz2-1)/2.0;

   q=0;
   float z0, y0;
   for(int k=0;k<nz2;k++) 
   {
      z0 = k*dz2 - zc2 + zc1;      
      for(int j=0;j<ny2;j++) 
      {
         y0 = j*dy2 - yc2 + yc1;
         for(int i=0;i<nx2;i++) 
         {
            z = (z0 + Zwarp[q])/dz1;
            y = (y0 + Ywarp[q])/dy1;
            x = (i*dx2 - xc2 + Xwarp[q] + xc1) /dx1;

            if(opt_cubicspline)
               im2[q++] = (short)(cubicSplineSynthesis(c, nx1, ny1, nz1, x, y, z, beta, del)+0.5);
            else
               im2[q++]=linearInterpolator(x,y,z,im1,nx1,ny1,nz1,np1);
         }
      }
   }

	if(opt_cubicspline)
	{
		free(beta);
		free(c);
	}

	return( im2 );
}

void approximate_affine(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
	int np;
	int q,N;
	float xc, yc, zc;
	float rx,ry,rz;
	float sx,sy,sz;
  	float  x,y,z;   
	float *AFF, *invAFF;	// affine transform

	double Mrx, Mry, Mrz;
	double Msx, Msy, Msz;
	double SR[9], RR[9];
	double *invRR;
	double A[9],B[3];

   AFF = (float *)calloc(16,sizeof(float));
   if(AFF==NULL)
   {
      printf("\nMemory allocation error (AFF) ...\n");
      exit(1);
   }

	np = nx*ny;

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;
	zc=dz*(nz-1)/2.0;

	Mrx=Mry=Mrz=0.0;
	Msx=Msy=Msz=0.0;

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes averages of the s and r vectors defined in
	// the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	N = 0;
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]!=0.0 && Ywarp[q]!=0.0 && Zwarp[q]!=0.0)
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = x;
			sy = y;
			sz = z;

			Mrx += rx; Mry += ry; Mrz += rz;
			Msx += sx; Msy += sy; Msz += sz;

			N++;
		}
	}

	if(N!=0)
	{
		Mrx /= N; Mry /= N; Mrz /= N;
		Msx /= N; Msy /= N; Msz /= N;
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// This portion of the code computes the two 3x3 matrix in Eq. (2) of the
	// tech. notes.
	/////////////////////////////////////////////////////////////////////////////////
	for(int i=0; i<9; i++) SR[i]=RR[i]=0.0;

	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]!=0.0 && Ywarp[q]!=0.0 && Zwarp[q]!=0.0)
		{
			rx = i*dx - xc;
			ry = j*dy - yc;
			rz = k*dz - zc;

			x = rx + Xwarp[q];
			y = ry + Ywarp[q];
			z = rz + Zwarp[q];

			sx = x;
			sy = y;
			sz = z;

			rx -= Mrx; ry -= Mry; rz -= Mrz;
			sx -= Msx; sy -= Msy; sz -= Msz;

			SR[0]+=sx*rx; SR[1]+=sx*ry; SR[2]+=sx*rz;
			SR[3]+=sy*rx; SR[4]+=sy*ry; SR[5]+=sy*rz;
			SR[6]+=sz*rx; SR[7]+=sz*ry; SR[8]+=sz*rz;

			RR[0]+=rx*rx; RR[1]+=rx*ry; RR[2]+=rx*rz;
			RR[3]+=ry*rx; RR[4]+=ry*ry; RR[5]+=ry*rz;
			RR[6]+=rz*rx; RR[7]+=rz*ry; RR[8]+=rz*rz;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	// estimate A according to Eq. (2) of the technical notes.
	/////////////////////////////////////////////////////////////////////////////////
	invRR = inv3(RR);
	multi(SR,3,3,invRR,3,3,A);
	free(invRR);
	/////////////////////////////////////////////////////////////////////////////////

	// estimate B according to Eq. (1) of the technical notes
	B[0] = Msx - A[0]*Mrx - A[1]*Mry - A[2]*Mrz;
	B[1] = Msy - A[3]*Mrx - A[4]*Mry - A[5]*Mrz;
	B[2] = Msz - A[6]*Mrx - A[7]*Mry - A[8]*Mrz;

	// Eq. (3) of tech. notes
	AFF[0]=(float)A[0]; AFF[1]=(float)A[1];  AFF[2]=(float)A[2];  AFF[3]=(float)B[0];
	AFF[4]=(float)A[3]; AFF[5]=(float)A[4];  AFF[6]=(float)A[5];  AFF[7]=(float)B[1];
	AFF[8]=(float)A[6]; AFF[9]=(float)A[7]; AFF[10]=(float)A[8]; AFF[11]=(float)B[2];
	AFF[12]=0.0; AFF[13]=0.0; AFF[14]=0.0; AFF[15]=1.0;

	// Eq. (3.5) of tech. notes
	invAFF = inv4(AFF);
	delete AFF;

	///////////////////////////////////////////////////////////////////////////////

	// replace T with invAFF
	for(int i=0; i<16; i++) T[i] = invAFF[i];

	delete invAFF;
}

void fillzerovectors(int nx, int ny, int nz, float dx, float dy, float dz, float *Xwarp, float *Ywarp, float *Zwarp, float *T)
{
  	float  x,y,z;   
  	float  xx,yy,zz;   
	float xc,yc,zc;
	float *invT;		
	int q, np;

    np = nx*ny;

	invT=inv4(T);

   	xc=dx*(nx-1)/2.0;     /* +---+---+ */
	yc=dy*(ny-1)/2.0;
	zc=dz*(nz-1)/2.0;

	q=0;
	for(int k=0;k<nz;k++) 
	for(int j=0;j<ny;j++) 
  	for(int i=0;i<nx;i++) 
	{
		q = k*np + j*nx + i;

		if(Xwarp[q]==0.0 && Ywarp[q]==0.0 && Zwarp[q]==0.0)
		{
		   xx = (i*dx - xc);
		   yy = (j*dy - yc);
		   zz = (k*dz - zc);

		   x = ( invT[0]*xx +invT[1]*yy +invT[2]*zz  +invT[3]  );
		   y = ( invT[4]*xx +invT[5]*yy +invT[6]*zz  +invT[7]  );
		   z = ( invT[8]*xx +invT[9]*yy +invT[10]*zz +invT[11] );

		   Xwarp[q] = x - xx;
		   Ywarp[q] = y - yy;
		   Zwarp[q] = z - zz;
	   }
	}

	free(invT);
}

void ivf(float *Xwarp1, float *Ywarp1, float *Zwarp1, DIM dim1, float *Xwarp2, float *Ywarp2, float *Zwarp2, DIM dim2)
{
   int i2, j2, k2;
   int ii, jj, kk;
   int np1, nv1;
   int np2, nv2;
   int wx, wy, wz;
   int v1, v1_part1, v1_part2;
   int v2, v2_part1, v2_part2;

   float sigma, K;
   float *W, w;
   float xc1, yc1, zc1;
   float xc2, yc2, zc2;
   float x,y,z;
   float xx, yy, zz;
   float x2, y2, z2;
   float x1, y1, z1;

   int nx1, ny1, nz1;
   int nx2, ny2, nz2;
   float dx1, dy1, dz1;
   float dx2, dy2, dz2;

   ////////////////////////////////////////////////////////////////////////////////

   nx1 = dim1.nx; ny1 = dim1.ny; nz1 = dim1.nz;
   nx2 = dim2.nx; ny2 = dim2.ny; nz2 = dim2.nz;

   dx1 = dim1.dx; dy1 = dim1.dy; dz1 = dim1.dz;
   dx2 = dim2.dx; dy2 = dim2.dy; dz2 = dim2.dz;

   ////////////////////////////////////////////////////////////////////////////////

   nv1 = nx1 * ny1 * nz1;
   np1 = nx1 * ny1;

   nv2 = nx2 * ny2 * nz2;
   np2 = nx2 * ny2;

   ////////////////////////////////////////////////////////////////////////////////

   sigma = FWHM/2.35482; // converts FWHM to standard deviation
   K = 2.0*sigma*sigma;

   // Twice the standard deviation is greater than 95% of the area
   wx = (int)ceilf(2.0*sigma/dx2);
   wy = (int)ceilf(2.0*sigma/dy2);
   wz = (int)ceilf(2.0*sigma/dz2);

   ////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("Standard deviation of the Gaussian kernel = %f\n\n",sigma);
      printf("Support region size: %d x %d x %d voxels\n\n", 2*wx+1, 2*wy+1, 2*wz+1);

      printf("Input matrix size = %d x %d x %d\n",nx1, ny1, nz1);
      printf("Input voxel size = %7.5f x %7.5f x %7.5f mm^3\n\n",dx1, dy1, dz1);

      printf("Output matrix size = %d x %d x %d\n",nx2, ny2, nz2);
      printf("Output voxel size = %7.5f x %7.5f x %7.5f mm^3\n",dx2, dy2, dz2);
   }

   ////////////////////////////////////////////////////////////////////////////////

   W = (float *)calloc(nv2,sizeof(float));
   if(W==NULL)
   {
      printf("\nMemory allocation error (W) ...\n");
      exit(1);
   }

   xc1 = dx1 * (nx1-1.0)/2.0;
   yc1 = dy1 * (ny1-1.0)/2.0;
   zc1 = dz1 * (nz1-1.0)/2.0;

   xc2 = dx2 * (nx2-1.0)/2.0;
   yc2 = dy2 * (ny2-1.0)/2.0;
   zc2 = dz2 * (nz2-1.0)/2.0;

   for(int k1=0; k1<nz1; k1++)
   {
      v1_part1 = k1*np1;
      z1 = dz1 * k1 - zc1;
      
      for(int j1=0; j1<ny1; j1++)
      {
         v1_part2 = j1*nx1;
         y1 = dy1 * j1 - yc1;

         for(int i1=0; i1<nx1; i1++)
         {
            v1 = v1_part1 + v1_part2 + i1;
            x1 = dx1 * i1 - xc1;

            if( Xwarp1[v1]!=0.0 && Ywarp1[v1]!=0.0 && Zwarp1[v1]!=0.0 )
            {
               x2 = x1 + Xwarp1[v1];
               y2 = y1 + Ywarp1[v1];
               z2 = z1 + Zwarp1[v1];

               i2 = (int)nearbyintf( (x2 + xc2)/dx2 );
               j2 = (int)nearbyintf( (y2 + yc2)/dy2 );
               k2 = (int)nearbyintf( (z2 + zc2)/dz2 );

               for(int k=-wz; k<=wz; k++)
               {
                  kk = k+k2;
                  v2_part1 = kk*np2;

                  for(int j=-wy; j<=wy; j++)
                  {
                     jj = j+j2;
                     v2_part2 = jj*nx2;

                     for(int i=-wx; i<=wx; i++)
                     {
                        ii = i+i2;

                        if( ii>0 && jj>0 && kk>0 && ii<nx2 && jj<ny2 && kk<nz2 )
                        {
                           x = dx2 * ii - xc2;
                           y = dy2 * jj - yc2;
                           z = dz2 * kk - zc2;

                           v2 = v2_part1 + v2_part2 + ii;

                           xx = x-x2;
                           yy = y-y2;
                           zz = z-z2;

                           w = expf( -(xx*xx + yy*yy + zz*zz)/K);

                           W[v2] += w;

                           Xwarp2[v2] += w*( x1 - x);
                           Ywarp2[v2] += w*( y1 - y);
                           Zwarp2[v2] += w*( z1 - z);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   for(int i=0; i<nv2; i++)
   if(W[i]>0.0)
   {
      Xwarp2[i] /= W[i];
      Ywarp2[i] /= W[i];
      Zwarp2[i] /= W[i];
   }

   {
      float T[16];
      approximate_affine(nx2, ny2, nz2, dx2, dy2, dz2, Xwarp2, Ywarp2, Zwarp2, T);
      fillzerovectors(nx2, ny2, nz2, dx2, dy2, dz2, Xwarp2, Ywarp2, Zwarp2, T);
   }

   free(W);
}

////////////////////////////////////////////////////////////////////////////
void *warpEngine( void *ptr )
{
   int kPz, jPy, iPx;
   float Sx, Sx2, Sy, Sxy;
   float num,den;
   float CC, CCMAX;
   float V;
   short *ARsub;    // array extracted from subject image
   short *ARtrg;   // array extracted from target image
   int xopt, yopt, zopt;
   int N3;
   int kHRnp, jHRnx, voxel;
   int ki, kf; // initial and final k values
   int worker;
   int ns; // number of slices per thread

   N3 = (2*Lx+1)*(2*Ly+1)*(2*Lz+1);

   ARsub=(short *)calloc(N3,sizeof(short));
   ARtrg=(short *)calloc(N3,sizeof(short));
   if(ARsub==NULL || ARtrg==NULL)
   {
      printf("\nMemory allocation error (ARsub, ARtrg) ...\n");
      exit(1);
   }

   worker = * ((int *) ptr);
   ns = HRnz/nthreads;
   ki = worker*ns;
   if(worker == (nthreads-1))
      kf = HRnz;
   else
      kf = ki + ns;

/*
   if(opt_v)
   {
      printf("Thread %d doing slices %d through %d\n",worker, ki+1, kf);
   }
*/

   for(int k=ki; k<kf; k++)
   {
      kHRnp = k*HRnp;

      for(int j=0; j<HRny; j++)
      {
         jHRnx = j*HRnx;

         for(int i=0; i<HRnx; i++)
         {
            voxel = kHRnp + jHRnx + i;

            if(HRtrgbrain[voxel]==0)
            {
               Xw[voxel] = 0.0;
               Yw[voxel] = 0.0;
               Zw[voxel] = 0.0;
               continue;
            }

            extractArray(HRtrghead, HRnx, HRny, HRnz, HRnp, i, j, k, Lx, Ly, Lz, ARtrg);

            Sy=0.0;
            for(int n=0; n<N3; n++) Sy += ARtrg[n];

            if( Sy == 0.0 )
            {
               Xw[voxel] = 0.0;
               Yw[voxel] = 0.0;
               Zw[voxel] = 0.0;
               continue;
            }

            // IMPORTANT: we are not interested in -tive correlations
            // if CMAX is set to -1, the program fails 
            CCMAX=0.0; 	
   
            xopt=yopt=zopt=0;
      
            for(int z=-Wz; z<=Wz; z++)
            {
               kPz = k+z;
               if( kPz<0 || kPz>=HRnz) continue;

               for(int y=-Wy; y<=Wy; y++)
               {
                  jPy = j+y;
                  if( jPy<0 || jPy>=HRny) continue;

                  for(int x=-Wx; x<=Wx; x++)
                  {
                     iPx = i+x;
                     if( iPx<0 || iPx>=HRnx) continue;

                     extractArray(HRsubhead, HRnx, HRny, HRnz, HRnp, iPx, jPy, kPz, Lx, Ly, Lz, ARsub);
      
                     Sx=Sx2=Sxy=0.0;
                     for(int n=0; n<N3; n++)
                     {
                        V = ARsub[n];
                        Sx += V;
                        Sx2 += (V*V);
                        Sxy += (V*ARtrg[n]);
                     }
         
                     num = Sxy-Sx*Sy/N3;
                     if(num<0.0) continue; // important to avoid valuing -tive correlations
         
                     den = sqrtf(Sx2-Sx*Sx/N3);
                     if(den==0.0) continue;
            
                     CC = num/den;
            
                     if( CC>CCMAX ) { CCMAX=CC; xopt=x; yopt=y; zopt=z; }
                  }
               }
            }
   
            Xw[voxel] = xopt*HRdx;
            Yw[voxel] = yopt*HRdy;
            Zw[voxel] = zopt*HRdz;
         }
      }
   }

   free(ARsub);
   free(ARtrg);
}

////////////////////////////////////////////////////////////////////////////
void computeWarpField(int resFactor, short *subhead, DIM subdim, short *trghead, short *trgbrain, DIM trgdim, 
float *Xwarp, float *Ywarp, float *Zwarp, float sd, float *T)
{
   pthread_t thread[NTHREADSMAX];
   int  iret[NTHREADSMAX];
   int worker[NTHREADSMAX];

   float *invT;

   int Tnx, Tny, Tnz, Tnv;
   float Tdx, Tdy, Tdz;

   int Snx, Sny, Snz;
   float Sdx, Sdy, Sdz;

   float *Xwarp_tmp, *Ywarp_tmp, *Zwarp_tmp;

   Tnx = trgdim.nx;
   Tny = trgdim.ny;
   Tnz = trgdim.nz;
   Tdx = trgdim.dx;
   Tdy = trgdim.dy;
   Tdz = trgdim.dz;

   Snx = subdim.nx;
   Sny = subdim.ny;
   Snz = subdim.nz;
   Sdx = subdim.dx;
   Sdy = subdim.dy;
   Sdz = subdim.dz;

   Tnv=Tnx*Tny*Tnz;

   if( resFactor != 1)
   {
      HRdx = (float)( Tdx * resFactor );
      HRdy = (float)( Tdy * resFactor );
      HRdz = (float)( Tdz * resFactor );

      HRnx = (int)(1.0*Tnx/resFactor + 0.5 );
      HRny = (int)(1.0*Tny/resFactor + 0.5 );
      HRnz = (int)(1.0*Tnz/resFactor + 0.5 );
   }
   else
   {
      HRdx = Tdx;
      HRdy = Tdy;
      HRdz = Tdz;

      HRnx = Tnx;
      HRny = Tny;
      HRnz = Tnz;
   }

   if(opt_v)
   {
      printf("\nResolution reduction factor = %d",resFactor);
      printf("\nMatrix size = %d x %d x %d (voxels)", HRnx, HRny, HRnz);
      printf("\nVoxel size = %8.6f x %8.6f x %8.6f (mm3)\n", HRdx,HRdy,HRdz);
   }

   HRnp = HRnx * HRny;

   if( resFactor != 1)
   {
      Xw=resizeXYZ(Xwarp, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
      Yw=resizeXYZ(Ywarp, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
      Zw=resizeXYZ(Zwarp, Tnx, Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
   }
   else
   {
      Xw = (float *)calloc( Tnv, sizeof(float) );
      Yw = (float *)calloc( Tnv, sizeof(float) );
      Zw = (float *)calloc( Tnv, sizeof(float) );
      if(Xw==NULL || Yw==NULL || Zw==NULL)
      {
         printf("\nMemory allocation error (Xw, Yw, Zw) ...\n");
         exit(1);
      }

      for(int i=0; i<Tnv; i++)
      {
         Xw[i] = Xwarp[i];
         Yw[i] = Ywarp[i];
         Zw[i] = Zwarp[i];
      }
   }

   if( resFactor != 1)
   {
      float *tmp;
      float sdx=0.0, sdy=0.0, sdz=0.0;

      if( ( HRdx*HRdx - Sdx*Sdx )>0.0 )
      {
         sdx=(float)( sqrt( (0.5/log(2.0)) * ( HRdx*HRdx - Sdx*Sdx ) )/Sdx );
      }

      if( ( HRdy*HRdy - Sdy*Sdy )>0.0 )
      {
         sdy=(float)( sqrt( (0.5/log(2.0)) * ( HRdy*HRdy - Sdy*Sdy ) )/Sdy );
      }

      if( ( HRdz*HRdz - Sdz*Sdz )>0.0 )
      {
         sdz=(float)( sqrt( (0.5/log(2.0)) * ( HRdz*HRdz - Sdz*Sdz ) )/Sdz );
      }

      tmp = smoothXYZ(subhead,Snx,Sny,Snz,sdx,sdy,sdz);
      invT=inv4(T);
      HRsubhead=computeReslicedImage2(tmp,Snx,Sny,Snz,Sdx,Sdy,Sdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Xw, Yw, Zw, invT);
      free(invT);
      free(tmp);
   }
   else
   {
      HRsubhead = subhead;
   }

   if(resFactor != 1 )
   {
      HRtrghead=resizeXYZ(trghead, Tnx ,Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
      HRtrgbrain=resizeXYZ(trgbrain, Tnx ,Tny, Tnz, Tdx, Tdy, Tdz, HRnx, HRny, HRnz, HRdx, HRdy, HRdz);
   }
   else
   {
      HRtrghead = trghead;
      HRtrgbrain = trgbrain;
   }

   for(int i=0; i<nthreads; i++)
   {
      worker[i]=i;
      iret[i] = pthread_create( &thread[i], NULL, warpEngine, (void *)(worker+i));
   }

   for(int i=0; i<nthreads; i++)
   {
      pthread_join( thread[i], NULL);
   }

   medianFilter(Xw, HRnx, HRny, HRnz, Wx, Wy, Wz);
   medianFilter(Yw, HRnx, HRny, HRnz, Wx, Wy, Wz);
   medianFilter(Zw, HRnx, HRny, HRnz, Wx, Wy, Wz);

   if(resFactor != 1 )
   {
      Xwarp_tmp=resizeXYZ(Xw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz);
      Ywarp_tmp=resizeXYZ(Yw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz);
      Zwarp_tmp=resizeXYZ(Zw, HRnx, HRny, HRnz, HRdx, HRdy, HRdz, Tnx, Tny, Tnz, Tdx, Tdy, Tdz);
      free(Xw); free(Yw); free(Zw);
   }
   else
   {
      Xwarp_tmp = Xw;
      Ywarp_tmp = Yw;
      Zwarp_tmp = Zw;
   }

   for(int n=0; n<Tnv; n++)
   {
      Xwarp[n] += Xwarp_tmp[n];
      Ywarp[n] += Ywarp_tmp[n];
      Zwarp[n] += Zwarp_tmp[n];
   }
   free(Xwarp_tmp); free(Ywarp_tmp); free(Zwarp_tmp);

   Xwarp_tmp=smoothXYZ(Xwarp, Tnx, Tny, Tnz, sd, sd*Wy/Wx, sd*Wz/Wx);
   Ywarp_tmp=smoothXYZ(Ywarp, Tnx, Tny, Tnz, sd, sd*Wy/Wx, sd*Wz/Wx);
   Zwarp_tmp=smoothXYZ(Zwarp, Tnx, Tny, Tnz, sd, sd*Wy/Wx, sd*Wz/Wx);

   for(int i=0; i<Tnv; i++)
   {
      Xwarp[i] = Xwarp_tmp[i];
      Ywarp[i] = Ywarp_tmp[i];
      Zwarp[i] = Zwarp_tmp[i];
   }
   free(Xwarp_tmp); free(Ywarp_tmp); free(Zwarp_tmp);

   if(resFactor != 1 )
   {
      free(HRsubhead);
      free(HRtrghead);
      free(HRtrgbrain);
   }
}

////////////////////////////////////////////////////////////////////////////

DIM check_atlas(char *atlasVolumePath, char *SSatlasVolumePath)
{
   nifti_1_header hdr1, hdr2;
   DIM dim;

   if( read_NIFTI_hdr(atlasVolumePath, &hdr1) == 0 ) 
   {
      exit(0);
   }

   if( hdr1.datatype != DT_SIGNED_SHORT ) 
   {
      printf("Sorry, but this program requires %s to be of datatype 4 (i.e., DT_SIGNED_SHORT).\n"
      "The current datatype is %d.\n", atlasVolumePath, hdr1.datatype);
      exit(0);
   }

   if( read_NIFTI_hdr(SSatlasVolumePath, &hdr2) == 0 ) 
   {
      exit(0);
   }

   if( hdr2.datatype != DT_SIGNED_SHORT ) 
   {
      printf("Sorry, but this program requires %s to be of datatype 4 (i.e., DT_SIGNED_SHORT).\n"
      "The current datatype is %d.\n", SSatlasVolumePath, hdr2.datatype);
      exit(0);
   }

   if(hdr1.dim[1]!=hdr2.dim[1] || hdr1.dim[2]!=hdr2.dim[2] || hdr1.dim[3]!=hdr2.dim[3])
   {
      printf("Images %s and %s have incompatible matrix dimensions.\n", atlasVolumePath, SSatlasVolumePath);
      exit(0);
   }

/*
// this was triggering errors
// removed 08/01/2011
   if( (int)(1000*hdr1.pixdim[1]) != (int)(1000*hdr2.pixdim[1]) || 
   (int)(1000*hdr1.pixdim[2]) != (int)(1000*hdr2.pixdim[2]) ||
   (int)(1000*hdr1.pixdim[3]) != (int)(1000*hdr2.pixdim[3]) )
   {
      printf("Images %s and %s have incompatible voxel dimensions.\n", atlasVolumePath, SSatlasVolumePath);
      exit(0);
   }
*/

   if(opt_v)
   {
      printf("Atlas matrix size = %d x %d x %d\n",hdr1.dim[1], hdr1.dim[2], hdr1.dim[3]);
      printf("Atlas voxel size = %f x %f x %f\n",hdr1.pixdim[1], hdr1.pixdim[2], hdr1.pixdim[3]);
   }

   dim.nx = hdr1.dim[1]; 
   dim.ny = hdr1.dim[2]; 
   dim.nz = hdr1.dim[3];

   dim.dx = hdr1.pixdim[1]; 
   dim.dy = hdr1.pixdim[2]; 
   dim.dz = hdr1.pixdim[3];

   return(dim);
}

////////////////////////////////////////////////////////////////////////////

DIM check_subject_image(char *subImagePath)
{
   nifti_1_header hdr;
   DIM dim;

   if( read_NIFTI_hdr(subImagePath, &hdr) == 0 ) 
   {
      exit(0);
   }

   if( hdr.datatype != DT_SIGNED_SHORT && hdr.datatype != DT_UINT16) 
   {
      printf("Sorry, but this program requires %s to be of datatype 4 or 512.\n"
      "The current datatype is %d.\n", subImagePath, hdr.datatype);
      exit(0);
   }

   if(opt_v)
   {
      printf("Subject image matrix size = %d x %d x %d\n",hdr.dim[1], hdr.dim[2], hdr.dim[3]);
      printf("Subject image voxel size = %f x %f x %f\n",hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);
   }

   dim.nx = hdr.dim[1]; 
   dim.ny = hdr.dim[2]; 
   dim.nz = hdr.dim[3];

   dim.dx = hdr.pixdim[1]; 
   dim.dy = hdr.pixdim[2]; 
   dim.dz = hdr.pixdim[3];

   return(dim);
}

////////////////////////////////////////////////////////////////////////////

void find_initial_rigid_body_transformation(char *subImagePath,char *trgfile,char *sstrgfile,float *sub2trg)
{
   if(opt_I)
   {
       for(int i=0; i<16; i++) sub2trg[i]=0.0;

       sub2trg[0] = sub2trg[5] = sub2trg[10] = sub2trg[15] = 1.0;

       return;
   }

   float AC_trg[4]={0.0, 0.0, 0.0, 1.0};
   float PC_trg[4]={0.0, 0.0, 0.0, 1.0};
   float RP_trg[4]={0.0, 0.0, 0.0, 1.0};
   float Tmsp_sub[16]; // transforms the input subject volume to MSP aligned PIL orientation
   float Tmsp_trg[16]; // transforms the input target volume to MSP aligned PIL orientation
   float Tacpc_sub[16];
   float Tacpc_trg[16];
   float *invT;		

   if( searchradius[0] < 0.0) searchradius[0] = 50.0;
   if( searchradius[1] < 0.0) searchradius[1] = 15.0;
   if( searchradius[2] < 0.0) searchradius[2] = 15.0;

   detect_AC_PC_MSP(subImagePath,subOrient,subacpcmodelfile,searchradius,AC_sub,PC_sub,RP_sub,Tmsp_sub,0,opt_v,0);
   orig_ijk_to_pil_xyz(Tmsp_sub, subdim, AC_sub, PC_sub);
   ACPCtransform(Tacpc_sub, Tmsp_sub, AC_sub, PC_sub, 1);

   // these are important
   if(!opt_AC) opt_AC=YES;
   if(!opt_PC) opt_PC=YES;
   if(!opt_RP) opt_RP=YES;

   // this is necessary because sstrgfile may contain incorrect orientation information
   if(trgOrient[0]=='\0')
   {
      //readOrientationFromFile(trgfile, trgOrient);
      getNiftiImageOrientation(trgfile, trgOrient);
   }

   opt_ppm = NO; 
   opt_txt = NO; 
   detect_AC_PC_MSP(sstrgfile,trgOrient,trgacpcmodelfile,searchradius,AC_trg,PC_trg,RP_trg,Tmsp_trg,0,opt_v,0);
   orig_ijk_to_pil_xyz(Tmsp_trg, trgdim, AC_trg, PC_trg);
   ACPCtransform(Tacpc_trg, Tmsp_trg, AC_trg, PC_trg, 1);

   invT = inv4(Tacpc_trg);
   multi(invT, 4,4, Tacpc_sub, 4, 4, sub2trg);
   free(invT);

   return;
}

////////////////////////////////////////////////////////////////////////////

void find_affine_transformation(float *sub2trg, int niter, char *subfile, char *atlasfile, char *ssatlasfile, int N)
{
   short *subhead;
   short *trgbrain, *trghead;
   float sd;
   float *Xwarp, *Ywarp, *Zwarp;
   int Tnv;
   float T[16];
   short *tmp;
   char dum[1024];
   FILE *fp;
   double X[3];
   float R[9], detR;
   double d; // displacement

   Lx = (N-1)/2; if(Lx==0) Lx++;
   Ly = (int)(Lx*trgdim.dx/trgdim.dy + 0.5); if(Ly==0) Ly++;
   Lz = (int)(Lx*trgdim.dx/trgdim.dz + 0.5); if(Lz==0) Lz++;

   Wx = (search_win-1)/2; if(Wx==0) Wx++;
   Wy = (int)(Wx*trgdim.dx/trgdim.dy + 0.5); if(Wy==0) Wy++;
   Wz = (int)(Wx*trgdim.dx/trgdim.dz + 0.5); if(Wz==0) Wz++;

   sd = 2.0*Wx;

   // ensure niter is positive and less than MAXITER
   if(niter<0)
   {
      niter=4;
   }
   else if (niter>MAXITER) 
   {
      niter=MAXITER;
   }

   Tnv = trgdim.nx*trgdim.ny*trgdim.nz;
   
   // allocate memory for and initialize warp field
   Xwarp=(float *)calloc(Tnv,sizeof(float));
   Ywarp=(float *)calloc(Tnv,sizeof(float));
   Zwarp=(float *)calloc(Tnv,sizeof(float));
   if(Xwarp==NULL || Ywarp==NULL || Zwarp==NULL)
   {
      printf("\nMemory allocation error (Xwarp, Ywarp, Zwarp) ...\n");
      exit(1);
   }

   if(opt_v)
   {
      printf("\n------------------------------------------------------------------------\n");
      printf("Estimating affine transformation ...\n");
      printf("Number of iterations = %d\n",niter);
      printf("Search window = %d x %d x %d\n",2*Wx+1, 2*Wy+1, 2*Wz+1);
      printf("Correlation window size = %d x %d x %d\n",2*Lx+1, 2*Ly+1, 2*Lz+1);
   }

   for(int i=0; i<16; i++) T[i]=sub2trg[i];

   subhead = readNiftiImage(subfile, &subdim, opt_v);

   if(subhead == NULL) 
   {
      printf("\nError: Reading volume %s failed, aborting ...\n", subfile);
      exit(1);
   }

   trgbrain = readNiftiImage(ssatlasfile, &trgdim, opt_v);

   if(trgbrain == NULL) 
   {
      printf("\nError: Reading volume %s failed, aborting ...\n", ssatlasfile);
      exit(1);
   }

   trghead = readNiftiImage(atlasfile, &trgdim, opt_v);

   if(trghead == NULL) 
   {
      printf("\nError: Reading volume %s failed, aborting ...\n", atlasfile);
      exit(1);
   }

   for(int i=0; i<Tnv; i++)
   {
      if( trgbrain[i]>0 ) 
      {
         trgbrain[i] = trghead[i];
      }
      else
      {
         trgbrain[i] = 0;
      }
   }

   for(int i=1; i<=niter; i++)
   {
      if(opt_v)
      {
         printf("\nIteration %d:",i);
      }

      X[0]=T[3]; X[1]=T[7]; X[2]=T[11]; // store translations from the previous iteration

      for(int n=0; n<Tnv; n++) Xwarp[n]=Ywarp[n]=Zwarp[n]=0.0;

      if(i==1)
      computeWarpField(8, subhead, subdim, trghead, trgbrain, trgdim, Xwarp, Ywarp, Zwarp, sd, T);

      if(i==1)
      computeWarpField(4, subhead, subdim, trghead, trgbrain, trgdim, Xwarp, Ywarp, Zwarp, sd, T);

      computeWarpField(2, subhead, subdim, trghead, trgbrain, trgdim, Xwarp, Ywarp, Zwarp, sd, T);

      affineLSE(trgbrain,trgdim.nx,trgdim.ny,trgdim.nz,trgdim.dx,trgdim.dy,trgdim.dz, Xwarp,Ywarp,Zwarp,T);

      R[0]=T[0]; R[1]=T[1]; R[2]=T[2];
      R[3]=T[4]; R[4]=T[5]; R[5]=T[6];
      R[6]=T[8]; R[7]=T[9]; R[8]=T[10];
		
      d = sqrt( (X[0]-T[3])*(X[0]-T[3]) + (X[1]-T[7])*(X[1]-T[7]) + (X[2]-T[11])*(X[2]-T[11]) );
      detR = det3(R);

      if(opt_v)
      {
         print_matrix("subject -----> target affine transformation",T);
         printf("\nScale change = %f",detR);
         printf("\nTranslation relative to previous iteration = %lf mm\n",d);
      }
   }

   for(int i=0; i<16; i++) sub2trg[i]=T[i];

   if(opt_v)
   {
      printf("\n------------------------------------------------------------------------\n");
   }

   if(Tfile[0] != '\0')
   {
      // save the final affine transformation in the specified outputfile
      sprintf(dum, "%s to %s affine transformation:", subfile, ssatlasfile);
      fp = fopen(Tfile, "w");
      if(fp==NULL)
      {
         printf("\nError opening %s!\n", Tfile);
      }
      printMatrix(sub2trg, 4, 4, dum, fp);
      fclose(fp);
   }

   free(subhead);
   free(trgbrain);
   free(trghead);
   free(Xwarp);
   free(Ywarp);
   free(Zwarp);
}

////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   // holds the atlasID specified at the command line
   char atlasID[512]="";  // intialization to NULL is important
   char atlasVolumePath[512]="";   // this will be set to $ARTHOME/atlas/<atlasID>.nii
   char SSatlasVolumePath[512]=""; // this will be set to $ARTHOME/atlas/<atlasID>_manual.nii
   char subImagePath[512]=""; 
   char outputfilename[512]="";
   char outputpath[512]=""; // where the final mask image will be written
   char subject_dir[512]=""; 

   float sub2trg[16];

   int niter=4; // number of iterations used in finding the initial affine transformation 
   int iter8=0;
   int iter4=1;
   int iter2=1;
   int iter1=1;

   float T[16];		// overall rigid body transformation
   float *Xwarp, *Ywarp, *Zwarp;

   short *subhead;
   int Snv;

   short *trgbrain, *trghead;
   int Tnv;

   float sd;
   int N;	// N=2*L+1

   /////////////////////////////////////////////////////////////////////
   // maximize the stack size
   /////////////////////////////////////////////////////////////////////
   {
      // define a rlimit structure
      struct rlimit rlpty;

      // get the resource limit of type RLIMIT_STACK for stack size
      (void) getrlimit(RLIMIT_STACK,&rlpty);

      // set the stack size to its maximum
      rlpty.rlim_cur=rlpty.rlim_max;
      (void) setrlimit(RLIMIT_STACK,&rlpty);
   }
   /////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////
   // important initializations

   searchradius[0] = 50.0;
   searchradius[1] = 15.0;
   searchradius[2] = 15.0;

   opt_ppm = NO;
   opt_txt = NO;

   /////////////////////////////////////////////////////////////////////

   if(argc==1)
   {
      print_help_and_exit();
   }

   /////////////////////////////////////////////////////////////////////

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'm':
            sprintf(subacpcmodelfile,"%s",optarg);
            break;
         case 'M':
            sprintf(trgacpcmodelfile,"%s",optarg);
            break;
         case 'v':
            opt_v=YES;
            break;
         case 'V':
            printf("Version: July 28, 2011\n");
            exit(0);
         case 'r':
            sprintf(atlasID,"%s",optarg);
            break;
         case 's':
            sprintf(subImagePath,"%s",optarg);
            break;
         case 'A':
            AC_sub[0] = atoi(argv[optind-1]);
            AC_sub[1] = atoi(argv[optind+0]);
            AC_sub[2] = atoi(argv[optind+1]);
            opt_AC = NO;
            break;
         case 'I':
            opt_I = YES;
            break;
         case 'P':
            PC_sub[0] = atoi(argv[optind-1]);
            PC_sub[1] = atoi(argv[optind+0]);
            PC_sub[2] = atoi(argv[optind+1]);
            opt_PC = NO;
            break;
         case 'S':
            RP_sub[0] = atoi(argv[optind-1]);
            RP_sub[1] = atoi(argv[optind+0]);
            RP_sub[2] = atoi(argv[optind+1]);
            opt_RP = NO;
            break;
         case 'p':
            opt_ppm = YES;
            break;
         case 't':
            opt_txt = YES;
            break;
         case '8':
            iter8=atoi(optarg);
            break;
         case '4':
            iter4=atoi(optarg);
            break;
         case '2':
            iter2=atoi(optarg);
            break;
         case '1':
            iter1=atoi(optarg);
            break;
         case 'h':
            print_help_and_exit();
         case 'F':
            FWHM=atof(optarg);
            break;
         case 'H':
            nthreads=atoi(optarg);
            break;
         case 'i':
            niter=atoi(optarg);
            break;
         case 'T':
            sprintf(Tfile,"%s",optarg);
            break;
         case 'B':
            sprintf(subOrient,"%s",optarg);
            break;
         case 'R':
            sprintf(trgOrient,"%s",optarg);
            break;
         case '0':
            searchradius[0] = atof(optarg);
            break;
         case 'a':
            searchradius[1] = atof(optarg);
            break;
         case 'c':
            searchradius[2] = atof(optarg);
            break;
         case '7':
            N=atoi(optarg);
            opt_w=YES;
            break;
         case '9':
            search_win=atoi(optarg);
            break;
         case 'd':
            sd = atof(optarg);
            opt_sd=YES;
            break;
         case 'o':
            sprintf(outputfilename,"%s",optarg);
            break;
         case '?':
            print_help_and_exit();
		}
	}

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // some initializations necessary for affine and non-linear transformation finding

   if(!opt_w || N<=0) N=5;
   if( (N%2)==0 ) N+=1;			// make sure it's odd
   if(N==1) N=3;

   if(search_win<=0) search_win=5;
   if( (search_win%2)==0 ) search_win+=1;      // make sure it's odd
   if(search_win==1) search_win=3;             // make sure it's at least 3

   if(nthreads<0) nthreads=1;
   if(nthreads>NTHREADSMAX) nthreads=NTHREADSMAX;

   if(nthreads>1 && opt_v)
   {
      printf("Number of threads = %d\n",nthreads);
   }
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // get the value of the ARTHOME environment variable
   // The getenv() function searches the environment list for a string that matches "ARTHOME".
   // It returns a pointer to the value in the environment, or NULL if there is no match.

   char *ARTHOME;  // full path of the directory of the ART software

   ARTHOME=getenv("ARTHOME");

   if(ARTHOME == NULL)
   {
      printf("The ARTHOME environment variable is not defined.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("ARTHOME = %s\n",ARTHOME);
   }
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////
   // Ensure that an atlas volume has been specified at the command line.
   if( atlasID[0] == '\0')
   {
      printf("Please specify an atlas ID using -atlas argument.\n");
      exit(0);
   }

   sprintf(atlasVolumePath,"%s/atlas/%s.nii",ARTHOME, atlasID);
   sprintf(SSatlasVolumePath,"%s/atlas/%s_manual.nii",ARTHOME, atlasID);

   if(opt_v)
   {
      printf("Atlas ID = %s\n",atlasID);
      printf("Atlas volume = %s\n",atlasVolumePath);
      printf("Skull-stripped atlas volume = %s\n",SSatlasVolumePath);
   }

   trgdim = check_atlas(atlasVolumePath, SSatlasVolumePath);
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////
   if( subImagePath[0] == '\0')
   {
      printf("Please specify a subject image using -sub argument.\n");		
      exit(0);
   }

   if(opt_v)
   {
      printf("Subject image (MRI volume to be skull-stripped) = %s\n",subImagePath);
   }

   subdim = check_subject_image(subImagePath);
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////////////////
   // if outputfilename is not specified, use the defult: "subID_trgID_mask.nii"
   if( outputfilename[0] == '\0')
   {
      char subjectID[512];

      niftiFilename(subjectID, subImagePath);

      sprintf(outputfilename,"%s_%s_mask.nii",subjectID,atlasID);
   }

   getDirectoryName(subImagePath, subject_dir);
   sprintf(outputpath,"%s/%s",subject_dir, outputfilename);

   if(opt_v)
   {
      printf("Output file = %s\n", outputpath);
   }
   ////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////
   find_initial_rigid_body_transformation(subImagePath, atlasVolumePath, SSatlasVolumePath, sub2trg);

   if(opt_v)
   {
      print_matrix("Initial subject --> target affine transformation",sub2trg);
   }

   find_affine_transformation(sub2trg, niter, subImagePath, atlasVolumePath, SSatlasVolumePath, N);
   ////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////
   if(opt_v)
   {
      printf("\n------------------------------------------------------------------------\n");
      printf("Estimating non-linear transformation ...\n");
      printMatrix(sub2trg, 4, 4, "Affine transformation:", NULL);
   }

   Lx = (N-1)/2; if(Lx==0) Lx++;
   Ly = (int)(Lx*trgdim.dx/trgdim.dy + 0.5); if(Ly==0) Ly++;
   Lz = (int)(Lx*trgdim.dx/trgdim.dz + 0.5); if(Lz==0) Lz++;

   Wx = (search_win-1)/2; if(Wx==0) Wx++;
   Wy = (int)(Wx*trgdim.dx/trgdim.dy + 0.5); if(Wy==0) Wy++;
   Wz = (int)(Wx*trgdim.dx/trgdim.dz + 0.5); if(Wz==0) Wz++;

   if(!opt_sd || sd<0.0)
      sd = 2.0*Wx;

   subhead = readNiftiImage(subImagePath, &subdim, opt_v);
   if(subhead == NULL) 
   {
      printf("\nError: Reading volume %s failed, aborting ...\n\n", subImagePath);
      exit(0);
   }

   trgbrain = readNiftiImage(SSatlasVolumePath, &trgdim, opt_v);
   if(trgbrain == NULL) 
   {
      printf("\nError: Reading volume %s failed, aborting ...\n\n", SSatlasVolumePath);
      exit(0);
   }

   trghead = readNiftiImage(atlasVolumePath, &trgdim, opt_v);
   if(trghead == NULL) 
   {
      printf("\nError: Reading volume %s failed, aborting ...\n\n", atlasVolumePath);
      exit(0);
   }

   Tnv = trgdim.nx*trgdim.ny*trgdim.nz;

   for(int i=0; i<Tnv; i++)
   {
      if( trgbrain[i]>0 ) 
      {
         trgbrain[i] = trghead[i];
      }
      else
      {
         trgbrain[i] = 0;
      }
   }

   // allocate memory for and initialize warp field
   Xwarp=(float *)calloc(Tnv,sizeof(float));
   Ywarp=(float *)calloc(Tnv,sizeof(float));
   Zwarp=(float *)calloc(Tnv,sizeof(float));
   if(Xwarp==NULL || Ywarp==NULL || Zwarp==NULL)
   {
      printf("\nMemory allocation error (Xwarp, Ywarp, Zwarp) ...\n");
      exit(1);
   }

   for(int n=0; n<Tnv; n++) Xwarp[n]=Ywarp[n]=Zwarp[n]=0.0;

   // Removing this (i.e. setting iter8=0) substantially improved the program
   for(int i=0; i<iter8; i++)
   {
      sd = 8.0;
      computeWarpField(8,subhead,subdim,trghead,trgbrain,trgdim,Xwarp,Ywarp,Zwarp,sd,sub2trg);
   }

   for(int i=0; i<iter4; i++)
   {
      sd = 6.0;
      computeWarpField(4,subhead,subdim,trghead,trgbrain,trgdim,Xwarp,Ywarp,Zwarp,sd,sub2trg);
   }

   for(int i=0; i<iter2; i++)
   {
      sd = 3.0;
      computeWarpField(2,subhead,subdim,trghead,trgbrain,trgdim,Xwarp,Ywarp,Zwarp,sd,sub2trg);
   }

   for(int i=0; i<iter1; i++)
   {
      sd = 2.0;
      computeWarpField(1,subhead,subdim,trghead,trgbrain,trgdim,Xwarp,Ywarp,Zwarp,sd,sub2trg);
   }

   combine_warps_and_trans(trgdim.nx,trgdim.ny,trgdim.nz,trgdim.dx,trgdim.dy,trgdim.dz,Xwarp,Ywarp,Zwarp,sub2trg);

   if(opt_v)
   {
      printf("\n------------------------------------------------------------------------\n");
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////

   float *iXwarp, *iYwarp, *iZwarp;
   short *subbrain;
   unsigned char *uc_tmp;
   nifti_1_header hdr;
   nifti1_extender ext;
   FILE *fp;

   Snv = subdim.nx*subdim.ny*subdim.nz;

   iXwarp = (float *)calloc(Snv,sizeof(float));
   iYwarp = (float *)calloc(Snv,sizeof(float));
   iZwarp = (float *)calloc(Snv,sizeof(float));
   if(iXwarp==NULL || iYwarp==NULL || iZwarp==NULL)
   {
      printf("\nMemory allocation error (iXwarp, iYwarp, iZwarp) ...\n");
      exit(1);
   }

   ivf(Xwarp, Ywarp, Zwarp, trgdim, iXwarp, iYwarp, iZwarp, subdim);

   for(int i=0; i<Tnv; i++) 
   if(trgbrain[i]>0) trgbrain[i]=100;

   subbrain=computeReslicedImage2(trgbrain,trgdim,subdim,iXwarp,iYwarp,iZwarp);

   if( read_NIFTI_hdr(subImagePath, &hdr) == 0 ) 
   {
      exit(0);
   }

   sprintf(hdr.descrip,"Created by ART `brainwash' program.");
   hdr.vox_offset = 352.0; // important
   ext.extension[0] = 0; // important
   hdr.datatype = DT_UINT8;  // sets the datatype to unsigned char
   hdr.bitpix = 8;  // sets the number of bits per voxel

   uc_tmp = (unsigned char *)calloc(Snv,1);
   if(uc_tmp==NULL)
   {
      printf("\nMemory allocation error (uc_tmp)!\n");
      exit(1);
   }

   // copy contents of subbrain to uc_tmp
   for(int i=0; i<Snv; i++) uc_tmp[i] = (unsigned char)subbrain[i];

   fp = fopen(outputpath,"w");
   if(fp==NULL)
   {
      printf("\nError opening %s!\n", outputpath);
      exit(1);
   }

   if( fwrite(&hdr, sizeof(nifti_1_header), 1, fp) != 1 )
   {
      printf("\nError writing header to %s!\n", outputpath);
      exit(1);
   }

   if( fwrite(&ext, sizeof(nifti1_extender), 1, fp) != 1 ) 
   {
      printf("\nError writing header extension to %s!\n", outputpath);
      exit(1);
   }

   if( fwrite(uc_tmp, 1, Snv, fp) != Snv )
   {
      printf("\nError writing voxel data to %s!\n", outputpath);
      exit(1);
   }
   fclose(fp);

   free(uc_tmp);
   free(subbrain);
   free(iXwarp); free(iYwarp); free(iZwarp);
   free(Xwarp); free(Ywarp); free(Zwarp);
   free(subhead);
   free(trgbrain);
   free(trghead);
}
