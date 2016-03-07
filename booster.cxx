/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  brainwash2_merge.cxx                                      *
*  Copyright 2010 by Babak A. Ardekani                       *
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
#include "permutation.h"

#define YES 1
#define NO 0
   
//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-n", 1, 'n'},

   {"-t", 1, 't'},
   {"-threshold", 1, 't'},

   {"-output", 1, 'o'},
   {"-o", 1, 'o'},

   {"-subject", 1, 's'},
   {"-sub", 1, 's'},
   {"-s", 1, 's'},

   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},

   {0, 0, 0}
};

char opt_v=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
float jaccardindex(unsigned char *setA, unsigned char *setB, int n)
{
   int sizeAorB=0;
   int sizeAandB=0;

   if(n<=0) 
   {
      printf("\n\nWarning: A non-positive array dimension passed to the diceindex() function.\n\n");
      return(0.0);
   }

   if( setA==NULL || setB==NULL) 
   {
      printf("\n\nWarning: A NULL array passed to the diceindex() function.\n\n");
      return(0.0);
   }

   for(int i=0; i<n; i++)
   {
      if( setA[i] > 0 || setB[i]>0)
         sizeAorB++;

      if( setA[i] > 0 && setB[i] > 0)
         sizeAandB++;
   }

   if( (sizeAorB) == 0 )
      return(0.0);

   return( sizeAandB/(1.0*sizeAorB) );
}
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
float diceindex(unsigned char *setA, unsigned char *setB, int n)
{
   int sizeA=0, sizeB=0, sizeAandB=0;

   if(n<=0) 
   {
      printf("\n\nWarning: A non-positive array dimension passed to the diceindex() function.\n\n");
      return(0.0);
   }

   if( setA==NULL || setB==NULL) 
   {
      printf("\n\nWarning: A NULL array passed to the diceindex() function.\n\n");
      return(0.0);
   }

   for(int i=0; i<n; i++)
   {
      if( setA[i] > 0 )
         sizeA++;

      if( setB[i] > 0 )
         sizeB++;

      if( setA[i] > 0 && setB[i] > 0)
         sizeAandB++;
   }

   if( (sizeA + sizeB) == 0 )
      return(0.0);

   return( 2.0*sizeAandB/(sizeA+sizeB) );
}
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
void tt(PERMUTATION *maskIndx, int n, char **maskList, unsigned char *consensusmask, float *d)
{
   nifti_1_header hdr;
   unsigned char *mask;
   int nv;

   for(int i=0; i<n; i++)
   {
      mask = (unsigned char *)read_nifti_image(maskList[ maskIndx->vec[i] ], &hdr);

      nv = hdr.dim[1]*hdr.dim[2]*hdr.dim[3];

      for(int v=0; v<nv; v++)
      if(mask[v] >= 50) mask[v]=1; else mask[v]=0;

      //d[i]=diceindex(mask, consensusmask, nv);
      d[i]=jaccardindex(mask, consensusmask, nv);

      d[i] = d[i]*d[i];

      printf("\tn=%d w=%f\n",i+1,d[i]);
      delete mask;
   }
}
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

float *compute_avg_mask(PERMUTATION *maskIndx, int n, char **maskList, float *d)
{
   nifti_1_header hdr;
   int nx, ny, nz;
   float dx, dy, dz;
   int nv;
   int type;
   char *mask;
   float *avg_mask;
   float sum_d=0.0;

/*
   if(opt_v)
   {
      printf("Masks used for skull-stripping:\n");
      for(int i=0; i<n; i++)
      {
         printf("%d\t%s\n", i+1, maskList[ maskIndx->vec[i] ] );
      }
   }
*/

   if(n==0) 
   {
      return(NULL);
   }

/*
   if(opt_v)
   {
      printf("Reading mask %d: %s ...\n",1, maskList[ maskIndx->vec[0] ]);
   }
*/

   mask = read_nifti_image(maskList[ maskIndx->vec[0] ], &hdr);
   if(mask==NULL)
   {
      return(NULL);
   }
   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];
   type=hdr.datatype;

/*
   if(opt_v)
   {
      printf("\tMatrix size = %d x %d x %d\n", nx, ny, nz);
      printf("\tVoxe size = %f x %f x %f\n", dx, dy, dz);
   }
*/

   nv = nx*ny*nz;
   avg_mask = (float *)calloc(nv, sizeof(float));
   if(avg_mask == NULL)
   {
      return(NULL);
   }

   sum_d=0.0;
   for(int i=0; i<n; i++)
      sum_d += d[i];

printf("\n%f\n",sum_d);

   switch(type) {
      case 2:
         for(int i=0; i<nv; i++) avg_mask[i] = ((unsigned char *)mask)[i]*d[0]/sum_d;
         break;
      case 4:
         for(int i=0; i<nv; i++) avg_mask[i] = ((short *)mask)[i]*d[0]/sum_d;
         break;
      case 8:
         for(int i=0; i<nv; i++) avg_mask[i] = ((int *)mask)[i]*d[0]/sum_d;
         break;
      case 16:
         for(int i=0; i<nv; i++) avg_mask[i] = ((float *)mask)[i]*d[0]/sum_d;
         break;
   }

   delete mask;

   for(int i=1; i<n; i++)
   {
/*
      if(opt_v)
      {
         printf("Reading mask %d: %s ...\n",i+1, maskList[ maskIndx->vec[i] ]);
      }
*/

      mask = read_nifti_image(maskList[ maskIndx->vec[i] ], &hdr);
      if(mask==NULL) 
      {
         return(NULL);
      }

      dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

/*
      if(opt_v)
      {
         printf("\tMatrix size = %d x %d x %d\n", hdr.dim[1], hdr.dim[2], hdr.dim[3]);
         printf("\tVoxe size = %f x %f x %f\n", hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);
      }
*/

      if( nv != hdr.dim[1]*hdr.dim[2]*hdr.dim[3] )
      {
         printf("Error: Input masks must all have the same matrix size.\n"); 
         return(NULL);
      }

      switch(type) {
         case 2:
            for(int j=0; j<nv; j++) avg_mask[j] += ((unsigned char *)mask)[j]*d[i]/sum_d;
            break;
         case 4:
            for(int j=0; j<nv; j++) avg_mask[j] += ((short *)mask)[j]*d[i]/sum_d;
            break;
         case 8:
            for(int j=0; j<nv; j++) avg_mask[j] += ((int *)mask)[j]*d[i]/sum_d;
            break;
         case 16:
            for(int j=0; j<nv; j++) avg_mask[j] += ((float *)mask)[j]*d[i]/sum_d;
            break;
      }

      delete mask;
   }

   //for(int i=0; i<nv; i++) avg_mask[i] /= n;

   return(avg_mask);
}

#if 0
float *compute_avg_mask(PERMUTATION *maskIndx, int n, char **maskList)
{
   int *nvec;
   int *increment;
   nifti_1_header hdr;
   int nx, ny, nz;
   float dx, dy, dz;
   int np, nv;
   int type;
   char *mask;
   float *avg_mask;
   float *mask_as_float;

   if(opt_v)
   {
      printf("Masks used for skull-stripping:\n");
      for(int i=0; i<n; i++)
      {
         printf("%d\t%s\n", i+1, maskList[ maskIndx->vec[i] ] );
      }
   }

   if(n==0) 
   {
      return(NULL);
   }

/*
   if(opt_v)
   {
      printf("Reading mask %d: %s ...\n",1, maskList[ maskIndx->vec[0] ]);
   }
*/

   mask = read_nifti_image(maskList[ maskIndx->vec[0] ], &hdr);
   if(mask==NULL)
   {
      return(NULL);
   }
   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];
   type=hdr.datatype;

/*
   if(opt_v)
   {
      printf("\tMatrix size = %d x %d x %d\n", nx, ny, nz);
      printf("\tVoxe size = %f x %f x %f\n", dx, dy, dz);
   }
*/

   nv = nx*ny*nz;
   np = nx*ny;
   avg_mask = (float *)calloc(nv, sizeof(float));
   mask_as_float = (float *)calloc(nv, sizeof(float));
   if(avg_mask == NULL)
   {
      return(NULL);
   }

   nvec = (int *)calloc(nz, sizeof(int));
   increment = (int *)calloc(nz, sizeof(int));
   for(int k=0; k<nz; k++) nvec[k]=0;

   switch(type) {
      case 2:
         for(int i=0; i<nv; i++) mask_as_float[i] = ((unsigned char *)mask)[i];
         break;
      case 4:
         for(int i=0; i<nv; i++) mask_as_float[i] = ((short *)mask)[i];
         break;
      case 8:
         for(int i=0; i<nv; i++) mask_as_float[i] = ((int *)mask)[i];
         break;
      case 16:
         for(int i=0; i<nv; i++) mask_as_float[i] = ((float *)mask)[i];
         break;
   }

   for(int i=0; i<nv; i++) avg_mask[i] = mask_as_float[i];

   for(int k=0; k<nz; k++)
   {
      increment[k]=0;
      for(int v=0; v<np; v++)
      {
         if(mask_as_float[k*np+v] > 0.0) 
         {
            increment[k] = 1;
         }
      }
      nvec[k] += increment[k];
   }

   delete mask;

   for(int i=1; i<n; i++)
   {
/*
      if(opt_v)
      {
         printf("Reading mask %d: %s ...\n",i+1, maskList[ maskIndx->vec[i] ]);
      }
*/

      mask = read_nifti_image(maskList[ maskIndx->vec[i] ], &hdr);
      if(mask==NULL) 
      {
         return(NULL);
      }

      dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

/*
      if(opt_v)
      {
         printf("\tMatrix size = %d x %d x %d\n", hdr.dim[1], hdr.dim[2], hdr.dim[3]);
         printf("\tVoxe size = %f x %f x %f\n", hdr.pixdim[1], hdr.pixdim[2], hdr.pixdim[3]);
      }
*/

      if( nv != hdr.dim[1]*hdr.dim[2]*hdr.dim[3] )
      {
         printf("Error: Input masks must all have the same matrix size.\n"); 
         return(NULL);
      }

      switch(type) {
         case 2:
            for(int j=0; j<nv; j++) mask_as_float[j] = ((unsigned char *)mask)[j];
            break;
         case 4:
            for(int j=0; j<nv; j++) mask_as_float[j] = ((short *)mask)[j];
            break;
         case 8:
            for(int j=0; j<nv; j++) mask_as_float[j] = ((int *)mask)[j];
            break;
         case 16:
            for(int j=0; j<nv; j++) mask_as_float[j] = ((float *)mask)[j];
            break;
      }

      for(int j=0; j<nv; j++) avg_mask[j] += mask_as_float[j];

      for(int k=0; k<nz; k++)
      {
         increment[k]=0;
         for(int v=0; v<np; v++)
         {
            if(mask_as_float[k*np+v] > 0.0) 
            {
               increment[k] = 1;
            }
         }
         nvec[k] += increment[k];
      }

      delete mask;
   }

   // to ensure division by zero is avoided
   for(int k=0; k<nz; k++)
   if(nvec[k]==0) nvec[k]=1;

   for(int k=0; k<nz; k++) 
   {
      for(int v=0; v<np; v++) 
      {
         avg_mask[k*np + v] /= nvec[k];
      }
   }

   return(avg_mask);
}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   float thresh=50.0;
   char subjectImage[1024];
   char outputfilename[1024];
   char outputpath[1024];
   int number_of_input_masks; // number of input masks
   int n=0; // number of masks used for skull-stripping the subjectImage
   char **maskList;

   subjectImage[0]='\0';
   outputfilename[0]='\0';
   outputpath[0]='\0';

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'n':
            n = atoi(optarg);
            break;
         case 't':
            thresh = atof(optarg);
            break;
         case 'o':
            sprintf(outputfilename, "%s", optarg);
            break;
         case 's':
            sprintf(subjectImage, "%s", optarg);
            break;
         case 'v':
            opt_v = YES;
            break;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   if( subjectImage[0]=='\0' )
   {
      printf("Please specify the subject image using the -s argument.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("Subject image = %s\n", subjectImage);
   }
   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   if( outputfilename[0]=='\0' )
   {
      char sub_filename[512];

      niftiFilename(sub_filename, subjectImage);

      sprintf(outputfilename,"%s_brainwash.nii",sub_filename);
   }

   {
      char subject_dir[512];
      getDirectoryName(subjectImage, subject_dir);
      sprintf(outputpath,"%s/%s",subject_dir, outputfilename);
   }

   if(opt_v)
   {
      printf("Output file = %s\n", outputpath);
   }
   ////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("Thresholding at %4.1f%% level\n",thresh);
   }

   ////////////////////////////////////////////////////////////////////////////////
   number_of_input_masks = argc-optind;

   if(opt_v)
   {
      printf("Number of input masks = %d\n", number_of_input_masks);
   }

   if( number_of_input_masks==0 )
   {
      exit(0);
   }
   maskList = argv + optind;

   if(opt_v)
   {
      printf("Input masks:\n");
      for(int i=0; i<number_of_input_masks; i++)
      {
         printf("%d\t%s\n",i+1, maskList[i]);
      }
   }
   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   if(n<=0 || n>number_of_input_masks)
   {
      n = number_of_input_masks;
   }

   if(opt_v) 
   {
      if(n==number_of_input_masks)
      {
         printf("Using all %d masks for skull-stripping\n", n);
      }
      else
      {
         printf("Using %d randomly selected masks for skull-stripping\n", n);
      }
   }

   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   PERMUTATION maskIndx(number_of_input_masks);
   //maskIndx.permute();

   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   float *avg_mask;
   nifti_1_header hdr;
   int nx, ny, nz;
   float dx, dy, dz;
   int nv;
   int type;
   float *d;
   unsigned char *consensusmask;

   hdr = read_NIFTI_hdr(subjectImage);
   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];
   type=hdr.datatype;
   nv = nx*ny*nz;

   consensusmask = (unsigned char *)calloc(nv,1);
   d = (float *)calloc(n,sizeof(float));

   for(int i=0; i<n; i++) d[i]=1.0;

   avg_mask = compute_avg_mask(&maskIndx, n, maskList, d);
   if(avg_mask==NULL)  // update later to print an error message
   {
      exit(1);
   }

   for(int i=0; i<nv; i++)
   {
      if(avg_mask[i]<thresh)
      {
         consensusmask[i]=0;
      }
      else
      {
         consensusmask[i]=1;
      }
   }

   for(int iter=0; iter<0; iter++)
   {
      printf("iteration %d\n",iter+1);

      tt(&maskIndx, n, maskList, consensusmask, d);

      free(avg_mask);
      avg_mask = compute_avg_mask(&maskIndx, n, maskList, d);

      for(int i=0; i<nv; i++)
      {
         if(avg_mask[i]<thresh)
         {
            consensusmask[i]=0;
         }
         else
         {
            consensusmask[i]=1;
         }
      }
   }

   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   short *subim;
   nifti1_extender ext;
   FILE *fp;

/*
   if(opt_v)
   {
      printf("Reading subject image: %s ...\n", subjectImage);
   }
*/
   subim = (short *) read_nifti_image( subjectImage, &hdr);
   if(subim==NULL) // update later to print an error message
   {
      exit(1);
   }

   for(int i=0; i<nv; i++)
   {
      if(consensusmask[i]==0)
      {
         subim[i]=0;
      }
   }

   fp = fopen(subjectImage,"r");
   fread(&hdr, sizeof(nifti_1_header), 1, fp);
   // very important (fixed bug reported by Jay 4/21/2011)
   if(hdr.dim[0]<1 || hdr.dim[0]>7)
   {
      swapniftiheader(&hdr);
   }
   fread(&ext, sizeof(nifti1_extender), 1, fp);
   fclose(fp);

   sprintf(hdr.descrip,"Created by ART `brainwash' program.");
   hdr.vox_offset = 352; // important
   ext.extension[0] = 0; // important

   if(opt_v)
   {
      printf("Writing output image: %s ...\n", outputpath);
   }

   fp = fopen(outputpath,"w");
   fwrite(&hdr, sizeof(nifti_1_header), 1, fp);
   fwrite(&ext, sizeof(nifti1_extender), 1, fp);
   fwrite(subim, sizeof(short), nv, fp);
   fclose(fp);
   ////////////////////////////////////////////////////////////////////////////////
}
