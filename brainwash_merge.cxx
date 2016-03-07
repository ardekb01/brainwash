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
   
int number_of_input_masks; // number of input masks

//////////////////////////////////////////////////////////////////////////////////////////////////
int opt;

static struct option options[] =
{
   {"-rat",0,'r'},

   {"-version",0,'V'},
   {"-Version",0,'V'},
   {"-V",0,'V'},

   {"-n", 1, 'n'},

   {"-T", 1, 't'},
   {"-t", 1, 't'},
   {"-threshold", 1, 't'},
   {"-thres", 1, 't'},

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
char opt_rat=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

float *compute_avg_mask(PERMUTATION *maskIndx, int n, char **maskList)
{
   nifti_1_header hdr;
   int nv;
   unsigned char *mask;
   float *avg_mask;

   if(opt_v && n<number_of_input_masks)
   {
      printf("Masks used for skull-stripping:\n");
      for(int i=0; i<n; i++)
      {
         printf("%d\t%s\n", i+1, maskList[ maskIndx->vec[i] ] );
      }
   }

   for(int i=0; i<n; i++)
   {
      if(opt_v)
      {
         printf("Reading mask %d: %s ...\n",i+1, maskList[ maskIndx->vec[i] ]);
      }

      mask = (unsigned char *)read_nifti_image(maskList[ maskIndx->vec[i] ], &hdr);
      if(mask==NULL) 
      {
         printf("\nError reading %s\n\n",maskList[ maskIndx->vec[i] ]);
         exit(1);
      }

      if(i==0)
      {
         nv = hdr.dim[1] * hdr.dim[2] * hdr.dim[3];

         avg_mask = (float *)calloc(nv, sizeof(float));
         if(avg_mask==NULL) 
         {
            printf("\nMemory allocation error (avg_mask)\n\n");
            exit(1);
         }

         for(int v=0; v<nv; v++) avg_mask[v]=0.0;
      }


      for(int v=0; v<nv; v++) avg_mask[v] += mask[v];

      delete mask;
   }

   for(int v=0; v<nv; v++) avg_mask[v] /= n;

   return(avg_mask);
}

float *compute_avg_mask2(PERMUTATION *maskIndx, int n, char **maskList)
{
   int *nvec;
   int *increment;
   nifti_1_header hdr;
   int nx, ny, nz;
   float dx, dy, dz;
   int np, nv;
   unsigned char *mask;
   float *avg_mask;
   float *mask_as_float;

   if(opt_v && n<number_of_input_masks)
   {
      printf("Masks used for skull-stripping:\n");
      for(int i=0; i<n; i++)
      {
         printf("%d\t%s\n", i+1, maskList[ maskIndx->vec[i] ] );
      }
   }

   mask = (unsigned char *)read_nifti_image(maskList[ maskIndx->vec[0] ], &hdr);
   if(mask==NULL)
   {
      printf("Error reading %s, aborting ...\n", maskList[ maskIndx->vec[0] ]);
      exit(1);
   }

   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];

   nv = nx*ny*nz;
   np = nx*ny;
   avg_mask = (float *)calloc(nv, sizeof(float));
   if(avg_mask == NULL)
   {
      printf("\nMemory allocation error (avg_mask), aborting ...\n\n");
      exit(1);
   }

   mask_as_float = (float *)calloc(nv, sizeof(float));
   if(mask_as_float == NULL)
   {
      printf("\nMemory allocation error (mask_as_float), aborting ...\n\n");
      exit(1);
   }

   nvec = (int *)calloc(nz, sizeof(int));
   increment = (int *)calloc(nz, sizeof(int));
   for(int k=0; k<nz; k++) nvec[k]=0;

   for(int v=0; v<nv; v++) 
   {
      mask_as_float[v] = mask[v];
      avg_mask[v] = mask_as_float[v];
   }

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

      mask = (unsigned char *)read_nifti_image(maskList[ maskIndx->vec[i] ], &hdr);
      if(mask==NULL)
      {
         printf("Error reading %s, aborting ...\n", maskList[ maskIndx->vec[i] ]);
         exit(1);
      }

      for(int v=0; v<nv; v++) 
      {
         mask_as_float[v] = mask[v];
         avg_mask[v] += mask_as_float[v];
      }

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

//////////////////////////////////////////////////////////////////////////////////////////////////
void print_help_and_exit()
{
   printf("\nUsage: brainwash_merge [-verbose -version -rat -o <filename>] -sub <subjectID.nii> <masks ...>" 
   "\n\n"
 
   "-sub/-s <subjectID.nii>: specifies the MRI volume to be skull-stripped.\n\n"

   "<mask ...> are a number of mask volumes obtained using the `brainwash' program.\n"
   "These are typically named as: <subjectID_atlasID_mask.nii>\n\n"

   "-verbose/-v: enables verbose mode.\n\n"

   "-version/-V: prints software version.\n\n"

   "-rat: required for skull-stripping rat images.\n\n"

   "-o <filename>: saves the output (skull-stripped <subjectID.nii>) in this file\n"
   "(default=<subjectID>_brainwash.nii).\n\n"
   );

   exit(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////////////////

void check_masks(char **maskList, int number_of_masks, DIM subdim)
{
   nifti_1_header hdr;

   for(int i=0; i<number_of_masks; i++)
   {
      if( read_NIFTI_hdr(maskList[i], &hdr) == 0 ) 
      {
         exit(0);
      }

      if( hdr.datatype != DT_UNSIGNED_CHAR) 
      {
         printf("Sorry, but this program requires %s to be of datatype 2 (i.e., DT_UNSIGNED_CHAR).\n"
         "The current datatype is %d.\n", maskList[i], hdr.datatype);
         exit(0);
      }

      if(hdr.dim[1]!=subdim.nx || hdr.dim[2]!=subdim.ny || hdr.dim[3]!=subdim.nz)
      {
         printf("Mask %s has incompatible matrix dimensions.\n", maskList[i]);
         exit(0);
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   float icv=0.0;

   nifti_1_header hdr;
   nifti1_extender ext;

   int nx, ny, nz;
   float dx, dy, dz;
   int nv;

   short *subim;

   DIM subdim;

   char subImagePath[512]="";
   char outputfilename[512]="";

   float thresh=50.0;
   char outputpath[1024];
   int n=0; // number of masks used for skull-stripping the subjectImage
   char **maskList;

   /////////////////////////////////////////////////////////////////////
   if(argc==1)
   {
      print_help_and_exit();
   }
   /////////////////////////////////////////////////////////////////////

   outputpath[0]='\0';

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'V':
            printf("Version: July 28, 2011\n");
            exit(0);
         case 'r':
            opt_rat = YES;
            break;
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
            sprintf(subImagePath, "%s", optarg);
            break;
         case 'v':
            opt_v = YES;
            break;
      }
   }

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

   ////////////////////////////////////////////////////////////////////////////////
   if( outputfilename[0]=='\0' )
   {
      char subjectID[512];

      niftiFilename(subjectID, subImagePath);

      sprintf(outputfilename,"%s_brainwash.nii",subjectID);
   }

   {
      char subject_dir[512];
      getDirectoryName(subImagePath, subject_dir);
      sprintf(outputpath,"%s/%s",subject_dir, outputfilename);
   }

   if(opt_v)
   {
      printf("Output (skull-stripped) volume = %s\n", outputpath);
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


   if(n<=0 || n>number_of_input_masks)
   {
      n = number_of_input_masks;
   }

   if(opt_v && n<number_of_input_masks) 
   {
      printf("Using %d randomly selected masks for skull-stripping\n", n);
   }

   ////////////////////////////////////////////////////////////////////////////////

   check_masks(maskList, number_of_input_masks, subdim);

   ////////////////////////////////////////////////////////////////////////////////
   float *avg_mask;

   PERMUTATION maskIndx(number_of_input_masks);

   maskIndx.permute();

   if(!opt_rat)
   {
      avg_mask = compute_avg_mask(&maskIndx, n, maskList);
   }
   else
   {
      avg_mask = compute_avg_mask2(&maskIndx, n, maskList);
   }
   if(avg_mask==NULL)  // update later to print an error message
   {
      exit(1);
   }
   ////////////////////////////////////////////////////////////////////////////////

   if(opt_v)
   {
      printf("Reading subject image: %s ...\n", subImagePath);
   }
   subim = (short *) read_nifti_image( subImagePath, &hdr);
   if(subim==NULL) // update later to print an error message
   {
      printf("\nError reading %s, aborting ...\n\n",subImagePath);
      exit(1);
   }
   nx=hdr.dim[1]; ny=hdr.dim[2]; nz=hdr.dim[3];
   dx=hdr.pixdim[1]; dy=hdr.pixdim[2]; dz=hdr.pixdim[3];
   nv = nx*ny*nz;

   for(int i=0; i<nv; i++)
   {
      if(avg_mask[i]<thresh)
      {
         subim[i]=0;
      }
      else
      {
         icv += 1.0;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   FILE *fp;

   if( read_NIFTI_hdr(subImagePath, &hdr) == 0 ) 
   {
      exit(0);
   }

   sprintf(hdr.descrip,"Created by ART `brainwash' program.");
   hdr.vox_offset = 352; // important
   ext.extension[0] = 0; // important

   if(opt_v)
   {
      printf("Writing output image: %s ...\n", outputpath);
   }

   fp = fopen(outputpath,"w");
   if(fp==NULL)
   {
      printf("\nError writing to %s, aborting ...\n\n",outputpath);
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

   if( fwrite(subim, sizeof(short), nv, fp) != nv )
   {
      printf("\nError writing voxel data to %s!\n", outputpath);
      exit(1);
   }

   fclose(fp);

   if(opt_v)
   {
      printf("Final brain mask volume = %d mm^3\n", (int)(icv*dx*dy*dz + 0.5));
   }

   free(subim);
   free(avg_mask);
}
