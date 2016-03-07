/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \
*  brainwash2.cxx                                            *
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
   {"-output", 1, 'o'},
   {"-o", 1, 'o'},

   {"-v", 0, 'v'},
   {"-verbose", 0, 'v'},

   {"-n", 1, 'n'},

   {0, 0, 0}
};

char opt_v=NO;
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   char outputfile[1024];
   int number_of_curves; // number of input curves
   int n=0; // number points on the curve
   float *dif;

   char **curveList;

   outputfile[0]='\0';

   while( (opt=getoption(argc, argv, options)) != -1)
   {
      switch (opt) {
         case 'o':
            sprintf(outputfile, "%s", optarg);
            break;
         case 'v':
            opt_v = YES;
            break;
         case 'n':
            n = atoi(optarg);
            break;
      }
   }

   ////////////////////////////////////////////////////////////////////////////////
   if( outputfile[0]=='\0' )
   {
      printf("Please specify an output file using the -o argument.\n");
      exit(0);
   }

   if(opt_v)
   {
      printf("Output file = %s\n", outputfile);
   }
   ////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////////////
   number_of_curves = argc-optind;

   if(opt_v)
   {
      printf("Number of input curves = %d\n", number_of_curves);
   }

   if( number_of_curves==0 )
   {
      exit(0);
   }
   curveList = argv + optind;

   if(opt_v)
   {
      printf("Input curves:\n");
      for(int i=0; i<number_of_curves; i++)
      {
         printf("%d\t%s\n",i+1, curveList[i]);
      }
   }
   ////////////////////////////////////////////////////////////////////////////////

   if(n<=0)
   {
      exit(0);
   }

   if(opt_v)
   {
      printf("Number of points = %d\n", n);
   }

   ////////////////////////////////////////////////////////////////////////////////
   float *avgD, *avgJ, *sdD, *sdJ, D, J;
   FILE *fpin, *fpout;

   avgD = (float *)calloc(n, sizeof(float));
   avgJ = (float *)calloc(n, sizeof(float));
   sdD = (float *)calloc(n, sizeof(float));
   sdJ = (float *)calloc(n, sizeof(float));
   dif = (float *)calloc(n, sizeof(float));

   for(int i=0; i<n; i++) avgD[i]=avgJ[i]=0.0;

   for(int i=0; i<n; i++) sdJ[i]=sdD[i]=0.0;

   for(int c=0; c<number_of_curves; c++)
   {
      fpin = fopen(curveList[c],"r");
      for(int i=0; i<n; i++)
      {
         fscanf(fpin,"%f\t%f\n",&D,&J);
         avgD[i] += D; 
         avgJ[i] += J; 
      }
      fclose(fpin);
   }


   for(int i=0; i<n; i++) 
   {
      avgD[i] /= number_of_curves;
      avgJ[i] /= number_of_curves;
   }

   for(int c=0; c<number_of_curves; c++)
   {
      fpin = fopen(curveList[c],"r");
      for(int i=0; i<n; i++)
      {
         fscanf(fpin,"%f\t%f\n",&D,&J);
         sdD[i] += (avgD[i]-D)*(avgD[i]-D); 
         sdJ[i] += (avgJ[i]-J)*(avgJ[i]-J); 
      }
      fclose(fpin);
   }

   for(int i=0; i<n; i++) 
   {
      sdD[i] /= (number_of_curves-1.0);
      sdJ[i] /= (number_of_curves-1.0);

      sdD[i] = sqrtf(sdD[i]);

      sdJ[i] = sqrtf(sdJ[i]);
   }

   for(int i=0; i<n; i++) 
   {
      printf("%d\t%f\t%f\t%f\t%f\n", i+1, avgD[i], avgJ[i], sdD[i], sdJ[i]);
   }

/*
   fpout=fopen(outputfile,"w");
   for(int c=0; c<number_of_curves; c++)
   {
      fpin = fopen(curveList[c],"r");
      for(int i=0; i<n; i++)
      {
         fscanf(fpin,"%f\n",&dum);
         fprintf(fpout,"%f\t",dum);
      }
      fclose(fpin);
      fprintf(fpout,"\n");
   }
   fclose(fpout);

   printf("Difference:\n");
   for(int i=0; i<n; i++)
   {
      dif[i] = avgD[n-1]-avgD[i];
      printf("%d %d\n",i+1, (int)(dif[i]*1000.0));
   }
*/
}
