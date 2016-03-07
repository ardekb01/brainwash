#include <stdlib.h>
#include <malloc.h>  
#include <math.h>
#include <strings.h>
#include <string.h>		// required by strlen()
#include <stdio.h>
#include <sys/types.h>  //      required by time(), stat()
#include <time.h>       //      required by time()
#include <sys/stat.h>   //      required by stat() 
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "spm_analyze.h"

#include <f2c.h>
#include <blaswrap.h>
extern "C" {
#include <clapack.h>
}

#include "babak_lib.h"
#include "niftiimage.h"

#define YES 1
#define NO 0

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

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
   NIFTIIMAGE im1, im2;
   int nv;
   float dx, dy, dz;
   short data_type;
   double sum;
   unsigned char *setA; 
   unsigned char *setB;

   im1.read( argv[1] );
   nv = im1.nv();
   data_type = im1.datatype();
   dx = im1.dx();
   dy = im1.dy();
   dz = im1.dz();

   setA = (unsigned char *)calloc(nv,1);

   if(data_type == 4 || data_type == 512)
   {
      short *image;
      image = (short *)im1.getdata();

      for(int i=0; i<nv; i++)
      {
         if(image[i]>0)
         {
            setA[i] = 1; 
         }
      }
   }

   if(data_type == 2)
   {
      char *image;
      image = im1.getdata();

      for(int i=0; i<nv; i++)
      if(image[i]>0)
      {
         setA[i] = 1; 
      }
   }

   im2.read( argv[2] );
   nv = im2.nv();
   data_type = im2.datatype();
   dx = im2.dx();
   dy = im2.dy();
   dz = im2.dz();

   setB = (unsigned char *)calloc(nv,1);

   if(data_type == 4 || data_type == 512)
   {
      short *image;
      image = (short *)im2.getdata();

      for(int i=0; i<nv; i++)
      {
         if(image[i]>0)
         {
            setB[i] = 1; 
         }
      }
   }

   if(data_type == 2)
   {
      char *image;
      image = im2.getdata();

      for(int i=0; i<nv; i++)
      if(image[i]>0)
      {
         setB[i] = 1; 
      }
   }

   float D, J;
   D = diceindex(setA, setB, nv);
   // J = jaccardindex(setA, setB, nv);
   J = D/(2.0-D);

   //printf("%f\t%f\n", D, J);
   printf("%f\n", J);
}
