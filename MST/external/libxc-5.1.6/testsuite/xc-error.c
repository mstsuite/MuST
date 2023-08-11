/*
 Copyright (C) 2014 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Buffer size */
#define BUFSIZE 4096

/* Maximum number of columns in data */
#define MAXCOL 100
/* Legend entry length */
#define LEGLEN 20

#define error_exit() fclose(in);		\
  fclose(ref);					\
  exit(1);

#define in_line()   cp=fgets(buf,BUFSIZE,in);			\
  if(cp!=buf) {							\
  fprintf(stderr,"Error reading line from input file.\n");	\
  error_exit();							\
  }  

#define ref_line()   cp=fgets(buf,BUFSIZE,ref);			\
  if(cp!=buf) {							\
  fprintf(stderr,"Error reading line from reference file.\n");	\
  error_exit();							\
  }  

double maxabs(double x, double y) {
  return fabs(x)>fabs(y) ? fabs(x) : fabs(y);
}

double error(double x, double y) {
  return fabs(x-y)/(1.0+maxabs(x,y));
}

int main(int argc, char **argv) {
  /* Input file and reference */
  FILE *in;
  FILE *ref;

  /* Functional IDs */
  int fidin, fidref;

  /* Sizes of input and reference */
  int nin, nref;

  /* Number of columns read in */
  int cin, cref;

  /* Input buffer */
  char buf[BUFSIZE];
  char *cp;
  int cur, nread;

  /* Loop indices */
  int i, j;
  
  /* Input and reference data */
  double din[MAXCOL], dref[MAXCOL];
  /* Column legends */
  char legin[MAXCOL][LEGLEN], legref[MAXCOL][LEGLEN];

  /* Maximum difference between input and output */
  double maxdiff[MAXCOL];
  double l_err;

  /* Input tolerance */
  double tol;
  /* Verbose operation? */
  int verbose;

  if(argc!=4 && argc!=5) {
    printf("Usage: %s file reference tolerance (verbose)\n",argv[0]);
    return 1;
  }

  /* Open files */
  in=fopen(argv[1],"r");
  if(!in) {
    fprintf(stderr,"Error opening input file.\n");
    exit(1);
  }

  ref=fopen(argv[2],"r");
  if(!ref) {
    fprintf(stderr,"Error opening reference file.\n");
    exit(1);
  }

  /* Set tolerance and verbosity */
  tol=atof(argv[3]);
  verbose=(argc==5) ? atoi(argv[4]) : 0;

  /* Read first line: functional id and file length */
  in_line();
  if(sscanf(buf,"%i %i",&fidin,&nin)!=2) {
    fprintf(stderr,"Error reading func_id and file size from input file.\n");
    error_exit();
  }

  ref_line();
  if(sscanf(buf,"%i %i",&fidref,&nref)!=2) {
    fprintf(stderr,"Error reading func_id and file size from reference file.\n");
    error_exit();
  }
  
  if(fidin!=fidref) {
    fprintf(stderr,"Functional ids %i and %i don't match!\n",fidin,fidref);
    error_exit();
  }

  if(nin!=nref) {
    fprintf(stderr,"Sizes of files %i and %i don't match!\n",nin,nref);
    error_exit();
  }

  /* Read in legends */
  in_line();
  cin=0;
  cur=0;
  while(sscanf(buf+cur,"%s%n",legin[cin],&nread)==1) {
    cin++;
    cur+=nread;

    if(cin==MAXCOL) {
      fprintf(stderr,"Array overflow. Increase MAXCOL.\n");
      error_exit();
    }
  }
  
  ref_line();
  cref=0;
  cur=0;
  while(sscanf(buf+cur,"%s%n",legref[cref],&nread)==1) {
    cref++;
    cur+=nread;

    if(cref==MAXCOL) {
      fprintf(stderr,"Array overflow. Increase MAXCOL.\n");
      error_exit();
    }
  }

  /* Compare legends */
  if(cin != cref) {
    fprintf(stderr,"Number of columns doesn't match: input %i, reference %i.\n",cin,cref);
    error_exit();
  }
  for(i=0;i<cin;i++)
    if(strcmp(legin[i],legref[i])) {
      fprintf(stderr,"Legends of column %i don't match: %s vs %s\n",i,legin[i],legref[i]);
      error_exit();
    }

  /* Initialize difference data */
  for(i=0;i<MAXCOL;i++)
    maxdiff[i]=0.0;
    
  /* Read in data */
  for(i=0;i<nin;i++) {
    static const char fmt[]="%lf%n";

    /* Input line */
    in_line();
    cur=0;
    j=0;
    while(sscanf(buf+cur,fmt,&din[j],&nread)==1) {
      j++;
      cur+=nread;
      
      if(j==MAXCOL) {
        fprintf(stderr,"Array overflow. Increase MAXCOL.\n");
	      error_exit();
      }	
    }

    /* Reference line */
    ref_line();
    cur=0;
    j=0;
    while(sscanf(buf+cur,fmt,&dref[j],&nread)==1) {
      j++;
      cur+=nread;
      
      if(j==MAXCOL) {
	      fprintf(stderr,"Array overflow. Increase MAXCOL.\n");
	      error_exit();
      }
    }

    /* Compute error */
    for(j=0;j<cin;j++){
      l_err = error(din[j], dref[j]);
      if(l_err > maxdiff[j])
	      maxdiff[j] = l_err;
      if(l_err > tol){
        fprintf(stderr, "\n%i %i %14.10le %14.10le %le\n", 
                i+2, j, din[j], dref[j], error(din[j], dref[j]));
      }
    }

  }

  fclose(in);
  fclose(ref);

  if(verbose) {
    /* Verbose operation */
    for(i=0;i<cin;i++)
      printf(" %13s",legin[i]);
    printf("\n");
    for(i=0;i<cin;i++)
      printf(" % e",maxdiff[i]);
    printf("\n");

  } else {
    /* Silent operation */
    double max=0.0;
    for(j=0;j<cin;j++)
      if(maxdiff[j]>max)
      	max=maxdiff[j];

    printf("%i\n",max<=tol);
  }

  return 0;
}
