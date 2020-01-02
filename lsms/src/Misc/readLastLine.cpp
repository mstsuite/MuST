#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readLastLine.hpp"


// return the size of the file. This does not change the current file position
long getFileSize(FILE *fp)
{
  long fsize = 0;
  long fpos=ftell(fp);

  fseek(fp,0,SEEK_END);
  fsize = ftell(fp); 
  fseek(fp,fpos,SEEK_SET);//reset stream position!!

  return fsize;
}

// read the last line in the file, returning the start in the buffer 
int readLastLine(FILE *fp, char *buffer, int bufferLength)
{
  int i;

  int size=bufferLength;
  long fileSize= getFileSize(fp);
  if(fileSize<size) size=fileSize;

  if(fseek(fp, -size, SEEK_END))
  {
    perror("cannot seek");
    exit(1);
  }

  size=fread(buffer, sizeof(char), bufferLength, fp);
  buffer[size] = '\0';
  i=size-1;
  if(buffer[i]=='\n'){
    buffer[i] = '\0';
  }
  while(i >=0 && buffer[i] != '\n')
    --i;
  ++i;
  return i;
}

int readNextIterationNumber(const char *fname)
{
  char buffer[256];
  int lastIterationNumber,pos;

  FILE *fp=fopen(fname,"r");
  if(fp==NULL) return 0;
  pos=readLastLine(fp,buffer,255);
  fclose(fp);
  // printf("last line: '%s'\n",&buffer[pos]); exit(0);
  sscanf(&buffer[pos],"%d",&lastIterationNumber);
  return lastIterationNumber+1;
}
