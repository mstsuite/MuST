/** gaps.h 
 *
 * George Bargoud: <gbargoud@gmail.com>
 **/
#include "gapslib.h"

#ifndef GAPS_H
#define GAPS_H

/******************\
* global variables *
\******************/
/* atomic radii for most elements */
extern double radii[];
/* symbols for most elements */
extern char * symbols[];
/* conversion factor */
double unit = 1.00;
/* Will be the transformation matrix to go from cartesian to a basis of the 
 * crystal's bounding vectors. */
double transform[][3] = {{1.0,0.0,0.0},
						 {0.0,1.0,0.0},
						 {0.0,0.0,1.0}};
extern struct {
	char * name;
	fill_t func;
} funcmap[];

/* Set in main or the conf file. */
/* Input and output files */
char * inFileName = NULL;
char * outFileName = NULL;
/* The filling algorithm */
fill_t fill;

/* whether to guess the radii */
int g = 0;
/* The type for virtual atoms */
int type = -1;
/* direct coordinates on input and output */
int direct = 0;
int direct_out = 0;
/* whether to print warning messages */
int warnings = 0;

/***********\
* functions *
\***********/
int main(int argc,char** argv);
/** Print help text and exit. */
void help(void);
/** redefine an atomic radius. */
void redefine(char* s);
/** open filename and read the crystal. If filename is null, it reads stdin. */
crystal_t parse(char* filename);
/** output crystal to filename. If filename is null, it outputs to stdout. */
void out(char * filename, crystal_t crystal, redef_t * redefs, int count);
/** takes the algorithm name and returns a pointer to the one that will be 
  * used  */
fill_t funcfromname(char * name);
/** Sets the transform matrix */
void set_transform(crystal_t * crystal);
/** Changed floating point numbers from fortran to C syntax */
static inline void ftoc(char * s);
/** A helper unction for Gauss-Jordan elimination */
void eliminate(double matrix[][3],int row,int col);
/** Set units. takes a string which is either "ang", "au" or a floating point 
  * conversion factor from angstroms. */
void units(char * s);
/** Parse a configuration file */
void parse_conf(char * filename);
/** Given a symbol, look up the atomic number */
int lookup(const char * symbol);
/** Guess proper radii for the atoms based on their distances and using the 
  * preloaded radii as a guide. */
redef_t * guess(crystal_t crystal,int * count);

#endif /* GAPS_H */
