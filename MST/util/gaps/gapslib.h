/** gapslib.h 
 * 
 * Useful macros, inline functions and structs. 
 * 
 * George Bargoud: <gbargoud@gmail.com>
 **/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#ifndef GAPS_LIB_H
#define GAPS_LIB_H

/*********\
* Structs *
\*********/
typedef struct 
{
	double x,y,z;
} vector_t;

typedef struct
{
	/* cartesian coordinates */
	double x,y,z;
	/* coordinates using a,b,c as the unit vectors */
	double Ca,Cb,Cc;
	/* Stores atomic number for real particles and radius for virtual ones */
	double rad;
	int Z;
} atom_t;

/* the type for the crystal */
typedef struct 
{
	vector_t a,b,c;
	atom_t * atoms;
	/* atomc is the count of atoms actually in the array and atomsiz is the
	size of the array itself. */
	int atomc,atomsiz;
} crystal_t;

/* fill_t is the type for the filling algorithms. */
typedef void (*fill_t)(crystal_t * crystal);

/* A redefinition of an atom. Used when guessing widths. */
typedef struct
{
	double A;
	int Z;
} redef_t;

/********************\
* External variables *
\********************/

/* The type of the virtual atoms */
extern int type;
/* Transformation matrix for changing from cartesian to the bounding vectors */
extern double transform[3][3];
/* Atomic radii */
extern double radii[];
/* If true, print warnings */
extern int warnings;

/********\
* Macros *
\********/
/* Define NAN if it is not defined */
#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* Define INFINITY if it is not defined */
#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

/* only behaves as a print if compiled with "-D debug" as an option */
#ifdef debug
#define dbg_printf(...) fprintf(stderr,__VA_ARGS__)
#else
#define dbg_printf(...)
#endif

/* print this message to stderr and then exit with a failure status */
#define fail(...) do {\
	fprintf(stderr,__VA_ARGS__);\
	fprintf(stderr,"\n");\
	exit(EXIT_FAILURE);\
} while(0)

/* get the part of i after the decimal point */
#define dec(i) ((i-((int)i))+(i<0?1:0))

/* get the magnitude of a vector */
#define MAGNITUDE(v) sqrt(SQMAG(v))

/* get the square of the magnitude of v */
#define SQMAG(v) DOT(v,v)

/* the dot product of u and v */
#define DOT(u,v) (u.x*v.x+u.y*v.y+u.z*v.z)

/* Check the MPI error code and fail if it is wrong. */
#define CHECK(rc) if (rc != MPI_SUCCESS)\
{\
	MPI_Abort(MPI_COMM_WORLD,rc);\
}

/* Macro for allocating enough space to copy a string into. */
#define stralloc(str) calloc(strlen(str)+1,sizeof(char))

/* Print if the warning flag is 1 */
#define warn_printf(...) if (warnings) fprintf(stderr,__VA_ARGS__)

/******************\
* Helper functions *
\******************/
/** Extend the list of atoms to all adjacent crystal blocks. */
crystal_t * extend(crystal_t crystal);
/** Append the extended version of atom to the crystal. buf is how close to the
  * edge it has to be to be appended. */
void appendall(crystal_t * crystal, atom_t atom);
/** Return a copy of atom moved along vector by distance */
inline atom_t move(atom_t atom,vector_t vector,double distance);
/** Get the distance between the outermost edges of the atoms. */
double dist(atom_t one,atom_t two);
/** Get the radius of atom in Angstroms */
inline double radius(atom_t atom);
/** Insert the atom in crystal's atoms array. Expands the array if needed. */
void append(crystal_t * crystal,atom_t atom);
/** Get the distance to the closest atom in the crystal to atom. */
double closest(crystal_t crystal,atom_t atom);
/** Get the distance to the closest atom as measured from the center of the 
  * atoms */
double closest_c(crystal_t crystal,atom_t atom);
/** Finds the vectors describing the voronoi polyhedron around an atom */
vector_t * voronoi(crystal_t crystal,atom_t atom,int * c);
/** Gets the vertices of the voronoi polyhedron */
vector_t * vertices(crystal_t crystal,atom_t atom,int * c);
/** Gets the radius of the circumscribed sphere of the Voronoi polyhedron */
double Rcc(crystal_t crystal,atom_t atom);
/** Gets the radius of the inscribed sphere of the Voronoi polyhedron */
double Ric(crystal_t crystal,atom_t atom);
/** The ratio of Rcc and closest_c */
double r_c_ratio(crystal_t crystal,atom_t atom);
/** Cross product */
inline vector_t cross(vector_t one,vector_t two);
/** Fill the Ca, Cb and Cc fields in atom. 
  * Returns 0 if in bounds and 1 if out of bounds. */
inline int get_coefficients(atom_t * atom);
/** Bring atom in the bounds */
void bringin(atom_t * atom,crystal_t crys);
/** return the gradient of closest at the point where atom is. */
vector_t grad(atom_t atom, crystal_t crystal,double (*func)(crystal_t,atom_t),
	double epsilon);
/** find a nearby spot which is better */
atom_t drift(atom_t atom,crystal_t crystal,double (*func)(crystal_t,atom_t),
	double d,double dmin);
/** Make the types for passing atoms, vectors and crystals through MPI. */
inline void make_mpi_types(MPI_Datatype newtypes[3]);
/** get the floor of the binary logarithm of n */
int floor_log2(int i);
#endif /* GAPS_LIB_H */
