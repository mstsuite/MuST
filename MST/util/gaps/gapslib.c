/** gapslib.c 
 * 
 * helper functions for fill.c and gaps.c 
 * 
 * George Bargoud: <gbargoud@gmail.com>
 **/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "gapslib.h"

#define TOL 1.0e-15

/** Extend the list of atoms into the neighboring shapes, any atoms close to 
  * edge are copied over */
crystal_t * extend(crystal_t crystal)
{
	crystal_t * new;
	int i;
	
	new = malloc(sizeof(crystal_t));
	if (new == NULL)
		fail("malloc failed");

	/* copy over the fields of crystal. */
	new->a = crystal.a;
	new->b = crystal.b;
	new->c = crystal.c;

	new->atoms = calloc(crystal.atomsiz,sizeof(atom_t));
	if (new->atoms == NULL)
		fail("calloc failed");
	new->atomsiz = crystal.atomsiz;
	new->atomc = 0;
	

	for (i=0;i<crystal.atomc;i++)
	{
		appendall(new,crystal.atoms[i]);
	}
	return new;
}

/** Append the atom and it's nearby versions in neighboring cells to the 
  * crystal. */
void appendall(crystal_t * crystal, atom_t atom)
{
/* Check if two atoms are equal */
#define EQ(atom1,atom2) ((fabs(atom1.x - atom2.x) < TOL) &&\
						 (fabs(atom1.y - atom2.y) < TOL) &&\
						 (fabs(atom1.z - atom2.z) < TOL))

	atom_t atom1,atom2,atom3;
	double rad = (radius(atom)+radii[type]);
	double buf = rad/(MAGNITUDE(crystal->a));
	append(crystal,atom);
	
	/* Check in the a direction. atom1 is the original atom after being moved
	 * in that direction */
	if (atom.Ca < buf)
	{
		atom1=move(atom,crystal->a,1.0);
		append(crystal,atom1);
	}
	else if (atom.Ca > (1-buf))
	{
		atom1=move(atom,crystal->a,-1.0);
		append(crystal,atom1);
	}
	else
	{
		atom1 = atom;
	}
	
	/* Similarly, check the b direction. atom2 is the original atom moved in 
	 * that direction and atom3 is atom1 moved in the b direction. */

	buf = rad/(MAGNITUDE(crystal->b));
	if (atom.Cb < buf)
	{
		atom2=move(atom,crystal->b,1.0);
		append(crystal,atom2);
		if (!EQ(atom,atom1))
		{
			atom3=move(atom1,crystal->b,1.0);
			append(crystal,atom3);
		}
		else
		{
			atom3 = atom;
		}
	}
	else if (atom.Cb > (1-buf))
	{
		atom2=move(atom,crystal->b,-1.0);
		append(crystal,atom2);

		if (!EQ(atom,atom1))
		{
			atom3=move(atom1,crystal->b,-1.0);
			append(crystal,atom3);
		}
		else
		{
			atom3 = atom;
		}
	}
	else
	{
		atom2 = atom;
		atom3 = atom;
	}
	
	/* Now all four of atom, atom1, atom2 and atom3 are moved in the c 
	 * direction if they exist and can be moved. This concludes all of the
	 * possible movements. */
	buf = rad/(MAGNITUDE(crystal->c));
	if (atom.Cc < buf)
	{
		atom_t tmp;
		
		tmp = move(atom,crystal->c,1.0);
		append(crystal,tmp);

		if (!EQ(atom,atom1))
		{
			tmp = move(atom1,crystal->c,1.0);
			append(crystal,tmp);
		}
		if (!EQ(atom,atom2))
		{
			tmp = move(atom2,crystal->c,1.0);
			append(crystal,tmp);
		}
		if (!EQ(atom,atom3))
		{
			tmp = move(atom3,crystal->c,1.0);
			append(crystal,tmp);
		}
	}
	else if (atom.Cc > (1-buf))
	{
		atom_t tmp;
		tmp = move(atom,crystal->c,-1.0);
		append(crystal,tmp);
		if (!EQ(atom,atom1))
		{
			tmp = move(atom1,crystal->c,-1.0);
			append(crystal,tmp);
		}
		if (!EQ(atom,atom2))
		{
			tmp = move(atom2,crystal->c,-1.0);
			append(crystal,tmp);
		}
		if (!EQ(atom,atom3))
		{
			tmp = move(atom3,crystal->c,-1.0);
			append(crystal,tmp);
		}
	}
#undef EQ
}

/** Move the atom. */
inline atom_t move(atom_t atom,vector_t vector,double distance)
{
	atom.x += (distance*vector.x);
	atom.y += (distance*vector.y);
	atom.z += (distance*vector.z);

	return atom;
}

/** Get the distance from the edge of one atom to the edge of another. */
double dist(atom_t one, atom_t two)
{
	double dx,dy,dz,centerdist;

	dx = one.x-two.x;
	dy = one.y-two.y;
	dz = one.z-two.z;
	
	/* get the distance between the centers */
	centerdist = sqrt((dx*dx)+(dy*dy)+(dz*dz));
	
	return centerdist;
}

/** Get the radius of an atom.
  * checks the value of Z, if it's a positive value, it is used as the atomic 
  * number. If it is a negative value then it is a virtual particle. */
inline double radius(atom_t atom)
{
	if (atom.Z == 0)
		return atom.rad;
	else if (atom.Z < 0)
		return radii[-atom.Z];
	else
		return radii[atom.Z];
}

/** Append an atom to the crystal, if the array is too small, it doubles it. */
void append(crystal_t * crystal,atom_t atom)
{
	if ((fabs(atom.x - atom.x) > TOL) || (fabs(atom.y - atom.y) > TOL) || (fabs(atom.z - atom.z) > TOL))
	{
		return;
	}
	if ((crystal->atomc) >= (crystal->atomsiz))
	{
		crystal->atomsiz *= 2;
		crystal->atoms=realloc(crystal->atoms,crystal->atomsiz*sizeof(atom_t));
		if (crystal->atoms == NULL)
		{
			fail("realloc error");
		}
	}
	(crystal->atoms)[crystal->atomc] = atom;
	(crystal->atomc)++;
}

/** Fill the Ca, Cb and Cc fields in atom. Return whether it is in bounds. */
inline int get_coefficients(atom_t * atom)
{
	atom->Ca = atom->x*transform[0][0] +
			   atom->y*transform[0][1] +
			   atom->z*transform[0][2];

	atom->Cb = atom->x*transform[1][0]+
			   atom->y*transform[1][1]+
			   atom->z*transform[1][2];

	atom->Cc = atom->x*transform[2][0]+
			   atom->y*transform[2][1]+
			   atom->z*transform[2][2];
	
	return ((atom->Ca < -0.0) || (atom->Ca > 1.0) ||
			(atom->Cb < -0.0) || (atom->Cb > 1.0) ||
			(atom->Cc < -0.0) || (atom->Cc > 1.0));
}

/* Bring atom into the bounds */
void bringin(atom_t * atom, crystal_t crys)
{
	dbg_printf("atom at %lf,%lf,%lf\n",atom->x,atom->y,atom->z);

	atom->Ca = dec(atom->Ca);
	atom->Cb = dec(atom->Cb);
	atom->Cc = dec(atom->Cc);

	atom->x = atom->Ca*crys.a.x+atom->Cb*crys.b.x+atom->Cc*crys.c.x;
	atom->y = atom->Ca*crys.a.y+atom->Cb*crys.b.y+atom->Cc*crys.c.y;
	atom->z = atom->Ca*crys.a.z+atom->Cb*crys.b.z+atom->Cc*crys.c.z;
	
	dbg_printf("moved to %lf,%lf,%lf\n",atom->x,atom->y,atom->z);
}

/** returns the distance to the closest atom. */
double closest(crystal_t crystal, atom_t atom)
{
	int i;
	double closest = 1.0/0.0;
	
	for(i=0;i<crystal.atomc;i++)
	{
		double d = dist(crystal.atoms[i],atom);
		d -= radius(crystal.atoms[i]);
		if (d < closest)
		{
			closest = d;
		}
	}
	closest -= radius(atom);
	return closest;
}

/** The closest distance as measured from the centers of the atoms */
double closest_c(crystal_t crystal, atom_t atom)
{
	int i;
	double closest = 1.0/0.0;
	
	for(i=0;i<crystal.atomc;i++)
	{
		double d = dist(crystal.atoms[i],atom);
		if (d < closest)
		{
			closest = d;
		}
	}
	return closest;
}

/** Gets the voronoi polyhedron around atom. Stores the number of vectors in *c
  */
vector_t * voronoi(crystal_t crystal,atom_t atom,int * c)
{
#define EQ(one,two) ((fabs(one.x - two.x) < TOL) && (fabs(one.y - two.y) < TOL) && (fabs(one.z - two.z) < TOL))
	vector_t * planes = NULL;
	int count=0,i=0,size=0;
	
	/* get the planes for the Voronoid polyhedron */
	for (i=0;i<crystal.atomc;i++)
	{
		int j;
		vector_t new;
		if (EQ(crystal.atoms[i],atom))
			continue;
		new.x = (crystal.atoms[i].x - atom.x)/2;
		new.y = (crystal.atoms[i].y - atom.y)/2;
		new.z = (crystal.atoms[i].z - atom.z)/2;
		for (j=0;j<count;j++)
		{
			double dot = DOT(planes[j],new);
			if (dot > SQMAG(planes[j]))
			{
				new.x = new.y = new.z = NAN;
				break;
			}
			else if (dot > SQMAG(new))
			{
				int k;
				count--;
				for (k=j;k<count;k++)
					planes[k] = planes[k+1];
			}
		}
		if (! (isnan(new.x) || isnan(new.y) || isnan(new.z)))
		{
			if (count >= size)
			{
				size++;
				planes = realloc(planes,size*sizeof(vector_t));
				if (planes == NULL) fail("realloc failed");
			}
			planes[count++] = new;
		}
	}
	*c = count;
	return planes;
#undef EQ
}

/** Get the vertices of the voronoi polyhedron around the atom */
vector_t * vertices(crystal_t crystal,atom_t atom,int * c)
{
	int count=0,vcount=0,vsize,i;
	vector_t * planes = voronoi(crystal,atom,&count);
	vector_t * verts = malloc(count*sizeof(vector_t));
	vsize = count;

	/* find the farthest vertex for the distance of Rcc */
	for (i=0;i<count;i++)
	{
		int j;
		double sqmag_i = SQMAG(planes[i]);
		/* normalised vector */
		for (j=i+1;j<count;j++)
		{
			int k;
			double sqmag_j = SQMAG(planes[j]);
			/* Normalised vector */
			vector_t i_x_j = cross(planes[i],planes[j]);
			for (k=j+1;k<count;k++)
			{
				vector_t vertex, j_x_k, k_x_i;
				double sqmag_k = SQMAG(planes[k]),div;
				/* normalised vector */
				int check;
				
				j_x_k = cross(planes[j],planes[k]);
				k_x_i = cross(planes[k],planes[i]);

				div = DOT(planes[i],j_x_k);
				
				vertex.x = sqmag_i*j_x_k.x + sqmag_j*k_x_i.x + sqmag_k*i_x_j.x;
				vertex.y = sqmag_i*j_x_k.y + sqmag_j*k_x_i.y + sqmag_k*i_x_j.y;
				vertex.z = sqmag_i*j_x_k.z + sqmag_j*k_x_i.z + sqmag_k*i_x_j.z;

				vertex.x /= div;
				vertex.y /= div;
				vertex.z /= div;

				for (check = 0;check<count;check++)
				{
					if ((check == i) || (check == j) || (check == k))
						continue;

					if (DOT(planes[check],vertex) > SQMAG(planes[check]))
					{
						vertex.x = vertex.y = vertex.z = NAN;
						break;
					}
				}
				if (isnan(vertex.x) || isnan(vertex.y) || isnan(vertex.z))
					continue;
				
				if (vcount >= vsize)
				{
					vsize = vcount + 1;
					verts = realloc(verts,vsize*sizeof(vector_t));
					if (verts == NULL) fail("realloc failed");
				}
				verts[vcount++] = vertex;
			}
		}
	}

	/* free the planes */
	free(planes);
	
	*c = vcount;
	return verts;
}

/** Get the radius of the circumscribing sphere of the Voronoi polyhedron of 
  * the atom */
double Rcc(crystal_t crystal,atom_t atom)
{
	int count=0,i;
	vector_t * verts = vertices(crystal,atom,&count);
	double r,mag;
	
	/* find the farthest vertex for the distance of Rcc */
	r=0;
	for (i=0;i<count;i++)
	{
		mag = MAGNITUDE(verts[i]);
		if (mag > r)
			r = mag;
	}

	/* free the planes */
	free(verts);
	return r;
}

/** Get the radius of the inscribing sphere of the Voronoi polyhedron of the
  * atom */
double Ric(crystal_t crystal,atom_t atom)
{
	int i,count=0;
	vector_t * planes = voronoi(crystal,atom,&count);
	double r = INFINITY;

	for (i=0;i<count;i++)
	{
		double mag = MAGNITUDE(planes[i]);
		if (mag < r)
			r = mag;
	}

	free(planes);
	return r;
}

/** The ratio of Rcc to closest_c */
double r_c_ratio(crystal_t crystal,atom_t atom)
{
	return (closest_c(crystal,atom)/Rcc(crystal,atom));
}

/** cross product */
inline vector_t cross(vector_t one,vector_t two)
{
	vector_t out;
	out.x = one.y*two.z - two.y*one.z;
	out.y = one.z*two.x - two.z*one.x;
	out.z = one.x*two.y - two.x*one.y;
	return out;
}

/** gets the unit vector gradient of the function func at the point of atom.
  */
vector_t grad(atom_t atom,crystal_t crystal,double (*func)(crystal_t, atom_t)
	,double epsilon)
{
	double val;
	vector_t grad;	

	/* get the partial derivative with respect to x */
	atom.x += epsilon;
	grad.x = func(crystal,atom);
	atom.x -= 2*epsilon;
	grad.x -= func(crystal,atom);
	atom.x += epsilon;

	/* get the partial derivative with respect to y */
	atom.y += epsilon;
	grad.y = func(crystal,atom);
	atom.y -= 2*epsilon;
	grad.y -= func(crystal,atom);
	atom.y += epsilon;

	/* get the partial derivative with respect to z */
	atom.z += epsilon;
	grad.z = func(crystal,atom);
	atom.z -= 2*epsilon;
	grad.z -= func(crystal,atom);
	atom.z += epsilon;

	/* get the magnitude of grad */
	val = MAGNITUDE(grad);
	
	if (fabs(val) < TOL)
	{
		grad.x = 0.0;
		grad.y = 0.0;
		grad.z = 0.0;
	}
	else
	{
		/* normalise it */
		grad.x /= val;
		grad.y /= val;
		grad.z /= val;
	}	
	return grad;
}

/** Moves the atom up the grad vector field to a maximum point on the function.
  */
atom_t drift(atom_t atom, crystal_t crystal, double (*func)(crystal_t, atom_t),
	double d,double dmin)
{
	atom_t tmp;
	double near,far,epsilon;
	int c = 0,limit = 1E4;

	d *= radius(atom);

	tmp=atom;
	far = func(crystal,atom);
	near = far - 1;

	while (d >= dmin)
	{
		dbg_printf("Step size: %lf\n",d);
		c=0;
		epsilon = 0.1*d;
		while ((far > near) && (c < limit))
		{
			dbg_printf("Step %d\n",c);
			vector_t direction;
			atom = tmp;	
			near = far;
			
			if (get_coefficients(&atom))
				bringin(&atom,crystal);
		
			direction = grad(atom,crystal,func,epsilon);
			
			if ((MAGNITUDE(direction)) < 0.9)
				return atom;
			
			tmp = move(atom,direction,d);
			tmp.Z = atom.Z;
			tmp.rad = atom.rad;
	
			far = func(crystal,tmp);
			c++;
		}
		near = -d*radius(atom);
		d /= 10;
	}

	if (get_coefficients(&atom))
			bringin(&atom,crystal);
	
	return atom;
}

/** Make MPI types for passing atom_t, vector_t and crystal_t structs. 
  * places them in newtypes in that order. */
inline void make_mpi_types(MPI_Datatype * newtypes)
{
#define MPI_ATOM_T		(newtypes[0])
#define MPI_VECTOR_T	(newtypes[1])
#define MPI_CRYSTAL_T	(newtypes[2])
	MPI_Datatype oldtypes[3];
	MPI_Aint offsets[3], lb[3];
	int lengths[3];

	/* Make the oldtypes and lengths arrays for atom_t */
	oldtypes[0] = MPI_DOUBLE;
	oldtypes[1] = MPI_INT;
	lengths[0] = 7;
	lengths[1] = 1;
	offsets[0] = 0;
	
	/* Get displacement of int */
	MPI_Type_get_extent(oldtypes[0],&(lb[1]),&(offsets[1]));
	offsets[1] *= lengths[0];
	offsets[1] += offsets[0];

	/* create the MPI_ATOM_T type for moving atoms around. */
	MPI_Type_create_struct(2,lengths,offsets,oldtypes,&MPI_ATOM_T);
	MPI_Type_commit(&MPI_ATOM_T);
	
	/* Fix lengths for creating the vector type */
	MPI_Type_contiguous(3,MPI_DOUBLE,&MPI_VECTOR_T);
	MPI_Type_commit(&MPI_VECTOR_T);

	/* Change variables to fit crystal_t */
	oldtypes[0] = MPI_VECTOR_T;
	oldtypes[1] = MPI_UNSIGNED_LONG;
	oldtypes[2] = MPI_INT;

	lengths[0] = 3;
	lengths[1] = 1;
	lengths[2] = 2;

	offsets[0] = 0;
	
	MPI_Type_get_extent(oldtypes[0],&(lb[1]),&(offsets[1]));
	offsets[1] *= lengths[0];
	offsets[1] += offsets[0];

	MPI_Type_get_extent(oldtypes[1],&(lb[2]),&(offsets[2]));
	offsets[2] *= lengths[1];
	offsets[2] += offsets[1];
	
	MPI_Type_create_struct(3,lengths,offsets,oldtypes,
		&MPI_CRYSTAL_T);
	MPI_Type_commit(&MPI_CRYSTAL_T);

#undef MPI_ATOM_T	
#undef MPI_VECTOR_T
#undef MPI_CRYSTAL_T
}

/** Rounded down log base 2. */
int floor_log2(int n)
{
	int log = 0, i;
	if (n <= 0)
		return -1;
	
	/* Start at i = 16 and go down to 1 by halfing each time */
	for (i=16;i<0;i >>= 1)
	{
		if ( n >> i)
		{
			n >>= i;
			log += i;
		}
	}
	return log;
}
