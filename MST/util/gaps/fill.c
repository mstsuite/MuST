/** fill.c 
 * 
 * contains all of the filling algorithms. 
 * 
 * George Bargoud: <gbargoud@gmail.com> 
 **/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "gapslib.h"
#include "fill.h"

#define COUNT	3
#define ATOMS	5

/** Builds up a grid with cells that can each fit a virtual atom and then tries
  * to fit one on each vertex after drifting it. Multithreaded. */
void gridfill(crystal_t * crystal)
{
#define MPI_ATOM_T		(types[0])
#define MPI_VECTOR_T	(types[1])
#define MPI_CRYSTAL_T	(types[2])
	MPI_Datatype types[3];
	/* full is the crystal extended on all sides, tmp stores new atoms found
	 * in this process for later */
	crystal_t * full;
	/* the rank of the process in MPI_COMM_WORLD and the size of MPI_COMM_WORLD
	 */
	int rank,size,i;
	
	/* Used for dividing up the crystal among the processes and then iterating 
	 * over it. */
	double ahi,alo,da,a;
	double bhi,blo,db,b;
	double chi,clo,dc,c;

	make_mpi_types(types);

	/* Broadcast the crystal to all threads */
	if (crystal)
	{
		MPI_Bcast(crystal,1,MPI_CRYSTAL_T,0,MPI_COMM_WORLD);
		MPI_Bcast(crystal->atoms,crystal->atomc,MPI_ATOM_T,0,MPI_COMM_WORLD);
	}
	else
	{
		crystal = malloc(sizeof(crystal_t));
		MPI_Bcast(crystal,1,MPI_CRYSTAL_T,0,MPI_COMM_WORLD);
		crystal->atoms = calloc(sizeof(atom_t),crystal->atomsiz);
		MPI_Bcast(crystal->atoms,crystal->atomc,MPI_ATOM_T,0,MPI_COMM_WORLD);
	}
	
	/* get the rank and total count */
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	/* set up full */
	full = extend(*crystal);
	dbg_printf("%d: full has %d atoms.\n",rank,full->atomc);	
	/* If rank is not 0, remove all atoms from crystal so as not to send 
	 * duplicates at the end. */
	if (rank)
	{
		crystal->atomc = 0;
		crystal->atomsiz = 1;
		free(crystal->atoms);
		crystal->atoms = calloc(crystal->atomsiz,sizeof(atom_t));
	}
	
	/* The step distance in each direction */
	da = 1.5*radii[type]/MAGNITUDE(crystal->a);
	db = 1.5*radii[type]/MAGNITUDE(crystal->b);
	dc = 1.5*radii[type]/MAGNITUDE(crystal->c);
	
	/* calculate the start and end point on the crystal for each process. */
	/* TODO: make a better system for splitting */
	
	alo   = ((double)(rank/2))/((double)((size+1)/2)) - 0.1*da;
	ahi   = ((double)((rank/2)+1))/((double)((size+1)/2)) + 0.1*da;
	
	blo   = ((rank%2)?0.5:0.0) - 0.1*db;
	bhi   = ((rank%2)?1.0:((rank+1 < size) ? 0.5:1.0)) + 0.1*db;
	
	clo   = 0.0 - 0.1*dc;
	chi   = 1.0 + 0.1*dc;
	
	dbg_printf("%d: assigned:\na\t%lf->%lf\nb\t%lf->%lf\nc\t%lf->%lf\n",
		rank,alo,ahi,blo,bhi,clo,chi);
	/* Loop over the search until it finds nothing. */
	do
	{
		i=0;
		/* Find the holes in between alo, blo, clo and ahi, bhi and chi */
		for (a=alo;a<=ahi;a+=da)
		{
			for (b=blo;b<=bhi;b+=db)
			{
				for (c=clo;c<=chi;c+=dc)
				{
					atom_t tmp;
					double close;
	
					tmp.x = (a*full->a.x)+(b*full->b.x)+(c*full->c.x);
					tmp.y = (a*full->a.y)+(b*full->b.y)+(c*full->c.y);
					tmp.z = (a*full->a.z)+(b*full->b.z)+(c*full->c.z);
					tmp.rad = radii[type];
					tmp.Z = 0;
	
					/* drift starting at intervals of 0.01 and shrinking to
					 * intervals of 1E-9. A smaller second value (dmin)
					 * means more accurate results. */
					dbg_printf("%d: Trying  %lf,%lf,%lf\n",
						rank,tmp.x,tmp.y,tmp.z);
					tmp=drift(tmp,*full,closest,0.1,1E-15);
					dbg_printf("%d: Reached %lf,%lf,%lf\n",
						rank,tmp.x,tmp.y,tmp.z);
					
					if (get_coefficients(&tmp))
					{
						dbg_printf("%d: Bringing it in.\n",rank);
						bringin(&tmp,*full);
					}
					close = closest(*full,tmp);
					if (close > -0.05*radius(tmp))
					{
						double ratio;
						dbg_printf("%d: passed proximity test at %lf,%lf,%lf\n"
							,rank,tmp.x,tmp.y,tmp.z);
						tmp.rad += close;
						ratio = Rcc(*full,tmp)/closest_c(*full,tmp);
						dbg_printf("%d: ratio is %lf\n",rank,ratio);
						if (ratio <= 1.00001)
						{
							dbg_printf("\t%d: Found an atom at %lf,%lf,%lf\n"
								,rank,tmp.x,tmp.y,tmp.z);
							append(crystal,tmp);
							appendall(full,tmp);
							dbg_printf("\t%d: Now has %d atoms.\n",rank,
									crystal->atomc);
							i++;
						}
					}	
				}
			}
		}
	} while (i);

	/* Put it all together */
	if (rank)
	{	
		/* send the count */
		MPI_Send(&(crystal->atomc),1,MPI_INT,0,COUNT,MPI_COMM_WORLD);
		/* send the atoms themselves */
		dbg_printf("%d: Sending %d atoms\n",rank,crystal->atomc);
		MPI_Send(crystal->atoms,crystal->atomc,MPI_ATOM_T,0,ATOMS,
			MPI_COMM_WORLD);
		dbg_printf("%d: Sent %d atoms\n",rank,crystal->atomc);
		free(crystal->atoms);
	}
	else
	{
		for (i=1;i<size;i++)
		{
			int count,j,source;
			MPI_Status status;
			atom_t * atoms;
			
			dbg_printf("Waiting for atom count\n");
			/* Receive a count of atoms */
			MPI_Recv(&count,1,MPI_INT,MPI_ANY_SOURCE,COUNT,MPI_COMM_WORLD,
				&status);
			source = status.MPI_SOURCE;
			/* Allocate the space */
			atoms = calloc(sizeof(atom_t),count);
			if (!atoms)
				fail("calloc failed");
			dbg_printf("Waiting for %d atoms from %d\n",count,source);
			/* Receive the atoms. */
			MPI_Recv(atoms,count,MPI_ATOM_T,source,ATOMS,MPI_COMM_WORLD,
				&status);
			dbg_printf("Received %d atoms from %d\n",count,source);	
			/* Check if they fit */
			for (j=0;j<count;j++)
			{
				double close;
				dbg_printf("Checking atom %d:%lf,%lf,%lf\n",atoms[j].Z,
					atoms[j].x,atoms[j].y,atoms[j].z);
				close = closest(*full,atoms[j]);
				if (close > -0.05*radius(atoms[j]))
				{
					atoms[j].rad += close;
					// atoms[j]=drift(atoms[j],*full,r_c_ratio,0.1,1E-4);
					if (Rcc(*full,atoms[j]) < 0.9*closest_c(*full,atoms[j]))
					{
						append(crystal,atoms[j]);
						appendall(full,atoms[j]);
					}
				}
			}
			/* Free the memory */
			free(atoms);
		}
	}
#undef MPI_ATOM_T
#undef MPI_VECTOR_T
#undef MPI_CRYSTAL_T
}
