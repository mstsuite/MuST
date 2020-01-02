/** gaps.c 
 * 
 * George Bargoud: <gbargoud@gmail.com>
 **/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <regex.h>
#include "mpi.h"

#include "gaps.h"

#define TOL 1.0e-15

int main (int argc, char ** argv)
{
	char c;
	crystal_t crystal;
	redef_t * redefs = NULL;
	int count = 0;
	/* MPI error code and rank */
	int rc,rank;

	/* Select default filling algorithm. */	
	fill = funcmap[0].func;

	/* Read the options and parse them */
	while ((c = getopt(argc,argv,"gwhi:o:r:u:a:c:t:d:")) != -1)
	{
		switch (c)
		{
			/* print help text */
			case 'h':
				help();
				break;

			/* define input file */
			case 'i':
				inFileName = optarg;
				break;

			/* define output file */
			case 'o':
				outFileName = optarg;
				break;
			
			/* What has direct coordinates */
			case 'd':
				switch (*optarg)
				{
					case 'i':
						direct = 1;
						break;
					case 'o':
						direct_out = 1;
						break;
					case 'n':
						direct = direct_out = 0;
						break;
					case 'a':
						direct = direct_out = 1;
						break;
					default:
						fail("%c not i o n or a.",*optarg);
				}
				break;

			/* redefine a radius, in format "Z=A" where Z is the 
			 * atomic number and A is the new radius in the last 
			 * selected units.*/
			case 'r':
				redefine(optarg);
				break;

			/* Choose the filling algorithm. */
			case 'a':
				fill = funcfromname(optarg);
				break;

			/* Choose the units */
			case 'u':
				units(optarg);
				break;

			/* Read a configuraton file. */
			case 'c':
				parse_conf(optarg);
				break;

			/* Guess radii. */
			case 'g':
				g = 1;
				break;
			
			/* Choose atom type for virtual particles */
			case 't':
				type = lookup(optarg);
				break;
			
			/* Activate warnings */
			case 'w':
				warnings = 1;
				break;

			/* Unknown option: print help text. */
			case '?':
				help();
				break;
		}
	}
	
	/* Start MPI */
	rc = MPI_Init(&argc,&argv);
	CHECK(rc);
	rc = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	CHECK(rc);

	if (optind < argc)
	{
		inFileName = argv[optind];
	}

	if (rank == 0)
	{
		crystal = parse(inFileName);
	}
	
	/* Broadcast the transformation matrix */
	MPI_Bcast(transform,9,MPI_DOUBLE,0,MPI_COMM_WORLD);

	/* If "-g" was in the options, do all the guessing stuff here. */
	if (g)
	{
		if (rank == 0)
			redefs = guess(crystal,&count);
		
		MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
		if (count)
		{
			/* Create an MPI_Datatype to send the redefs array */
			MPI_Aint offsets[2],lb;
			MPI_Datatype oldtypes[2],MPI_REDEF;
			int lengths[2],i;

			oldtypes[0] = MPI_DOUBLE;
			oldtypes[1] = MPI_INT;

			lengths[0] = 1;
			lengths[1] = 1;

			offsets[0] = 0;
			MPI_Type_get_extent(oldtypes[0],&lb,&(offsets[1]));
			
			MPI_Type_create_struct(2,lengths,offsets,oldtypes,&MPI_REDEF);
			MPI_Type_commit(&MPI_REDEF);
			
			/* In the threads of nonzero rank, create a new redefs array */
			if (rank > 0)
			{
				redefs = calloc(sizeof(redef_t),count);
				if (redefs == NULL)
					fail("calloc error");
			}
			MPI_Bcast(redefs,count,MPI_REDEF,0,MPI_COMM_WORLD);
			
			/* redefine the radii */
			for (i=0;i<count;i++)
			{
				if (redefs[i].Z < 0)
				{
					radii[-(redefs[i].Z)] = redefs[i].A;
				}
				else
				{
					radii[redefs[i].Z] = redefs[i].A;
				}
			}
		}
		/* Broadcast the type of virtual atom in case it changed */
		MPI_Bcast(&type,1,MPI_INT,0,MPI_COMM_WORLD);
	}

	if (type == -1)
		type = 0;
	
	if (rank == 0)
	{
		fill(&crystal);
		/* Output the full crystal */
		out(outFileName,crystal,redefs,count);
	}
	else
	{
		free(redefs);
		fill(NULL);
	}

	MPI_Finalize();
	exit(EXIT_SUCCESS);
}

/** takes in the string "Z=A" where Z is an atomic number and A is a radius
  * for that atom in Angstroms and sets A as the radius for Z in the radii 
  * array */
void redefine(char* s)
{
	int Z;
	double A;
	char * sym;

	/* parse the string and make sure the format is correct */
	if (sscanf(s,"%d=%lf",&Z,&A) != 2)
	{
		sym = strtok(s,"= \t");
		A = atof(strtok(NULL,"= \t\n"));

		Z = lookup(sym);
		if (Z < 0)
		{
			fail("%s is not a valid atomic symbol",sym);
		}
	}

	radii[Z] = A;

}

/** Guess radii for atoms and return an array of those guesses */
redef_t * guess(crystal_t crystal, int * count)
{
	redef_t * redefs = NULL;
	int i,tmptype=0;
	
	dbg_printf("Guessing radii\n");
	*count = 0;

	for (i=0;i<crystal.atomc;i++)
	{
		int j;
		double closest = 1.0/0.0;
		int closest_Z = 0;
		/* If type had not been defined, make it the smallest atom in the
		 * crystal */
		if (type < 0)
		{
			if ((!tmptype) || (radius(crystal.atoms[i]) < radii[tmptype]))
				tmptype = crystal.atoms[i].Z;
		}
		for (j=(i+1);j<crystal.atomc;j++)
		{
			double distance = dist(crystal.atoms[i],crystal.atoms[j]);
			distance -= radius(crystal.atoms[i])+radius(crystal.atoms[j]);
			if (distance < closest)
			{
				closest = distance;
				closest_Z = crystal.atoms[j].Z;
			}
		}
		if (closest < 0)
		{
			double factor;
			redef_t this, that;
			atom_t tmp;

			this.Z = crystal.atoms[i].Z;
			that.Z = closest_Z;
			
			tmp.Z = that.Z;

			this.A = radius(crystal.atoms[i]);
			that.A = radius(tmp);

			/* The factor each should be changed by */
			factor = 1 + (closest/(this.A + that.A));

			this.A *= 0.999*factor;
			that.A *= 0.999*factor;

			/* Check if it is already there. */
			for (j=0;j<(*count);j++)
			{
				/* If they refer to the same particle, real or virtual */
				if ((redefs[j].Z == this.Z) || (redefs[j].Z == -(this.Z)))
				{
					/* Only replace it if it is smaller */
					if (this.A < redefs[j].A)
					{
						/* Change the radius */
						if (this.Z < 0)
						{
							this.Z = -(this.Z);
						}
						redefs[j] = this;
						radii[this.Z] = this.A;
						dbg_printf("redefining %d to %lf\n",this.Z,this.A);
					}
					this.A = NAN;
					break;
				}
			}

			/* If they were not there already, add them in */
			/* if (this.A == this.A) */
			if (!isnan(this.A))
			{
				(*count)++;
				redefs = realloc(redefs,sizeof(redef_t)*(*count));
				if (redefs == NULL)
					fail("realloc failed");
				
				if (this.Z < 0)
				{
					this.Z = -(this.Z);
				}
				/* Append the new value to redefs */
				redefs[(*count)-1] = this;
				radii[this.Z] = this.A;
				/* Change the radius */
				dbg_printf("redefining %d to %lf\n",this.Z,this.A);
			}
			
			/* Do the same for the other redef. */
			for (j=0;j<(*count);j++)
			{
				if ((redefs[j].Z == that.Z) || (redefs[j].Z == -(that.Z)))
				{
					/* Only replace it if it is smaller */
					if (that.A < redefs[j].A)
					{
						if (that.Z < 0)
						{
							that.Z = -(that.Z);
						}
						redefs[j] = that;
						radii[that.Z] = that.A;
						dbg_printf("redefining %d to %lf\n",that.Z,that.A);
					}
					that.A = NAN;
					break;
				}
			}
			
			/* if (that.A == that.A) */
			if (!isnan(that.A))
			{
				/* (*count) cannot be 0 here. */
				(*count)++;
				redefs = realloc(redefs,sizeof(redef_t)*(*count));
				if (redefs == NULL)
					fail("realloc failed");
				
				if (that.Z < 0)
				{
					that.Z = -(that.Z);
				}
				redefs[(*count)-1] = that;
				radii[that.Z] = that.A;
				dbg_printf("redefining %d to %lf\n",that.Z,that.A);
			}
		}
	}
	/* If type was not defined, set it to tmptype */
	if (type < 0)
	{
		dbg_printf("Using %d for virtual type\n",tmptype);
		type = tmptype;
	}
	return redefs;
}

/** Parse the crystal from a file. While parsing the file it also fixes
 * any atoms that are not within the paralelipiped that is described. */
crystal_t parse(char* filename)
{
	FILE * fp;
	atom_t atom;
	crystal_t crys;
	char line[BUFSIZ];
	regex_t vector_re, atom_re;
	regmatch_t matches[6];
	vector_t a, b, c;
	
	/* This is to stop some compilers from sending warning flags. */
	a.x = 0.0; a.y = 0.0; a.z = 0.0;
	b.x = 0.0; b.y = 0.0; b.z = 0.0;
	c.x = 0.0; c.y = 0.0; c.z = 0.0;

	if (filename == NULL)
	{
		fp = stdin;
	}
	else
	{
		fp = fopen(filename,"r");
	}

	if (fp == NULL)
	{
		fail("file error while opening %s",filename);
	}

	crys.atomc = 0;
	crys.atomsiz = 2;

	/* Read the vectors for the dimensions
	 * (-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*) matches a floating point number
	 * in C or FORTRAN scientific format or regular floating point format. */
	regcomp(&vector_re,"^[ \t]*"
		"(-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*)[ \t]*[,; \t][ \t]*"
		"(-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*)[ \t]*[,; \t][ \t]*"
		"(-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*)",
		REG_EXTENDED);
	while (fgets(line,sizeof(line),fp))
	{
		if (!regexec(&vector_re,line,4,matches,0))
		{
			ftoc(line+matches[1].rm_so);
			a.x = atof(line+matches[1].rm_so);
			a.y = atof(line+matches[2].rm_so);
			a.z = atof(line+matches[3].rm_so);
			crys.a = a;
			break;
		}
	}
	dbg_printf("Read in a: %lf,%lf,%lf\n",crys.a.x,crys.a.y,crys.a.z);
	while (fgets(line,sizeof(line),fp))
	{
		if (!regexec(&vector_re,line,4,matches,0))
		{
			ftoc(line+matches[1].rm_so);
			b.x = atof(line+matches[1].rm_so);
			b.y = atof(line+matches[2].rm_so);
			b.z = atof(line+matches[3].rm_so);
			crys.b = b;
			break;
		}
	}
	dbg_printf("Read in b: %lf,%lf,%lf\n",crys.b.x,crys.b.y,crys.b.z);
	while (fgets(line,sizeof(line),fp))
	{
		if (!regexec(&vector_re,line,4,matches,0))
		{
			ftoc(line+matches[1].rm_so);
			c.x = atof(line+matches[1].rm_so);
			c.y = atof(line+matches[2].rm_so);
			c.z = atof(line+matches[3].rm_so);
			crys.c = c;
			break;
		}
	}
	dbg_printf("Read in c: %lf,%lf,%lf\n",crys.c.x,crys.c.y,crys.c.z);
	
	regfree(&vector_re);

	/* Set the transformation vector in the crystal. This may switch the 
	 * bounding vectors around if needed for inverting the matrix they form. */
	set_transform(&crys);

	/* Read in the atoms */
	crys.atoms = calloc(crys.atomsiz,sizeof(atom_t));
	if (crys.atoms == NULL)
	{
		fail("calloc error");
	}
	
	/* Atom regular expression. Matches:
	 *	index, atomic number, x, y, z
	 *	index, chemical symbol, x, y, z
	 *	atomic number, x, y, z
	 *	chemical symbol, x, y, z
	 * with delimiters being either whitespace commas or semicolons. 
	 * Ignores index.*/
	regcomp(&atom_re,"^[ \t]*([0-9]*[ \t]*[;, \t])?[ \t]*"
		"(-?[0-9A-Z][0-9a-z]*[.]?[0-9]*[eED]?[+-]?[0-9]*)[ \t]*[,; \t][ \t]*"
		"(-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*)[ \t]*[,; \t][ \t]*"
		"(-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*)[ \t]*[,; \t][ \t]*"
		"(-?[0-9]*[.]?[0-9]+[eED]?[+-]?[0-9]*)",
		REG_EXTENDED);
	while (fgets(line,sizeof(line),fp))
	{
		int d;
		/* char symbol[BUFSIZ]; */
		char symbol[3];
		
		/* If the line does not match atom coordinates */
		if (regexec(&atom_re,line,6,matches,0))
		{
			/* check for the word DirectCoordinates and make sure it does not
			 * come affter one of the comment marks. */
			char * dc, * comment;
			dc = strstr(line,"DirectCoordinates");
			if (dc != NULL)
			{
				comment = strpbrk(line,"!#");
				if ((comment == NULL) || (comment > dc))
				{
					dbg_printf("Using direct coordinates.\n");
					direct = 1;
				}
			}
			continue;
		}

		/* extract the chemical symbol or atomic number from line */
		strncpy(symbol,line+matches[2].rm_so,
                        2);
                symbol[2] = '\0';
	             /* matches[2].rm_eo-matches[2].rm_so); */

		if (strpbrk(symbol,"."))
		{
			ftoc(symbol);
			atom.Z = 0;
			atom.rad = atof(symbol);
			dbg_printf("Read in radius of %lf\n",atom.rad);
		}
		else
		{
			atom.Z = lookup(symbol);
			if (atom.Z == 0)
				atom.rad = radii[0];
		}

		ftoc(line+matches[3].rm_so);
		if (direct)
		{
			/* matches[3].rm_so is the starting offset for Ca */
			atom.Ca = atof(line+matches[3].rm_so);
			/* matches[4].rm_so is the starting offset for Cb */
			atom.Cb = atof(line+matches[4].rm_so);
			/* matches[5].rm_so is the starting offset for Cc */
			atom.Cc = atof(line+matches[5].rm_so);
			
			atom.x = atom.Ca*a.x + atom.Cb*b.x + atom.Cc*c.x;
			atom.y = atom.Ca*a.y + atom.Cb*b.y + atom.Cc*c.y;
			atom.z = atom.Ca*a.z + atom.Cb*b.z + atom.Cc*c.z;
		}
		else
		{
			/* matches[3].rm_so is the starting offset for x */
			atom.x = atof(line+matches[3].rm_so);
			/* matches[4].rm_so is the starting offset for y */
			atom.y = atof(line+matches[4].rm_so);
			/* matches[5].rm_so is the starting offset for z */
			atom.z = atof(line+matches[5].rm_so);
		}
		
		if (atom.Z == 0)
		{
			dbg_printf("radius: %lf\n",atom.rad);
		}
		else
		{
			dbg_printf("type %d\n",atom.Z);
		}
		dbg_printf("\t Position: %lf,%lf,%lf\n",atom.x,atom.y,atom.z);

		/* In case a, b and c were switched, the coefficients need to be 
		 * recalculated. */
		d = get_coefficients(&atom);

		dbg_printf("\tdirect coordinates: %lf,%lf,%lf\n",
			atom.Ca,atom.Cb,atom.Cc);

		if (d)
		{
			warn_printf("WARNING: atom %d of input file corrected.\n",
				crys.atomc);
			warn_printf("\tOld values: %lf, %lf, %lf\n",atom.x,atom.y,atom.z);
			
			bringin(&atom,crys);
			
			warn_printf("\tNew values: %lf, %lf, %lf\n",atom.x,atom.y,atom.z);
			fflush(stderr);
		}
		/* Append the new atom to the array */
		append(&crys,atom);
	}
	
	regfree(&atom_re);

	if (fp != stdin)
		fclose(fp);
	
	return crys;
}

/** Changes floating point numbers from fortran to c syntax */
static inline void ftoc(char * s)
{
	while(*s)
	{
		if (*s == 'D')
			*s = 'E';
		s++;
	}
}

/** Looks up the atomic number based on the symbol. */
int lookup(const char * symbol)
{
	int Z=0,sign=1;

	/* Check if it is virtual */
	if (*symbol == '-')
	{
		sign = -1;
		symbol++;
	}

	/* If symbol can be made into a nonzero integer or starts with a zero */
	if (atoi(symbol) || *symbol == '0')
		return sign*atoi(symbol);
	
	/* NULL is the end of the table */
	while(symbols[Z] != NULL)
	{
		if (!strcmp(symbol,symbols[Z]))
		{
			return sign*Z;
		}
		Z++;
	}
	fail("%s is not a valid atomic number or symbol.",symbol-(sign<0?1:0));
}
/** Gets the algorithm chosen based on its name or fails with an error if it
  * doesn't exist. returns a function pointer. */
fill_t funcfromname(char * name)
{
	unsigned int i;
	
	/* While the name is not NULL, go through the list */
	for (i=0u;funcmap[i].name;i++)
	{
		if (!strcmp(funcmap[i].name,name))
		{
			return funcmap[i].func;
		}
	}

	/* If it was not found, fail */
	fail("%s is not a valid algorithm name.",name);
}

/** parse configuration file for options in it */
void parse_conf(char * filename)
{
	/* line from the file */
	char line[BUFSIZ];
	FILE * fp;
	
	fp = fopen(filename,"r");
	if (fp == NULL)
		fail("invalid filename: %s",filename);
	while (fgets(line,BUFSIZ,fp))
	{
			char *set,*val;
			set = strtok(line," =\t\n");
			val = strtok(NULL," =\t\n");
			
			if ((set == NULL) || (val == NULL))
				continue;
			
			switch (*set)
			{
				/* set algorithm name */
				case 'a':
				case 'A':
					fill = funcfromname(val);
					break;
				/* set input file name */
				case 'i':
				case 'I':
					inFileName = stralloc(val);
					if (inFileName == NULL)
						fail("calloc error");
					strcpy(inFileName,val);
					break;
				/* set output file name */
				case 'o':
				case 'O':
					outFileName = stralloc(val);
					if (outFileName == NULL)
						fail("calloc error");
					strcpy(outFileName,val);
					break;
				/* set units */
				case 'u':
				case 'U':
					units(val);
					break;
				/* resize atom, other method */
				case 'r':
				case 'R':
					redefine(val);
					break;
				/* define the type for the atom */
				case 't':
				case 'T':
					type = lookup(val);
					break;
				/* use direct coordinates for input, output, all or nothing. */
				case 'd':
				case 'D':
					switch (*val)
					{
						case 'i':
						case 'I':
							direct = 1;
							break;
						case 'o':
						case 'O':
							direct_out = 1;
							break;
						case 'n':
						case 'N':
							direct = direct_out = 0;
							break;
						case 'a':
						case 'A':
							direct = direct_out = 1;
							break;
					}
					break;
				/* guess is on */
				case 'g':
				case 'G':
					g=1;
					break;
				/* resize atom */
				case '0':
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				case '7':
				case '8':
				case '9':
					sprintf(line,"%s=%s",set,val);
					redefine(line);
					break;
				default:
					continue;
			}
	}
	fclose(fp);
}

/** Output the crystal to the described file. If the filename is NULL, send it
  * to standard out. */
void out(char * filename,crystal_t crystal,redef_t * redefs, int count)
{
	int i;
	/* if no output file selected, output to stdout */
	FILE * fp = (filename == NULL? stdout:fopen(filename,"w"));
	char * format = (direct_out?"%.16lE, %.16lE, %.16lE\n":
		"%.16lf, %.16lf, %.16lf\n");

	if (fp == NULL)
	{
		fail("error in opening %s",filename);
	}

	fprintf(fp,"# scale=%lE Angstroms\n\n",unit);

	/* Print the vectors */
	
	fprintf(fp,format,crystal.a.x,crystal.a.y,crystal.a.z);
	fprintf(fp,format,crystal.b.x,crystal.b.y,crystal.b.z);
	fprintf(fp,format,crystal.c.x,crystal.c.y,crystal.c.z);
	fprintf(fp,"\n");
	
	if (count)
	{
		fprintf(fp,"# Redefined atomic radii:\n");
		for (i=0;i<count;i++)
			fprintf(fp,"# %s=%lf\n",symbols[redefs[i].Z],redefs[i].A);
		fprintf(fp,"\n");
	}

	if (direct_out)
		fprintf(fp,"DirectCoordinates\n");

	/* Print atoms, including the virtual ones */
	for (i=0; i<crystal.atomc; i++)
	{
		atom_t atom = crystal.atoms[i];

		if (atom.Z >= 0)
			fprintf(fp,"%s, ",symbols[atom.Z]);
		/* else if (atom.Z == 0)
			fprintf(fp,"%lf, ",atom.rad); */
		else
			fprintf(fp,"%d, ",atom.Z);


		/* If DirectCoordinates are enabled, print them. */
		if (direct_out)
			fprintf(fp,format,atom.Ca,atom.Cb,atom.Cc);
		else
			fprintf(fp,format,atom.x,atom.y,atom.z);
	}

	if (fp != stdout)
	{
		fclose(fp);
	}

	return;
}

/** Set the transform matrix from cartesian coordinates to a basis made up of
  * the edges of the crystal */
void set_transform(crystal_t * crys)
{
	/* Values that will be reused */
	double invert[3][3];
	
	/* indices for the transformation later */
	int i,j;

	/* Set transform to the identity for Gauss-Jordan elimination */
	transform[0][0] = 1.0;
	transform[1][1] = 1.0;
	transform[2][2] = 1.0;

	/* Change the vectors around to make an invertible matrix. */
	/* check the first vector */
	if (fabs(crys->a.x) < TOL)
	{
		vector_t tmp = crys->a;
		if (fabs(crys->b.x) < TOL)
		{
			crys->a = crys->c;
			crys->c = tmp;
		}
		else
		{
			crys->a = crys->b;
			crys->b = tmp;
		}
	}
	/* check the second vector */
	if (fabs(crys->b.y - (crys->b.x*crys->a.y)/crys->a.x) < TOL)
	{
		vector_t tmp = crys->b;
		crys->b = crys->c;
		crys->c = tmp;
	}
	
	/* The matrix to invert */
	invert[0][0]=crys->a.x;	invert[0][1]=crys->b.x;	invert[0][2]=crys->c.x;
	invert[1][0]=crys->a.y;	invert[1][1]=crys->b.y;	invert[1][2]=crys->c.y;
	invert[2][0]=crys->a.z;	invert[2][1]=crys->b.z;	invert[2][2]=crys->c.z;
	
	/* remove from the second row the first element */
	eliminate(invert,1,0);
	/* remove from the third row the first element */
	eliminate(invert,2,0);
	
	/* remove from the third row the second element */
	eliminate(invert,2,1);

	/* remove from the second row the third element */
	eliminate(invert,1,2);
	/* remove from the first row the third element */
	eliminate(invert,0,2);

	/* remove from the first row the second element */
	eliminate(invert,0,1);

	/* Now to divide the rows in the transformation matrix */
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			transform[i][j] /= invert[i][i];

			/* Make sure there are no NaNs */
			/* if (transform[i][j] != transform[i][j]) */
			if (isnan(transform[i][j]))
			{
				fail("Vectors are not linearly independent");
			}
		}
	}
}

/** Row subtraction for Gauss-Jordan elimination.
  * 	row and col are the cell to be eliminated
  * 	col is also the row to use for the elimination */
void eliminate(double matrix[][3],int row,int col)
{
	int i;
	double mult;
	
	mult = matrix[row][col]/matrix[col][col];
	for (i=0;i<3;i++)
	{
		matrix[row][i] -= mult*matrix[col][i];
		transform[row][i] -= mult*transform[col][i];
	}
}

/** Change the units used for input and output.
  * s should be either "ang" for Angstroms, "au" for atomic units or a floating
  * point conversion factor from Angstroms. */
void units(char * s)
{
	double new,fact;
	int i;
	/* Will change to angstroms */
	if (!strcmp(s,"ang"))
	{
		new = 1.00;
	}
	/* Will change to atomic units */
	else if (!strcmp(s,"au"))
	{
		new = 1.88972613;
	}
	/* Different conversion factor, make sure it reads something. If i is 0 or 
	 * EOF, there has been an error. */
	else if (!(i=sscanf(s,"%lf",&new)) || (i == EOF))
	{
		fail("%s not \"au\", \"ang\" or a floating point number.",s);
	}
	
	/* The conversion factor from what is currently being used to the new 
	 * units. */
	fact = new/unit;
	
	unit = new;

	/* convert the radii */
	i=0;

	/* List ends in NaN. */
	/* while(radii[i] == radii[i]) */
	while(!isnan(radii[i]))
	{
		radii[i] *= fact;
		i++;
	}
}

/** Print help text */
void help()
{
	fprintf(stderr,"\nArguments:\n");
	fprintf(stderr," -h\tPrint help text.\n");
	fprintf(stderr," -r\tChange the value of an atomic radius.\n");
	fprintf(stderr,"\tDefaults in radii.h\n");
	fprintf(stderr," -i\tSelect input  file. Defaults to stdin.\n");
	fprintf(stderr," -o\tSelect output file. Defaults to stdout.\n");
	fprintf(stderr," -a\tSelect algorithm to use.\n");
	fprintf(stderr,"\tDefaults to tetrafill.\n");
	fprintf(stderr," -v\tChange the radius for virtual particles.\n");
	fprintf(stderr,"\tDefault 1 Angstrom.\n");
	fprintf(stderr,"\nfor more information, see README\n\n");
	exit(EXIT_SUCCESS);
}
