/** fill.h 
 * 
 * prototypes for the filling algorithms and the mad from their names
 * 
 * George Bargoud: <gbargoud@gmail.com>
 **/

#ifndef FILL_H
#define FILL_H

/********************\
* Filling algorithms * 
\********************/
/** Builds a 3d grid and places a virtual atom on each vertex and then drifts 
  * it. */
void gridfill(crystal_t * crystal);
/***********************\
* Filling algorithm map *
\***********************/
struct {
	char * name;
	fill_t func;
} funcmap [] = {
	{"gridfill"		,gridfill		},
	{NULL,NULL}
};

#endif /* FILL_H */
