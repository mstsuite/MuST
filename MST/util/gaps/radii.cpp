/**	radii.c:
 *	
 * Empirical values for atomic radii in crystals. From: 
 * 	Slater, J C. "Atomic Radii in Crystals." Journal of Chemical Physics. 
 * 		41.10 (1964): Print.
 *  
 * Value of 0 means its radius is unknown in a crystalline structure.
 * it is usually because it is a noble gas or too unstable to measure.
 * 
 * The name in comments after each element refers to which of the above 
 * publications the value was obtained from.
 * 
 * George Bargoud: <gbargoud@gmail.com>
 **/
#include <math.h>

#ifndef NAN
#define NAN (0.0/0.0)
#endif
/* Atomic radii in Angstroms */
double radii[] = 
{
/* Virtual */  1.00, /* Virtual particle used in 'gaps' */
/* H  */ 0.25, /* Slater */
/* He */ 0.00,
/* Li */ 1.45, /* Slater */
/* Be */ 1.05, /* Slater */
/* B  */ 0.85, /* Slater */
/* C  */ 0.70, /* Slater */
/* N  */ 0.65, /* Slater */
/* O  */ 0.60, /* Slater */
/* F  */ 0.50, /* Slater */
/* Ne */ 0.00,
/* Na */ 1.80, /* Slater */
/* Mg */ 1.50, /* Slater */
/* Al */ 1.25, /* Slater */
/* Si */ 1.10, /* Slater */
/* P  */ 1.00, /* Slater */
/* S  */ 1.00, /* Slater */
/* Cl */ 1.00, /* Slater */
/* Ar */ 0.00,
/* K  */ 2.20, /* Slater */
/* Ca */ 1.80, /* Slater */
/* Sc */ 1.60, /* Slater */
/* Ti */ 1.40, /* Slater */
/* V  */ 1.35, /* Slater */
/* Cr */ 1.40, /* Slater */
/* Mn */ 1.40, /* Slater */
/* Fe */ 1.40, /* Slater */
/* Co */ 1.35, /* Slater */
/* Ni */ 1.35, /* Slater */
/* Cu */ 1.35, /* Slater */
/* Zn */ 1.35, /* Slater */
/* Ga */ 1.30, /* Slater */
/* Ge */ 1.25, /* Slater */
/* As */ 1.15, /* Slater */
/* Se */ 1.15, /* Slater */
/* Br */ 1.15, /* Slater */
/* Kr */ 0.00,
/* Rb */ 2.35, /* Slater */
/* Sr */ 2.00, /* Slater */
/* Y  */ 1.80, /* Slater */
/* Zr */ 1.55, /* Slater */
/* Nb */ 1.45, /* Slater */
/* Mo */ 1.45, /* Slater */
/* Tc */ 1.35, /* Slater */
/* Ru */ 1.30, /* Slater */
/* Rh */ 1.35, /* Slater */
/* Pd */ 1.40, /* Slater */
/* Ag */ 1.60, /* Slater */
/* Cd */ 1.55, /* Slater */
/* In */ 1.55, /* Slater */
/* Sn */ 1.45, /* Slater */
/* Sb */ 1.45, /* Slater */
/* Te */ 1.40, /* Slater */
/* I  */ 1.40, /* Slater */
/* Xe */ 0.00,
/* Cs */ 2.60, /* Slater */
/* Ba */ 2.15, /* Slater */
/* La */ 1.95, /* Slater */
/* Ce */ 1.85, /* Slater */
/* Pr */ 1.85, /* Slater */
/* Nd */ 1.85, /* Slater */
/* Pm */ 1.85, /* Slater */
/* Sm */ 1.85, /* Slater */
/* Eu */ 1.85, /* Slater */
/* Gd */ 1.80, /* Slater */
/* Tb */ 1.75, /* Slater */
/* Dy */ 1.75, /* Slater */
/* Ho */ 1.75, /* Slater */
/* Er */ 1.75, /* Slater */
/* Tu */ 1.75, /* Slater */
/* Yb */ 1.75, /* Slater */
/* Lu */ 1.75, /* Slater */
/* Hf */ 1.55, /* Slater */
/* Ta */ 1.45, /* Slater */
/* W  */ 1.35, /* Slater */
/* Re */ 1.35, /* Slater */
/* Os */ 1.30, /* Slater */
/* Ir */ 1.35, /* Slater */
/* Pt */ 1.35, /* Slater */
/* Au */ 1.35, /* Slater */
/* Hg */ 1.50, /* Slater */
/* Tl */ 1.90, /* Slater */
/* Pb */ 1.80, /* Slater */
/* Bi */ 1.60, /* Slater */
/* Po */ 1.90, /* Slater */
/* At */ 0.00,
/* Rn */ 0.00,
/* Fr */ 0.00,
/* Ra */ 2.15, /* Slater */
/* Ac */ 1.95, /* Slater */
/* Th */ 1.80, /* Slater */
/* Pa */ 1.80, /* Slater */
/* U  */ 1.75, /* Slater */
/* Np */ 1.75, /* Slater */
/* Pu */ 1.75, /* Slater */
/* Am */ 1.75, /* Slater */
NAN /* NAN represents the end of the list */
};
