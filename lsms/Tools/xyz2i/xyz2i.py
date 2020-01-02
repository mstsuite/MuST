#!/usr/bin/env python

def write_b1_header(f,bravais,id,n):
    f.write("*******************************************************************************\n\
 *                        Input parameter data file                            *\n\
 *                      BIGCELL code (version 1.7)                             *\n\
 *******************************************************************************\n\
                                                                               *\n\
 cpath  : path-name-> input : output : keep : info : etc : files:______________*\n\
./                                                                             F\n\
 cfspath: path-name-> v_gopen_xxx : w_gopen_xxx : files:_______________________*\n\
./                                                                             F\n\
 Input file: In general this should be the name of this file___________________*\n\
i_bigcell                                                                      F\n\
 Text to identify system: (used to construct file names)_______________________*\n")
    f.write("{0:<79}F\n".format(id))
    f.write(" subroutine stop level [a10 format] lower case_________________________________*\n\
main                                                                           F\n\
 print level___[ipr=0 gets min. O/P]___________________________________________*\n\
           0                                                                   U\n\
 nbasis___[# sublattices andi/or # atoms in system]____________________________*\n")
    f.write("{0:<79}U\n".format(n))
    f.write(" Spin index__[1,2,3 -> non-spin-polar;spin-polar;spin-canting]_________________*\n\
          3                                                                    U\n\
__Bx.123456789012345__By.123456789012345__Bz.123456789012345_Brav. Lat. [a.u.]_*\n")
    f.write(" {0[0]:>19.15f} {0[1]:>19.15f} {0[2]:>19.15f}\n".format(bravais[0]))
    f.write(" {0[0]:>19.15f} {0[1]:>19.15f} {0[2]:>19.15f}\n".format(bravais[1]))
    f.write(" {0[0]:>19.15f} {0[1]:>19.15f} {0[2]:>19.15f}\n".format(bravais[2]))

def write_b1_positions(f,xyz,Z,pot,lmax,rLIZ,rST1,rST2,rST3,rST4,n):
    f.write(" |nlines___|seed___|ranp_mode___|seed_cant_____________________________________*\n")
    f.write(  "{0:>8}     0.0            0          0.0                                     U\n".format(n))
    f.write(" Nm___x.123456789012345___y.123456789012345___z.123456789012340__c.1234_i2345__*\n")
    for i in range(n):
        f.write("{0:>3} {1[0]:>19.15f} {1[1]:>19.15f} {1[2]:>19.15f}  1.0000     0  U\n".format(Z,xyz[i]))
    f.write("_Nm___|vfname_(1a15)|__lmax__r_liz7__r_st1__r_st2__r_st3__r_st4__met_rad_______*\n")
    for i in range(n):
        f.write("{0:>3}   {1:<15}\n  {2:>4}  {3:>6.3f}  {4:>5.2f}  {5:>5.2f}  {6:>5.2f}  {7:>5.2f}    2.000   F\n".format(Z,pot,lmax,rLIZ,rST1,rST2,rST3,rST4))

def write_b1_tail(f):
    f.write("  N_x   N_y   N_z  | # repeats of basic unit cell in x,y,z directions__________*\n\
    1     1     1                                                              U\n\
 |# sublats remodeled  |temp_max  |num_temp_steps  |max_picks__________________*\n\
  0                     1000       10               25000                      U\n\
 |sub-lattice index  |# SRO params read  |SRO parameters_______________________*\n\
  0                   0                   0.000  0.000  0.000  0.000  0.000    U\n\
                                                                               %\n\
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\
                                                                               %\n\
 NOTES:....                                                                    %\n\
                                                                               %\n\
 1)    nlines:    no. of lines of data entries following the line              %\n\
                  \"Nm___x.1234...\".                                            %\n\
 2)      seed:    the initial seed value of the random number generator for    %\n\
                  generating random sample.                                    %\n\
 3) ranp_mode: = 0, generate random sample without going through               %\n\
                    sublattice by sublattice                                   %\n\
               = 1, generate random sample on the basis of sublattices         %\n\
 4)   seed_cant: the initial seed value of the random number generator for     %\n\
                  generating moment orientations. If < 0, a moment pointing    %\n\
                  to (0, 0, 1) direction is generated.                         %\n\
 5)   to re-modelling whole system at once rather than going through           %\n\
      each sub-lattice, set sub-lattice index to be zero                       %\n\
                                                                               %\n\
 6) Above F indicates a line that is formatted                                 %\n\
          U indicates a line that is unformatted                               %\n\
          * indicates a text line used to indicate function of following lines %\n\
          = indicates a line that is not read in : Used for coments            %\n\
                                                                               %\n\
1234567890123456789012345678901234567890123456789012345678901234567890123456789%\n\
 ==============================================================================%\n")

def write_lsms_1(filename,id,n):
    f = open(filename,"w")
    f.write(" ****************************************************************************1\n\
 *                        Input parameter data file                         *2\n\
 *                         LSMS code (version 1.5)                          *3\n\
 ****************************************************************************4\n\
                                                                            *5\n\
 cpath  : path-name-> input : output : keep : info : etc : files:___________*6\n\
./                                                                          *7\n\
 cfspath: path-name-> v_gopen_xxx : w_gopen_xxx : files:____________________*8\n\
./                                                                          *9\n\
 Input file for Real Space SCF Code_________________________________________*0\n\
i_lsms                                                                      *1\n\
 Text to identify system: (used to construct file names)____________________*2\n\
{0:<76}*3\n\
 output_to_screen [= y, yes; = n, no, will output to a file]________________*4\n\
n                                                                           *5\n\
 subroutine stop level [a10 format] lower case______________________________*6\n\
main                                                                        *7\n\
 node_print, print_instr, nprint____________________________________________*8\n\
  0 -1 -1                                                                   *9\n\
 specify number of atoms in the system______________________________________*0\n\
  {1:>7}                                                                   *1\n\
 nrelc,nrelv,mtasa [0(>1)=Rel(Non)-Rel; 0(>1)=Scalar(Non)-Rel; 0(1)=MT(ASA)]*2\n\
 10 10   1 0.00000   0.0000   0.0000                                        *3\n\
 nspin [1=>para, 2=>spin-pol, 3=>spin-cant] : i_vdif [0=>vdif=zero, 1=>vdif]*4\n\
  3                                            0                    1       *5\n\
 Text to identify system: (used to construct file names)____________________*6\n\
 Title Text Not Specified                                                   *7\n\
 read in ngaussr & ngaussq, the no. of Gaussian pnts for volume integration *8\n\
 10   40                                                                    *9\n\
|info_table name (a30)________|info_evec name (a30) [or \"default\"]__________*0\n\
info_table_{0:<15}    info_evec_{0:<15}                     *1\n\
 igrid : ebot :  etop  :  eitop  : eibot : npts : kelvin : nument : ialg :__*4\n\
   2   -0.30000 0.00000  0.82500 +0.00250  031     0000      00     02      *5\n\
 nscf : alpdv : alpma : alpev :  pot :  dga : iharris : i_potwrite : movie__*6\n\
   10   0.000   0.000   0.000     0      2       0           1         0    *7\n\
 ntstep : tstep :   etol     :  ptol   :  eftol  :  rmstol  : SD & SCF tols_*6 \n\
   01     1.000   0.00000005   0.00001   0.00001   0.0000001                *7\n\
 ctq : j_ij_(ctq: The coefficient for torque [ >= 0.0 ], j_ij = 0 or 1)_____*8\n\
 0.0   0                                                                    *9\n\
 No. of atoms whose mixing parameters and scheme need to be adjusted [>= 0]_*0\n\
 0                                                                          *1\n\
 atom name__node number___alphadv___alphama___alphaev___pot___dga___________*2\n\
 Co         000           0.90000   0.50000   0.00000    1     0            *3\n\
 No. of atoms whose Zc, Zs, and Zv need to be adjusted [>= 0]_______________*4\n\
 0                                                                          *5\n\
 atom name__node number____Zcore__Zsemi__Zvale [Zcore+Zsemi+Zvale = Zatom]__*6\n\
 Co         000            10.50  07.50  10.00                              *7\n\
     etol      ptol      eftol      rmstol  tolerances for selfconsistency:-*8\n\
 01 1.000 0.0000005   0.0001    0.0001    0.000001                          *9\n\
 if igrid=1 and npts>300, no. extra energy points to be read in: npts-300___*0\n\
 extra energy points(in complex format):____________________________________*1\n\
 (0.7520000D+00, 0.3000000D+00)                                             *2\n\
 ****************************************************************************9\n\
 *                           End of the File                                *0\n\
 ****************************************************************************1\n".format(id,n))
    f.close()

def write_bigcell_1(filename,bravais,xyz,id,Z,pot,lmax,rLIZ,rST1,rST2,rST3,rST4,n):
    f = open(filename,"w")
    write_b1_header(f,bravais,id,n)
    write_b1_positions(f,xyz,Z,pot,lmax,rLIZ,rST1,rST2,rST3,rST4,n)
    write_b1_tail(f)
    f.close()

def read_xyz(filename,scale,bravais):
    f=open(filename,"r")
    lines=f.readlines()
    f.close()
    xyz=[]
    for line in lines:
        recs=line.split()
        x=scale*float(recs[0])
        y=scale*float(recs[1])
        z=scale*float(recs[2])
        xyz.append((x,y,z))
    return xyz

def main():
    id="fe432"
    n=2
    scale=1.0/0.5291772
    a=2.8665
    nx=6
    ny=6
    nz=6
    bravais=[(nx*a*scale,0.0,0.0),
             (0.0,ny*a*scale,0.0),
             (0.0,0.0,nz*a*scale)]
    xyz=read_xyz("pos.dat",scale,bravais)
    n=len(xyz)

    write_bigcell_1("i_bigcell",bravais,xyz,id,26,"fe_v",3,12.5,89.0,90.0,91.0,92.0,n)
    write_lsms_1("i_lsms",id,n)

if __name__ == "__main__" :
  main()
