#include "calculateEvec.hpp"

void calculateEvec(LSMSSystemParameters &lsms, LocalTypeInfo &local)
{
  Real tolerance = 1.0e-8;

  for (int i=0; i<local.num_local; i++)
  {
    // Small part in mufind_c.f from LSMS 1.9
    if (lsms.n_spin_cant == 2)         // nspin >=3
    {
      // Calculate moment
      Real moment[3];
      moment[0] = local.atom[i].dosckint[1] + local.atom[i].evec[0] * local.atom[i].mcpsc_mt;
      moment[1] = local.atom[i].dosckint[2] + local.atom[i].evec[1] * local.atom[i].mcpsc_mt;
      moment[2] = local.atom[i].dosckint[3] + local.atom[i].evec[2] * local.atom[i].mcpsc_mt;
   
      // getevec.f from LSMS 1.9
      // Determine evecNew according to moment  
      Real evecMagnitude = std::sqrt(moment[0] * moment[0] + \
                                     moment[1] * moment[1] + \
                                     moment[2] * moment[2]);
      if (lsms.global.iprint > 0)
      {
        printf(" GETEVEC: moment = (%12.8f, %12.8f, %12.8f) magnitude = %12.8f\n", moment[0], moment[1], moment[2], evecMagnitude);
      }
  
      if (evecMagnitude > tolerance)
      {
        /*
        =================================================================
        evecNew is the new moment orientation:
             evecNew[0] = the x-component of e-vector
             evecNew[1] = the y-component of e-vector
             evecNew[2] = the z-component of e-vector
        it is determined by the total moment inside the muffin-tin sphere
        =================================================================
        */
        local.atom[i].evecNew[0] = moment[0] / evecMagnitude;
        local.atom[i].evecNew[1] = moment[1] / evecMagnitude;
        local.atom[i].evecNew[2] = moment[2] / evecMagnitude; 
      }
      else
      {
        local.atom[i].evecNew[0] = local.atom[i].evec[0];
        local.atom[i].evecNew[1] = local.atom[i].evec[1];
        local.atom[i].evecNew[2] = local.atom[i].evec[2];
      }
    }
    else                               // nspin = 1 or 2
    {
      local.atom[i].evecNew[0] = local.atom[i].evec[0];
      local.atom[i].evecNew[1] = local.atom[i].evec[1];
      local.atom[i].evecNew[2] = local.atom[i].evec[2];
    }

    /*
    ================================================================
    Store direction & mag. mom. corresponding to output chg. den....
    ================================================================
    Not yet implemented. (need to see where evec_out and mm_out are used for)
    */

    local.atom[i].evecOut[0] = local.atom[i].evecNew[0];
    local.atom[i].evecOut[1] = local.atom[i].evecNew[1];
    local.atom[i].evecOut[2] = local.atom[i].evecNew[2];

  }
  
  return;

}


void mixEvec(LSMSSystemParameters &lsms, LocalTypeInfo &local, Real alpev)
{

  /*
  ================================================================
  perform simple mixing of evec_new and evec_old, and redefine
  evec_new........................................................
  ================================================================
  !! This should be placed in mixing.hpp !!
  */

  Real tolerance = 1.0e-8;

  for (int i=0; i<local.num_local; i++)
  { 
 
    if (lsms.global.iprint > 0)
    {
      printf("Moment direction before mixing = (%12.8f, %12.8f, %12.8f)\n", local.atom[i].evecNew[0], local.atom[i].evecNew[1], local.atom[i].evecNew[2]);
    }
  
    local.atom[i].evecNew[0] = alpev * local.atom[i].evecNew[0] + (1.0-alpev) * local.atom[i].evec[0];
    local.atom[i].evecNew[1] = alpev * local.atom[i].evecNew[1] + (1.0-alpev) * local.atom[i].evec[1];
    local.atom[i].evecNew[2] = alpev * local.atom[i].evecNew[2] + (1.0-alpev) * local.atom[i].evec[2];
  
    Real evecMagnitude = std::sqrt(local.atom[i].evecNew[0] * local.atom[i].evecNew[0] + \
                                   local.atom[i].evecNew[1] * local.atom[i].evecNew[1] + \
                                   local.atom[i].evecNew[2] * local.atom[i].evecNew[2]);
  
    if (evecMagnitude < tolerance)
    {
      printf("GETEVEC: magnitude of evec too small. (= %35.25f)\n", evecMagnitude);
    }
  
    local.atom[i].evecNew[0] = local.atom[i].evecNew[0] / evecMagnitude;
    local.atom[i].evecNew[1] = local.atom[i].evecNew[1] / evecMagnitude;
    local.atom[i].evecNew[2] = local.atom[i].evecNew[2] / evecMagnitude; 
  
    if (lsms.global.iprint > 0)
    {
      printf("Moment direction after mixing = (%12.8f, %12.8f, %12.8f)\n", local.atom[i].evecNew[0], local.atom[i].evecNew[1], local.atom[i].evecNew[2]);
    } 

  }

  return;

}



