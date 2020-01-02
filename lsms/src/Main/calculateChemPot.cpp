/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "calculateChemPot.hpp"

// calculate the chemical potential and the eigenvalue sum
// see LSMS_1.9 mufind
void calculateChemPot(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local, Real &eigensum)
{
  Real wspace[4];

/*
c     ================================================================
c     calculates ::: integrated n(e) .................................
c     calculate WS valence charges up to etop => chem.pot.............
c     xtws is obtained from a global sum..............................
c     ================================================================
*/

  Real xtws = 0.0, tnen = 0.0;
  for (int i=0; i<local.num_local; i++)
  {
    if (lsms.n_spin_cant == lsms.n_spin_pola) // either non spin polarized or spin canted
    {
      xtws += local.atom[i].dosint[0] * Real(local.n_per_type[i]);
      tnen += local.atom[i].doslast[0] * Real(local.n_per_type[i]);
    }
    else // spin polarized not spin canted
    {
      xtws += (local.atom[i].dosint[0] + local.atom[i].dosint[1]) * Real(local.n_per_type[i]);
      tnen += (local.atom[i].doslast[0] + local.atom[i].doslast[1]) * Real(local.n_per_type[i]);
    }
  }

    wspace[0] = xtws;
    wspace[1] = tnen;

    globalSum(comm, wspace, 2);

    xtws = wspace[0];
    tnen = wspace[1];

// if(iharris<=1)
    Real etop = lsms.chempot;
    if (lsms.global.iprint >= 0)
    {
      printf("atom[0].dosint[0] = %lf  atom[0].dosint[1] = %lf\n",local.atom[0].dosint[0], local.atom[0].dosint[1]);
      if(lsms.n_spin_cant==2)
	printf("atom[0].dosint[2] = %lf  atom[0].dosint[3] = %lf\n",local.atom[0].dosint[2], local.atom[0].dosint[3]);
      printf("calculateChempot: xtws = %lf zvaltss = %lf tnen = %lf old chempot = %lf\n",
             xtws, lsms.zvaltss, tnen, lsms.chempot);
    }
//    lsms.chempot = lsms.energyContour.etop + (lsms.zvaltss-xtws)/tnen;
    lsms.chempot = etop + (lsms.zvaltss - xtws) / tnen;
    // restrict the change in the chemical potential and prevent it from dropping below ebot
    if(std::abs(lsms.chempot-etop)>0.1)
    {
      lsms.chempot = etop + 0.1 *  (lsms.zvaltss - xtws) / std::abs(lsms.zvaltss - xtws);
      if (lsms.global.iprint >= 0)
	printf("           chempot stepsize restrcited!\n");
    }
    if(lsms.chempot < lsms.energyContour.ebot)
    {
      lsms.chempot = 0.5*(lsms.energyContour.ebot + etop);
      if (lsms.global.iprint >= 0)
	printf("           chempot limited by botom of contour!\n");
    }
    
    if (lsms.global.iprint >= 0)
      printf("           new chempot = %lf\n", lsms.chempot);
    Real chempot = lsms.chempot;
//    Real etop = lsms.energyContour.etop;
/*
c     ================================================================
c     compute integrated densities of states, torques and stoner 
c     parameters up to fermi energy...................................
c     ================================================================
*/
    eigensum = 0.0;
    for (int i=0; i<local.num_local; i++)
    {
      for (int is=0; is<lsms.n_spin_pola+2*(lsms.n_spin_cant-1); is++)
      {
        local.atom[i].dosint[is] += local.atom[i].doslast[is] * (chempot - etop);
        local.atom[i].dosckint[is] += local.atom[i].doscklast[is] * (chempot - etop);
        local.atom[i].evalsum[is] += 0.5 * (chempot + etop) * local.atom[i].doslast[is] * (chempot - etop);
/*
//evalsum(is)=evalsum(is)+
//     &               (half*(chempot+etop)-(is-1)*vdif*i_vdif)*
//     &               doslast(is)*(chempot-etop)
*/
        for (int ir=0; ir < local.atom[i].jws; ir++)
          local.atom[i].greenint(ir,is) += local.atom[i].greenlast(ir,is) * (chempot - etop);
      }

/*
c     ================================================================
c     calculate contibution to eigenvalue sum from spin-up and spin-
c     down bands, and average D.O.S. at the last point on energy mesh.
c     ================================================================
 */
      if (lsms.n_spin_cant == 2)
      {
        Real d0 = local.atom[i].doslast[0];
        Real dm = local.atom[i].doslast[1] * local.atom[i].evec[0] +
                  local.atom[i].doslast[2] * local.atom[i].evec[1] +
                  local.atom[i].doslast[3] * local.atom[i].evec[2];
        Real ev0 = local.atom[i].evalsum[0];
        Real evm = local.atom[i].evalsum[1] * local.atom[i].evec[0]+
                   local.atom[i].evalsum[2] * local.atom[i].evec[1] +
                   local.atom[i].evalsum[3] * local.atom[i].evec[2];
//        avedos_ef(1)=0.5*(d0+dm)
//        avedos_ef(2)=0.5*(d0-dm)
        local.atom[i].evalsum[0] = 0.5 * (ev0 + evm);
        local.atom[i].evalsum[1] = 0.5 * (ev0 - evm);
        eigensum += ev0 * Real(local.n_per_type[i]);
      }
      else if (lsms.n_spin_pola == 2)
      {
//        avedos_ef(1)=doslast(1)
//        avedos_ef(2)=doslast(2)
        eigensum += (local.atom[i].evalsum[0] + local.atom[i].evalsum[1]) * Real(local.n_per_type[i]);
      }
      else
      {
//        avedos_ef(1)=doslast(1)
        eigensum += local.atom[i].evalsum[0] * 2.0 * Real(local.n_per_type[i]);
      }
    }
      globalSum(comm, eigensum);
/*
      // if(relativity==full) then
      if(nrel_rel.ne.0) then
        do is=1,3
          dosint_orb(is)=dosint_orb(is)+
     >                   doslast_orb(is)*(chempot-etop)
          dosckint_orb(is)=dosckint_orb(is)+
     >                   doscklast_orb(is)*(chempot-etop)
          do ir=1,jws
            densint_orb(ir,is)=densint_orb(ir,is)+
     >                   denslast_orb(ir,is)*(chempot-etop)
          end do
        end do
      end if
*/
  if(lsms.relativity==full)
  {
     // calculate orbital moments
     printf("Orbital moments not yet implemented in calulateChempPot.cpp!\n");
  }

}
