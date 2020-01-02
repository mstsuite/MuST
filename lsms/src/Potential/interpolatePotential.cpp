#include "interpolatePotential.hpp"

/*  This is equivalent to the new grid setup part of pot_adapt.f in LSMS 1.9. */
/*  However, the interpolations of potential and charge density will not be   */
/*  performed if the new mesh is the same as the old mesh.                    */


// Warning: ASA case is not taken care of yet!!

void interpolatePotential(LSMSSystemParameters &lsms, AtomData &atom)
{

  // Generate new radial mesh according to the new Muffin-Tin radius (inscribed radius calculated in setupVorpol)
  std::vector<Real> r_mesh_old;
  Matrix<Real> vr_old, rhotot_old;
  int jmt0 = atom.r_mesh.size();
  int jws0 = atom.jws;
  int f = 0;
  Real dvr;                        // Derivative of vr

  r_mesh_old.resize(jmt0);
  vr_old.resize(jmt0, 2);
  rhotot_old.resize(jmt0, 2);

  r_mesh_old = atom.r_mesh;
  vr_old = atom.vr;
  rhotot_old = atom.rhotot;

  atom.generateRadialMesh();

  // Find the new jws after new mesh is defined
  for (int ir = 0; ir < lsms.global.iprpts; ir++) {
    if (atom.r_mesh[ir] < atom.rws && lsms.mtasa == 0)
      atom.jws = ir + 2;
    if (atom.jws > lsms.global.iprpts)
    {
      printf("Problem! jws > iprpts. jws = %8d iprpts = %8d\n", atom.jws, lsms.global.iprpts);
      printf("         rmt = %lf  rws = %lf last r_mesh = %lf\n",
        atom.rmt, atom.rws, atom.r_mesh[lsms.global.iprpts-1]);
    }
  }

  // Interpolate potential onto the new mesh
  for (int is = 0; is < lsms.n_spin_pola; is++)
  {
    for (int ir = 0; ir < atom.jmt; ir++)
    {
      if (atom.r_mesh[ir] < r_mesh_old[jmt0-1])
        interp_(&r_mesh_old[0], &vr_old(0,is), &jmt0, &atom.r_mesh[ir], &atom.vr(ir,is), &dvr, &f);
      else
        atom.vr(ir,is) = vr_old(jmt0-1,is) * atom.r_mesh[ir] / r_mesh_old[jmt0-1];
        //YingWai's check
        //printf("ir = %8d r_mesh_old = %35.25f r_mesh = %35.25f vr_old = %35.25f vr = %35.25f\n", ir, r_mesh_old[ir], atom.r_mesh[ir], vr_old(ir,is), atom.vr(ir,is));
    }
  }

    // Interpolate charge density onto the new mesh
  if (jws0 > 0) {
    for (int is = 0; is < lsms.n_spin_pola; is++)
    {
      for (int ir = 0; ir < atom.jmt; ir++)
      {
        if (atom.r_mesh[ir] < r_mesh_old[jmt0-1])
          interp_(&r_mesh_old[0], &rhotot_old(0,is), &jmt0, &atom.r_mesh[ir], &atom.rhotot(ir,is), &dvr, &f);
        else
          atom.rhotot(ir,is) = rhotot_old(jmt0-1,is);
          //YingWai's check
          //printf("ir = %8d r_mesh_old = %30.25f r_mesh = %35.25f rhotot_old = %35.25f rhotot = %30.25f\n", ir, r_mesh_old[ir], atom.r_mesh[ir], rhotot_old(ir,is), atom.rhotot(ir,is));
      }
    }
  }

return;

}

