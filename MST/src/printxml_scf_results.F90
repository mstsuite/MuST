subroutine printxml_scf_results()
   use KindParamModule, only : IntKind, RealKind
   use PhysParamModule, only : Bohr2Angstrom, Ryd2eV
   use XMLDataModule, only : writeStartTag, writeEndTag, writeElement
   use SystemModule, only : getNumAtoms
   use ForceModule, only : isForceAvailable, getForce
   use TotalEnergyModule, only : getEnergyPerAtom
!
   implicit none
!
   character (len=54) :: force_vec
   character (len=18) :: energy_val
!
   integer (kind=IntKind) :: n, NumAtoms
!
   real (kind=RealKind) :: f(3), df(3)
   real (kind=RealKind) :: AU2SI
   real (kind=RealKind) :: en_band, en_coul, en_xc, en_total
!
   NumAtoms = getNumAtoms()
!
   AU2SI = Ryd2eV/Bohr2Angstrom
!
   call writeStartTag(tag='scf_results')
!     ===========================================
      if ( isForceAvailable() ) then
         call writeStartTag(tag='varray',ename='forces')
         do n = 1, NumAtoms
            f = getForce(df,global_id=n)
            write(force_vec,'(3f18.8)')(f(1:3)-df(1:3))*AU2SI
            call writeElement(tag='v',econtent=force_vec)
         enddo
         call writeEndTag(tag='varray')
      endif
!     ===========================================
      call writeStartTag(tag='energy')
         en_total = getEnergyPerAtom(e_band=en_band,e_coul=en_coul,e_xc=en_xc)
         write(energy_val,'(f18.8)')en_total
         call writeElement(tag='i',ename='energy_total',econtent=energy_val)
         write(energy_val,'(f18.8)')en_band
         call writeElement(tag='i',ename='energy_band',econtent=energy_val)
         write(energy_val,'(f18.8)')en_coul
         call writeElement(tag='i',ename='energy_Hartree',econtent=energy_val)
         write(energy_val,'(f18.8)')en_xc
         call writeElement(tag='i',ename='energy_xc',econtent=energy_val)
      call writeEndTag(tag='energy')
!     ===========================================
   call writeEndTag(tag='scf_results')
!
end subroutine printxml_scf_results
