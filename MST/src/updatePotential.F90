subroutine updatePotential(LocalNumAtoms,n_spin_pola)
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, SQRTm1, CZERO
!
   use RadialGridModule, only : getNumRmesh
!
   use PotentialTypeModule, only : isFullPotential
!
   use PotentialModule, only : setVdif, setSphPotr, setPotential
!
   use AtomModule, only : getLocalNumSpecies
!
   use PotentialGenerationModule, only : getNewSphPotr => getSphPotr
   use PotentialGenerationModule, only : getNewPotential => getPotential
   use PotentialGenerationModule, only : getVdif, getPotLmax
   use PotentialGenerationModule, only : getPotComponentFlag
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: n_spin_pola
!
   integer (kind=IntKind) :: id, ia, is, NumRs, jl, lmax, ir, jmax
   integer (kind=IntKind), pointer :: flag(:)
!
   real (kind=RealKind) :: pot_r
   real (kind=RealKind), pointer :: new_SphPotr(:)
   real (kind=RealKind), pointer :: p_vdif(:)
!
   complex (kind=CmplxKind), pointer :: new_PotentialL(:)
!
   if ( .not.isFullPotential() ) then
      do id = 1, LocalNumAtoms
         NumRs = getNumRmesh(id)
         do ia = 1, getLocalNumSpecies(id)
            do is = 1, n_spin_pola
               new_SphPotr => getNewSphPotr(id,ia,is)
!              -------------------------------------------------------
               call setSphPotr(id,ia,is,new_SphPotr)
!              -------------------------------------------------------
            enddo
         enddo
      enddo
!
      if (n_spin_pola == 2) then
!        -------------------------------------------------------------
         p_vdif => getVdif()
         call setVdif(p_vdif(1))
!        -------------------------------------------------------------
      endif
   else
      do id = 1, LocalNumAtoms
         lmax = getPotLmax(id)
         jmax = ((lmax+1)*(lmax+2))/2
         NumRs = getNumRmesh(id)
         flag => getPotComponentFlag(id)
         if ( maxval(flag) <=1 ) then
            do ia = 1, getLocalNumSpecies(id)
               do is = 1, n_spin_pola
                  do jl = 1,jmax
                     new_PotentialL => getNewPotential("Total",id,ia,is,jl)
                     if ( flag(jl)==0 ) then
                        new_PotentialL = CZERO
                     endif
!                    -------------------------------------------------
                     call setPotential(id,ia,is,jl,new_PotentialL,flag(jl))
!                    -------------------------------------------------
                  enddo
               enddo
            enddo
         else
            do ia = 1, getLocalNumSpecies(id)
               do is = 1, n_spin_pola
                  do jl = 1,jmax
                     new_PotentialL => getNewPotential("Total",id,ia,is,jl)
                     if ( flag(jl)==1 ) then
                        do ir =1,NumRs
                           pot_r = real(new_potentialL(ir),kind=RealKind)
                        enddo
                     else if ( flag(jl)==2 ) then
                        do ir =1,NumRs
                           pot_r = real(-sqrtm1*new_potentialL(ir),kind=RealKind)
                           new_potentialL(ir) = cmplx(ZERO,pot_r,kind=CmplxKind)
                        enddo
                     else if ( flag(jl)==0 ) then
                        new_PotentialL = CZERO
                     endif
!                    -------------------------------------------------
                     call setPotential(id,ia,is,jl,new_PotentialL,flag(jl))
!                    -------------------------------------------------
                  enddo
               enddo
            enddo
         endif
      enddo
!
      if (n_spin_pola == 2) then
!        -------------------------------------------------------------
         p_vdif => getVdif()
         call setVdif(p_vdif(1))
!        -------------------------------------------------------------
      endif
   endif
!
end subroutine updatePotential
