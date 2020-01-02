!  *******************************************************************
!  *  inserts a vector in a list of vectors such that they are in 
!  *  order of increasing length.
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine ord3v(v3out_1,v3out_2,v3out_3,vsqout,nv3,v3in,vsqin)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
   implicit  none
!
   integer (kind=IntKind), intent(inout) :: nv3
   integer (kind=IntKind) :: nv
   integer (kind=IntKind) :: nvm
!
   real (kind=RealKind), intent(inout) :: v3out_1(nv3+1)
   real (kind=RealKind), intent(inout) :: v3out_2(nv3+1)
   real (kind=RealKind), intent(inout) :: v3out_3(nv3+1)
   real (kind=RealKind), intent(inout) :: vsqout(nv3+1)
   real (kind=RealKind), intent(in) :: v3in(3)
   real (kind=RealKind), intent(in) :: vsqin
!
!  ===================================================================
   do nv=1,nv3
      if(vsqout(nv).ge.vsqin) then
         do nvm=nv3,nv,-1
            vsqout(nvm+1)=vsqout(nvm)
         enddo
         do nvm=nv3,nv,-1
            v3out_1(nvm+1) = v3out_1(nvm)
         enddo
         do nvm=nv3,nv,-1
            v3out_2(nvm+1) = v3out_2(nvm)
         enddo
         do nvm=nv3,nv,-1
            v3out_3(nvm+1) = v3out_3(nvm)
         enddo
         vsqout(nv)=vsqin
         v3out_1(nv) = v3in(1)
         v3out_2(nv) = v3in(2)
         v3out_3(nv) = v3in(3)
         nv3=nv3+1
         return
      endif
   enddo
   vsqout(nv3+1)=vsqin
   v3out_1(nv3+1) = v3in(1)
   v3out_2(nv3+1) = v3in(2)
   v3out_3(nv3+1) = v3in(3)
   nv3=nv3+1
!
   end subroutine ord3v
!  ===================================================================
