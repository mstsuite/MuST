!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBandStructure(is,eb,et,ne,kpts,nk,                   &
                               computeSingleScattering,getScatteringMatrix)
!  ===================================================================
   use StrConstModule, only : getStrConstMatrix
!
   implicit none
!
   character (len=20) :: sname = "calBandStructure"
!
   interface
      subroutine computeSingleScattering(is,site,e,vs,ss,ir,wr)
         use KindParamModule, only : IntKind, CmplxKind
         integer (kind=IntKind), intent(in) :: is,site
         complex (kind=CmplxKind), intent(in) :: e,vs
         character (len=*), intent(in), optional :: ir
         logical, intent(in), optional :: ss,wr
      end subroutine computeSingleScattering
   end interface
!
   interface
      function getSingleScatteringMatrix(smt,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in) :: site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface

!
   integer (kind=IntKind), intent(in) :: is, ne, nk
   integer (kind=IntKind) :: ik, col, row, ie, id, js, ns
!
   real (kind=RealKind), intent(in) :: kpts(3,nk)
   real (kind=RealKind) :: kvec(3), de
!
   real (kind=RealKind), intent(in) :: eb, et ! Find energy eigenvalues
                                              ! within (eb, et)
!
   complex (kind=CmplxKind), pointer :: scm(:,:)
!
   do ik = 1, nk    ! Loop over the k-point mesh
      kvec(1:3) = kpts(1:3,ik)
      write(6,'(a,3f12.6)')'kvec(1:3) = ',kvec(1:3)
!
      de = (et-eb)/real(ne-1,kind=RealKind)
      do ie = 1, ne ! Loopp over the energy mesh
         energy = eb + (ie-1)*de
         if (abs(energy) < TEN2m6) then
            energy = energy + de*HALF
         endif
         write(6,'(a,f12.6)')'energy = ',real(energy,kind=RealKind)
!
         if (isRelativistic) then !xianglin
            kappa = sqrt(2.d0*Me*energy + energy**2/LightSpeed**2)
         else
            kappa = sqrt(energy)
         endif
!
         do col = 1, LocalNumAtoms
            do row = 1, GlobalNumAtoms
!              -------------------------------------------------------
               scm => getStrConstMatrix(kvec,kappa,id_array(row),jd_array(col), &
                                        lmaxi_array(row),lmaxj_array(col))
!              -------------------------------------------------------
               sc_blocks(row,col)%strcon_matrix = scm
            enddo
         enddo
!
         do id = 1, LocalNumAtoms
            do js = 1, nSpinCant
               ns = max(js,is)
!              -------------------------------------------------------
               call computeSingleScattering(ns,id,energy,CZERO)
!              -------------------------------------------------------
            enddo
         enddo
!
!        =============================================================
!        setup S- and C- matrix in global frame
!        -------------------------------------------------------------
         call setupSCinGlobalFrame(getSingleScatteringMatrix)
!        -------------------------------------------------------------
!
         KKR_MatrixBand = CZERO
!        =============================================================
!        Compute the KKR matrix (kappa*C+B*S), which is stored in KKR_MatrixBand
!        -------------------------------------------------------------
         call computeKKRMatrix()
!        -------------------------------------------------------------
!
!        =============================================================
!        At this point, KKR_MatrixBand is calculated. Note:
!           det[ KKR_Matrix ] = CZERO gives rise to the band structure
!           KKR_Matrix is a KKRMatrixSizeCant x KKRMatrixSizeCant 
!           matrix, and is divided into multiple Bands of columns so
!           that it is disibuted on multiple processors. In other words,
!           each band of columns, called KKR_MatrixBand, is allocated 
!           on a processor as follows:
!               allocate ( KKR_MatrixBand(KKRMatrixSizeCant*BandSizeCant) )
!           Physically, KKR_MatrixBand is a matrix (2-D array) with 
!           row size = KKRMatrixSizeCant:
!
!                         [ ( 1, 1 )                 ( 1, 2 )               ...      ( 1, BandSizeCant )          ]
!                         | ( 2, 1 )                 ( 2, 2 )               ...      ( 2, BandSizeCant )          |
!                         |    .                        .                   ...         .                         |
!        KKR_MatrixBand = |    .                        .                   ...         .                         |
!                         |    .                        .                   ...         .                         |
!                         [ (KKRMatrixSizeCant,1)   (KKRMatrixSizeCant,2)   ...  (KKRMatrixSizeCant,BandSizeCant) ]
!
!           Of course, for thye number of processors = 1, BandSizeCant = KKRMatrixSizeCant
!        =============================================================
!
!        Leo: The following part needs to be worked on...
!        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 

!        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        End of the change.
!        =============================================================
      enddo
   enddo
!
   end subroutine calBandStructure
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupSCinGlobalFrame(getSingleScatteringMatrix)
!  ===================================================================
   use WriteMatrixModule,  only : writeMatrix
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SpinRotationModule, only : rotateLtoG
!
   implicit none
!
   integer (kind=IntKind) :: t0size, kkri_ns, i, kmax_kkr
!
   complex (kind=CmplxKind), pointer :: sm1(:,:), sm2(:,:)
   complex (kind=CmplxKind), pointer :: cm1(:,:), cm2(:,:)
   complex (kind=CmplxKind), pointer :: pm(:), gmat(:,:)
!
   interface
      function getSingleScatteringMatrix(smt,site,atom,dsize) result(sm)
         use KindParamModule, only : IntKind, CmplxKind
         character (len=*), intent(in) :: smt
         integer (kind=IntKind), intent(in) :: site, atom
         integer (kind=IntKind), intent(out), optional :: dsize
         complex (kind=CmplxKind), pointer :: sm(:,:)
      end function getSingleScatteringMatrix
   end interface
!
   cosine_g => stcm_g   ! Use stcm_g as the space for storing the cosine matrix
!
   if (isRelativistic) then !xianglin in between
      do i = 1, LocalNumAtoms
!        ================================================================
!        Obtain the Jinv-matrix, Sine-Matrix, and t-matrix in Global frame
!        ================================================================
         kmax_kkr = MatrixBand(i)%kmax_kkr
         kkri_ns =  kmax_kkr*nSpinCant
         t0size = kmax_kkr*kmax_kkr*nSpinCant*nSpinCant
         sm1 => getSingleScatteringMatrix('Sine-Matrix',i)
         cm1 => getSingleScatteringMatrix('Cosine-Matrix',i)
!        ----------------------------------------------------------------
         call zcopy( t0size, sm1, 1, sine_g(1,i), 1 )     !save s into sine_g
         call zcopy( t0size, cm1, 1, cosine_g(1,i), 1 )   !save jinv(Method 0) or s^t(Method 1) into jinv_g 
!        ----------------------------------------------------------------
      enddo
   else
      do i = 1, LocalNumAtoms
!        =============================================================
!        Obtain the Cosine-matrix, and Sine-Matrix t-matrix in Global frame
!        =============================================================
         kmax_kkr = MatrixBand(i)%kmax_kkr
         t0size = kmax_kkr*kmax_kkr
         if ( nSpinCant == 2 ) then
            kkri_ns =  kmax_kkr*nSpinCant
!           =============================================================
!           calculate sine_g and cosine_g in global frame of reference.
!           =============================================================
            sm1 => getSingleScatteringMatrix('Sine-Matrix',1,i)
            sm2 => getSingleScatteringMatrix('Sine-Matrix',2,i)
            pm => sine_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, sm1, sm2, gmat)
!           -------------------------------------------------------------
            cm1 => getSingleScatteringMatrix('Cosine-Matrix',1,i)
            cm2 => getSingleScatteringMatrix('Cosine-Matrix',2,i)
            pm => cosine_g(:,i)
            gmat => aliasArray2_c(pm,kkri_ns,kkri_ns)
!           -------------------------------------------------------------
            call rotateLtoG(i, kmax_kkr, kmax_kkr, cm1, cm2, gmat)
!           -------------------------------------------------------------
         else
            sm1 => getSingleScatteringMatrix('Sine-Matrix',1,i)
            cm1 => getSingleScatteringMatrix('Cosine-Matrix',1,i)
!           -------------------------------------------------------------
            call zcopy( t0size, sm1, 1, sine_g(1,i), 1 )
            call zcopy( t0size, cm1, 1, cosine_g(1,i), 1 )
!           -------------------------------------------------------------
         endif
      enddo
   endif
!
   nullify(sm1,sm2,cm1,cm2,pm,gmat)
!
   end subroutine setupSCinGlobalFrame
!  ===================================================================
