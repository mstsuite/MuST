module CurrentMatrixModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initCurrentMatrixModule, &
          endCurrentMatrixModule,  &
          calculateClebschGordanCoefficient, &
          calfx, &
          callsize, &
          kdelta, &
          populateClebschGordanTable, & 
          getClebschGordanCoefficient, &
          calSFsum, &
          calJxFromPhiLrIntegral, &
          calCurrentMatrix, &
          RadialIntegral, &
          calJyzFromJx, &
          calJtildeCPA, &
          getJMatrix
!

private
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: n_spin_cant
   integer (kind=IntKind) :: master_size

   integer (kind=IntKind), allocatable :: print_instruction(:)
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
   integer (kind=IntKind), allocatable :: kmax_kkr(:)
   integer (kind=IntKind), allocatable :: kmax_phi(:)
   integer (kind=Intkind), allocatable :: lofk(:), mofk(:), jofk(:), m1m(:)
   integer (kind=Intkind), allocatable :: lofj(:), mofj(:)
   integer (kind=IntKind), allocatable :: num_species(:)
   integer (kind=IntKind) :: lmax_phi_max, kmax_kkr_max, lmax_green_max, &
                             lmax_sigma, lmax_sigma_2, lmax_cg
   integer (kind=IntKind) :: kmax_phi_max, kmax_green_max, iend_max, &
                             kmax_sigma, kmax_sigma_2, kmax_max, kmax_cg
   integer (kind=IntKind) :: NumSpecies

   type CGCoeff
      integer (kind=IntKind) :: lsize
      integer (kind=IntKind), pointer :: l3_list(:)
      real (kind=RealKind), pointer :: cgl3(:)
   end type CGCoeff

   type (CGCoeff), allocatable :: cg(:,:)


   real (kind=RealKind), allocatable, target :: cgspace(:,:,:)
   integer (kind=IntKind), allocatable, target :: l3space(:,:,:)
   complex (kind=CmplxKind), allocatable, target :: gspacep(:)


   complex (kind=CmplxKind), allocatable :: iden(:,:)
   complex (kind=CmplxKind) :: ep1(3), em1(3), e0(3)

!  Current data stored here
!  --------------------------------------------------------------------
   complex (kind=CmplxKind), allocatable, target :: jspace(:,:,:,:,:,:),  &
       jspace2(:,:,:,:,:,:), jspace3(:,:,:,:,:,:), jspace4(:,:,:,:,:,:)

   complex (kind=CmplxKind), allocatable, target :: jtspace(:,:,:,:,:,:), &
       jtspace2(:,:,:,:,:,:), jtspace3(:,:,:,:,:,:), jtspace4(:,:,:,:,:,:)
!  ---------------------------------------------------------------------
   complex (kind=CmplxKind), pointer :: gaunt(:,:,:)
   complex (kind=CmplxKind), allocatable :: store_space(:)

   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, kGID
   logical :: Initialized = .false.
! 
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initCurrentMatrixModule(energy, num_atoms, lmaxkkr, lmaxphi, lmaxgreen,  &
                              pola, cant, rel, istop, iprint, mode)
!  ===================================================================
   use RadialGridModule, only : getNumRmesh, getMaxNumRmesh
   use AtomModule, only : getLocalNumSpecies
   use PolyhedraModule, only : getNumPolyhedra
   use WriteMatrixModule, only : writeMatrix
   use ScfDataModule, only : isKKR, isLSMS, isKKRCPA, isKKRCPASRO
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
   use GauntFactorsModule, only : getK3, getNumK3, getGauntFactor

   integer (kind=IntKind), intent(in) :: num_atoms
   integer (kind=IntKind), intent(in) :: lmaxkkr(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxphi(num_atoms)
   integer (kind=IntKind), intent(in) :: lmaxgreen(num_atoms)
   integer (kind=IntKind), intent(in) :: pola
   integer (kind=IntKind), intent(in) :: cant
   integer (kind=IntKind), intent(in) :: rel
   character (len=*), intent(in) :: istop
   integer (kind=IntKind), intent(in) :: iprint(num_atoms)
   integer (kind=IntKind), intent(in) :: mode
   
   complex (kind=CmplxKind), intent(in) :: energy

   integer (kind=LongIntKind) :: wspace_size, gspace_size
   integer (kind=IntKind) :: i, lmax_max, jmax, iend, kmax, NumSpecies, NumPolyhedra
   integer (kind=IntKind) :: lmax, kl, jl, m, n, l, j, jsize
   integer (kind=IntKind) :: klp1, klp2, i3, klg
   integer (kind=IntKind), pointer :: nj3(:,:), kj3(:,:,:)

   real (kind=RealKind), pointer :: cgnt(:,:,:)
   
   complex (kind=CmplxKind) :: tmp

   LocalNumAtoms = num_atoms
   n_spin_pola = pola
   n_spin_cant = cant
   NumPolyhedra = getNumPolyhedra()

   allocate( print_instruction(LocalNumAtoms) ) 
   allocate( lmax_kkr(LocalNumAtoms) )
   allocate( lmax_phi(LocalNumAtoms) )
   allocate( kmax_kkr(LocalNumAtoms) )
   allocate( kmax_phi(LocalNumAtoms) )
   allocate( num_species(LocalNumAtoms) )
!
   kmax_kkr_max = 1
   kmax_phi_max = 1
   kmax_green_max = 1
   lmax_phi_max = 0
   lmax_green_max = 0
   lmax_max = 0
   NumSpecies = 0

   do i=1, LocalNumAtoms
      lmax_kkr(i) = lmaxkkr(i)
      kmax_kkr(i) = (lmaxkkr(i)+1)**2
      lmax_phi(i) = lmaxphi(i)
      kmax_phi(i) = (lmaxphi(i)+1)**2
      print_instruction(i) = iprint(i)
      num_species(i) = getLocalNumSpecies(i)
      kmax_kkr_max = max(kmax_kkr_max, (lmaxkkr(i)+1)**2)
      kmax_phi_max = max(kmax_phi_max, (lmaxphi(i)+1)**2)
      kmax_green_max = max(kmax_green_max, (lmaxgreen(i)+1)**2)
      lmax_phi_max = max(lmax_phi_max, lmaxphi(i))
      lmax_green_max = max(lmax_green_max, lmaxgreen(i))
      lmax_max = max(lmax_max, lmaxgreen(i), lmaxkkr(i), lmaxphi(i), &
                     lmax_sigma, lmax_sigma_2)
      NumSpecies = max(NumSpecies, num_species(i))
   enddo
   
   lmax_sigma = lmax_max + 2
   lmax_sigma_2 = lmax_max + 2
   kmax_sigma = (lmax_sigma + 1)**2
   kmax_sigma_2 = (lmax_sigma_2 + 1)**2

   lmax_cg = lmax_sigma
   kmax_cg = (lmax_cg + 1)**2
   master_size = kmax_kkr_max
   
   allocate(gspacep(kmax_cg*kmax_cg*kmax_cg))

!  -------------------------------------------------------------------
   call initGauntFactors(lmax_max,istop,-1)
!  -------------------------------------------------------------------
   iend_max = getMaxNumRmesh()

!  -------------------------------------------------------------------
!  Initialize the arrays which store conductivity data
!  -------------------------------------------------------------------
   jsize = master_size

   allocate(lofk(kmax_cg), mofk(kmax_cg),  &
         jofk(kmax_cg), m1m(-lmax_cg:lmax_cg))
   allocate(lofj((lmax_cg+1)*(lmax_cg+2)/2), mofj((lmax_cg+1)*(lmax_cg+2)/2))
   allocate(cg(kmax_cg, kmax_cg))
   allocate(cgspace(kmax_cg, kmax_cg, 2*lmax_cg + 2))
   allocate(l3space(kmax_cg, kmax_cg, 2*lmax_cg + 2))
   allocate(jspace(jsize, jsize, &
             LocalNumAtoms, NumSpecies, n_spin_pola, 3), &
   jspace2(jsize, jsize, LocalNumAtoms, NumSpecies, n_spin_pola, 3), &
   jspace3(jsize, jsize, LocalNumAtoms, NumSpecies, n_spin_pola, 3), &
   jspace4(jsize, jsize, LocalNumAtoms, NumSpecies, n_spin_pola, 3)) 
   
   if (mode == 3 .or. mode == 4) then
     allocate(jtspace(jsize, jsize, &
             LocalNumAtoms, NumSpecies, n_spin_pola, 3), &
     jtspace2(jsize, jsize, LocalNumAtoms, NumSpecies, n_spin_pola, 3), &
     jtspace3(jsize, jsize, LocalNumAtoms, NumSpecies, n_spin_pola, 3), &
     jtspace4(jsize, jsize, LocalNumAtoms, NumSpecies, n_spin_pola, 3))
     jtspace = CZERO; jtspace2 = CZERO; jtspace3 = CZERO; jtspace4 = CZERO
   endif

   allocate(iden(jsize, jsize))

   cgspace = ZERO
   l3space = 0
   jspace = CZERO; jspace2 = CZERO; jspace3 = CZERO; jspace4 = CZERO
   iden = CZERO

   do i = 1, jsize
      iden(i, i) = CONE
   enddo

!  ===================================================================
!  calculate the factors: lofk, mofk, jofk, lofj, and mofj............
!  ===================================================================
   kl=0; jl = 0
   do l=0,lmax_cg
      n=(l+1)*(l+2)/2-l
      do m=-l,l
         kl=kl+1
         lofk(kl)=l
         mofk(kl)=m
         jofk(kl)=n+abs(m)
         if (m >= 0) then
            jl = jl + 1
            lofj(jl) = l
            mofj(jl) = m
         endif
      enddo
   enddo
!
!  ===================================================================
!  calculate the factor (-1)**m and store in m1m(-lmax:lmax)..........
!  ===================================================================
   m1m(0)=1
   do m=1,lmax_cg
      m1m(m)=-m1m(m-1)
   enddo
   do m=-1,-lmax_max,-1
      m1m(m)=-m1m(m+1)
   enddo

   nj3 => getNumK3()
   kj3 => getK3()
   cgnt => getGauntFactor()
   gaunt => aliasArray3_c(gspacep,kmax_cg,kmax_cg,kmax_cg)
   gaunt = CZERO

   do klp2 = 1, kmax_cg
     do klg = 1, kmax_cg
       do klp1 = 1, kmax_cg
          do i3 = 1, nj3(klp1,klg)
            if (kj3(i3,klp1,klg) == klp2) then
               gaunt(klp1,klg,klp2) = cgnt(i3,klp1,klg)
            endif
          enddo
       enddo
     enddo
   enddo

   nullify(nj3, kj3, cgnt)
!
!  -------------------------------------------------------------------
   call endGauntFactors()
!  -------------------------------------------------------------------

!  l = 1 spherical basis vectors expressed in a Cartesian basis

   ep1(1) = -CONE/sqrt(2.0)
   ep1(2) = -SQRTm1/sqrt(2.0)
   ep1(3) = CZERO
  
   em1(1) = CONE/sqrt(2.0)
   em1(2) = -SQRTm1/sqrt(2.0)
   em1(3) = CZERO

   e0(1) = CZERO
   e0(2) = CZERO
   e0(3) = CONE
    
!  --------------------------------------------------------------------
   call populateClebschGordanTable()
!  --------------------------------------------------------------------

   Initialized = .true.

   end subroutine initCurrentMatrixModule
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endCurrentMatrixModule()
!  ===================================================================

   integer (kind=IntKind) :: L_1, L_2

   deallocate(jspace, jspace2, jspace3, jspace4)
   deallocate(jtspace, jtspace2, jtspace3, jtspace4)

   do L_1 = 1, kmax_cg
     do L_2 = 1, kmax_cg
       nullify(cg(L_1, L_2)%cgl3, cg(L_1, L_2)%l3_list)
     enddo
   enddo
   
   deallocate(cg)
   deallocate(print_instruction, lmax_kkr, lmax_phi, kmax_kkr, kmax_phi, &
   lofk, mofk, jofk, m1m, lofj, mofj, num_species, iden)  
   deallocate(cgspace, l3space, gspacep)
   nullify(gaunt)

   end subroutine endCurrentMatrixModule
!  ===================================================================
   
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calculateClebschGordanCoefficient(l_1, m_1, l_2,  m_2)  result(cg)
!  ===================================================================
!
!  The function calculates  <l_1m_1l_2m_2|l(m_1+m_2)>
!  for |l_1 - l_2|  <= l <= (l_1 + l_2)
 
   integer (kind=IntKind), intent(in) :: l_1, l_2, m_1, m_2
   integer (kind=IntKind) :: i, temp, lsize
   real (kind=RealKind) :: coeff, log_coeff, first_term, &
                             second_term, third_term
   real (kind=RealKind) :: A_ml, A_pl, A_ol
   real (kind=RealKind) :: numerator, denominator
   real (kind=RealKind), allocatable :: cg(:)
   integer (kind=IntKind), allocatable :: l_3(:)
   
   lsize = callsize(l_1, m_1, l_2, m_2)
   allocate(cg(lsize), l_3(lsize))

   l_3(lsize) = l_1 + l_2 + 1
   do i = lsize-1, 1,-1
      l_3(i) = l_3(i + 1) - 1
   enddo
   
   cg(lsize) = ZERO
   
   coeff = (calProductSeries(2*l_1, l_1))/calProductSeries(2*l_1 + 2*l_2, 2*l_2)
   coeff = coeff*gamma(l_1 + 1.0)
   
   if (mod(l_1 - m_1, 2) == 0) then
      coeff = coeff*(calProductSeries(l_1 + l_2 + m_1 + m_2, l_1 + m_1)/&
                 calProductSeries(l_1 - m_1, (l_1 - m_1)/2))
      coeff = coeff*(1.0/gamma((l_1 - m_1)/2 + 1.0))
   else
      coeff = coeff*(calProductSeries(l_1 + l_2 + m_1 + m_2, l_1 + m_1)/&
                  calProductSeries(l_1 - m_1, (l_1 - m_1 + 1)/2))
      coeff = coeff*(1.0/gamma((l_1 - m_1 + 1)/2 + 1.0))
   endif

   if (mod(l_2 + m_2, 2) == 0) then 
      coeff = coeff*(calProductSeries(l_1 + l_2 - m_1 - m_2, l_2 - m_2)/&
                 calProductSeries(l_2 + m_2, (l_2 + m_2)/2))
      coeff = coeff*(1.0/gamma((l_2 + m_2)/2 + 1.0)) 
   else 
      coeff = coeff*(calProductSeries(l_1 + l_2 - m_1 - m_2, l_2 - m_2)/&
                  calProductSeries(l_2 +  m_2, (l_2 + m_2 + 1)/2))
      coeff = coeff*(1.0/gamma((l_2 + m_2 + 1)/2 + 1.0))
   endif
   
   cg(lsize - 1) = sqrt(coeff)
   
   do i = lsize - 1,2, -1
      temp = l_3(i)
      A_ol = m_1 - m_2 + (m_1 + m_2)*(l_2*(l_2 + 1.0) - l_1*(l_1 + 1.0))/(temp*(temp + 1.0))
      A_pl = calfx(temp + 1, l_1, m_1, l_2, m_2)
      A_ml = calfx(temp, l_1, m_1, l_2, m_2)
      cg(i - 1) = (A_ol/A_ml)*cg(i) - (A_pl/A_ml)*cg(i + 1)
   enddo

   end function calculateClebschGordanCoefficient
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calProductSeries(start_val, end_val) result(prod)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: start_val, end_val
   integer (kind=IntKind) :: i
   real (kind=RealKind) :: prod
   
   prod = 1.0
   do i = start_val, end_val + 1, -1
      prod = prod*i
   enddo

   end function calProductSeries
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calfx(x, l_1, m_1, l_2, m_2)  result(A)
!  ===================================================================
   implicit none

   integer (kind=IntKind), intent(in) :: x, l_1, m_1, l_2, m_2
   integer (kind=IntKind) :: m
   real (kind=RealKind) :: pre_A
   real (kind=RealKind) :: A

   m = m_1 + m_2
   pre_A = (x**2 - m**2)*((l_1 + l_2 + 1.0)**2.0 - x**2.0)*&
       ((x**2 - (l_1 - l_2)**2.0)/(x**2*(4*x**2.0 - 1.0)))
   
   A = sqrt(pre_A)
   
   end function calfx
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function callsize(l_1, m_1, l_2, m_2) result(lsize)
!  ===================================================================
   
   integer (kind=IntKind), intent(in) :: l_1, m_1, l_2, m_2
   integer (kind=IntKind) :: lsize

   if (abs(m_1 + m_2) <= abs(l_1 - l_2)) then
     lsize = 2*min(l_1,l_2) + 2 
   else
     lsize = l_1 + l_2 - abs(m_1 + m_2) + 2 
   endif

   end function callsize
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function kdelta(i, j) result(kdval)
!  ===================================================================

   integer (kind=IntKind), intent(in) :: i, j
   integer (kind=IntKind) :: kdval

   if (i == j) then
     kdval = 1
   else
     kdval = 0
   endif
   
   end function kdelta
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine populateClebschGordanTable()
!  ===================================================================

   integer (kind=IntKind) :: l_1, m_1, l_2, m_2, lsize, k_1, k_2
   integer (kind=IntKind) :: i
   lsize = callsize(l_1, m_1, l_2, m_2)
  
   do k_2 = 1, kmax_cg
     do k_1 = 1, kmax_cg
       l_1 = lofk(k_1); m_1 = mofk(k_1)
       l_2 = lofk(k_2); m_2 = mofk(k_2)
       lsize = callsize(l_1, m_1, l_2, m_2)
       cg(k_1, k_2)%lsize = lsize
       cgspace(k_1, k_2, 1:lsize) =  &
          calculateClebschGordanCoefficient(l_1, m_1, l_2, m_2)
       cg(k_1, k_2)%cgl3 => cgspace(k_1, k_2, 1:lsize)
       l3space(k_1, k_2, lsize) = l_1 + l_2 + 1
       do i = lsize - 1, 1, -1
          l3space(k_1, k_2, i) = l3space(k_1, k_2, i+1) - 1
       enddo
       cg(k_1, k_2)%l3_list => l3space(k_1, k_2, 1:lsize)
     enddo
   enddo
  
   end subroutine populateClebschGordanTable
!  ===================================================================
  
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getClebschGordanCoefficient(l_1, m_1, l_2, m_2, l_3, m_3) result(cg_coeff)
!  ===================================================================
!
!  returns the pre-calculated coefficient <l_1m_1l_2m_2|l_3m_3>
!
   integer (kind=IntKind), intent(in) :: l_1, m_1, l_2, m_2, l_3, m_3
   integer (kind=IntKind) :: k_1, k_2, lsize
   real (kind=RealKind) :: cg_coeff
   
   cg_coeff = CZERO

   if (abs(m_1) > l_1 .or. abs(m_2) > l_2 .or. abs(m_3) > l_3) then
      call ErrorHandler("getClebschGordanCoefficient", &
              "Illegal value of m_1 and/or m_2 and/or m_3")
   else if (m_1 + m_2 /= m_3) then
      cg_coeff = ZERO
   else if (l_3 > l_1 + l_2 + 1 .or. l_3 < abs(l_1 - l_2)) then
      cg_coeff = ZERO
   else
      k_1 = l_1**2 + l_1 + m_1 + 1
      k_2 = l_2**2 + l_2 + m_2 + 1
      lsize = cg(k_1, k_2)%lsize
      do i = 1, lsize
        if (l_3 == cg(k_1,k_2)%l3_list(i)) then
           cg_coeff = cg(k_1,k_2)%cgl3(i)
        endif
      enddo
   endif
   
   end function getClebschGordanCoefficient
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSFsum(n, r, iend, gL, gK2, dir, c, qterm) result(SFqsum)
!  ===================================================================
 
   use StepFunctionModule, only : interpolateStepFunction
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, iend, gL, gK2
   integer (kind=IntKind) ,intent(in) :: dir, c
   real (kind=RealKind), intent(in) :: r(iend)
   complex (kind=CmplxKind), intent(in) :: qterm(iend, kmax_sigma_2, kmax_sigma_2)

   integer (kind=IntKind) :: q, k_2, mk_2, gK2q, i
   real (kind=RealKind) :: cg_coeff
   complex (kind=CmplxKind) :: alpha
   complex (kind=CmplxKind) :: SFqsum(iend)

   k_2 = lofk(gK2)
   mk_2 = mofk(gK2)
   SFqsum = CZERO
   alpha = CZERO

   if (c == 1) then
      do q = -1, 1
        if (k_2 == lofk(kmax_cg)) then
           SFqsum = SFqsum + CZERO
        else
          gK2q = (k_2 + 1)**2 + (k_2 + 1) + (mk_2 + q) + 1
          cg_coeff = getClebschGordanCoefficient(k_2 + 1, mk_2 + q, 1, -q, k_2, mk_2)
          if (q == -1) then
            alpha = sqrt((k_2 + 1.0)/(2*k_2 + 1.0))*ep1(dir)*cg_coeff
!           ------------------------------------------------------
            call zaxpy(iend, alpha, qterm(:,gK2q,gL), 1, SFqsum, 1)
!           ------------------------------------------------------
          else if (q == 0) then
            alpha = sqrt((k_2 + 1.0)/(2*k_2 + 1.0))*e0(dir)*cg_coeff
!           ------------------------------------------------------
            call zaxpy(iend, alpha, qterm(:,gK2q,gL), 1, SFqsum, 1)
!           ------------------------------------------------------
          else
            alpha = sqrt((k_2 + 1.0)/(2*k_2 + 1.0))*em1(dir)*cg_coeff
!           ------------------------------------------------------
            call zaxpy(iend, alpha, qterm(:,gK2q,gL), 1, SFqsum, 1)
!           ------------------------------------------------------
          endif
        endif
      enddo
   else if (c == -1) then
      do q = -1, 1
         if (k_2 == mk_2 .and. (q == 1 .or. q == 0)) then
            SFqsum = SFqsum + CZERO
         else if (k_2 == -mk_2 .and. (q == -1 .or. q == 0)) then
            SFqsum = SFqsum + CZERO
         else if (k_2 == mk_2 + 1 .and. q == 1) then
            SFqsum = SFqsum + CZERO
         else if (k_2 == abs(mk_2 - 1) .and. q == -1) then
            SFqsum = SFqsum + CZERO
         else
            gK2q = (k_2 - 1)**2 + (k_2 - 1) + (mk_2 + q) + 1
            cg_coeff = getClebschGordanCoefficient(k_2 - 1, mk_2 + q, 1, -q, k_2, mk_2)
            if (q == -1) then
              alpha = sqrt(k_2/(2*k_2 + 1.0))*ep1(dir)*cg_coeff
!             ------------------------------------------------------
              call zaxpy(iend, alpha, qterm(:,gK2q,gL), 1, SFqsum, 1)
!             ------------------------------------------------------
            else if (q == 0) then
              alpha = sqrt(k_2/(2*k_2 + 1.0))*e0(dir)*cg_coeff
!             ------------------------------------------------------
              call zaxpy(iend, alpha, qterm(:,gK2q,gL), 1, SFqsum, 1)
!             ------------------------------------------------------
            else
              alpha = sqrt(k_2/(2*k_2 + 1.0))*em1(dir)*cg_coeff
!             ------------------------------------------------------
              call zaxpy(iend, alpha, qterm(:,gK2q,gL), 1, SFqsum, 1)
!             ------------------------------------------------------
            endif
         endif
      enddo
   else
      call ErrorHandler("calVSHcoeff", "Incorrect value for choice variable", c)
   endif

   end function calSFsum
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calCGCoeff(gK1, gK2, dir, choice) result(ri_prefactor)
!  ===================================================================
   
   integer (kind=IntKind), intent(in) :: gK1, gK2, dir, choice
   integer (kind=IntKind) :: k_1, mk_1, k_2, mk_2, q
   real (kind=RealKind) :: cg_coeff
   complex (kind=CmplxKind) :: e_vec, ri_prefactor

   k_1 = lofk(gK1); mk_1 = mofk(gK1)
   k_2 = lofk(gK2); mk_2 = mofk(gK2)

   ri_prefactor = CZERO
   if (choice == 0) then
      do q = -1, 1
         if (k_2 == mk_2 .and. (q == 0 .or. q == 1)) then
           cg_coeff = 0.0
         else if (k_2 == -mk_2 .and. (q == 0 .or. q == -1)) then
           cg_coeff = 0.0
         else if (k_2 == mk_2 + 1 .and. q == 1) then
           cg_coeff = 0.0
         else if (k_2 == -mk_2 + 1 .and. q == -1) then
           cg_coeff = 0.0
         else
           cg_coeff = getClebschGordanCoefficient(k_2 - 1, mk_2 + q, &
                    1, -q, k_2, mk_2)
           if (q == -1) then
             e_vec = ep1(dir)
           else if (q == 0) then
             e_vec = e0(dir)
           else
             e_vec = em1(dir)
           endif
           ri_prefactor = ri_prefactor + &
             sqrt(k_2/(2*k_2 + 1.0))*e_vec*cg_coeff*kdelta(k_1, k_2 - 1)*kdelta(mk_1, mk_2 + q)
         endif
      enddo
   else if (choice == 1) then
      do q = -1, 1
         if (k_2 == lofk(kmax_cg)) then
            cg_coeff = 0.0
         else
            cg_coeff = getClebschGordanCoefficient(k_2 + 1, mk_2 + q, &
                    1, -q, k_2, mk_2)
            if (q == -1) then
              e_vec = ep1(dir)
            else if (q == 0) then
              e_vec = e0(dir)
            else 
              e_vec = em1(dir)
            endif
         endif
         ri_prefactor = ri_prefactor + &
           sqrt((k_2 + 1.0)/(2*k_2 + 1.0))*(e_vec)*cg_coeff*kdelta(k_1,k_2+1)*kdelta(mk_1,mk_2+q)
      enddo
   endif 

   end function calCGCoeff
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function RadialIntegral(iend, n, ic, is, gL1, gL2, gL,  &
                             gK1, gK2, dir, choice, qterm, mt) result(rad_val)
!  ===================================================================
!  We are solving the integral 
!  
!  if choice == 0
!  \int_r r^2 chi_L_5(r) R^{m}_{K_1L_4} [d_r R^{m}_{K_1L_1}(r)  \
!                          + (k_2 + 1)/r R^{m}_{K_2L_1}(r)]
!  
!  if choice == 1
!  \int_r r^2 chi_L_5(r) R^{m}_{K_1L_4} [d_r R^{m}_{K_2L_1}(r)  \
!                          - k_2/r R^{m}_{K_2L_1}(r)]
!  -------------------------------------------------------------------
   use StepFunctionModule, only : getStepFunction
   use RadialGridModule, only : getRmesh, getRadialGridRadius
   use SSSolverModule, only : getRegSolution, getRegSolutionDerivative, &
                           getSolutionRmeshSize, getSineMatrix
   use IntegrationModule, only : calIntegration
   use WriteMatrixModule, only : writeMatrix   
   use MatrixInverseModule, only : MtxInv_LU

   integer (kind=IntKind), intent(in) :: iend, n, ic, is
   integer (kind=IntKind), intent(in) :: gK1, gK2, gL1, gL2, gL
   integer (kind=IntKind), intent(in) :: dir, choice, mt
   complex (kind=CmplxKind), intent(in) :: qterm(iend, kmax_sigma_2, kmax_sigma_2)

   integer (kind=IntKind) :: k_2
   integer (kind=IntKind) :: i, dsize
   real (kind=RealKind) :: rmt
   real (kind=RealKind), pointer :: radial_grid(:)
   
   complex (kind=CmplxKind) :: rad_val, mt_coeff
   complex (kind=CmplxKind), pointer :: phi(:,:,:), dphi(:,:,:)
   complex (kind=CmplxKind), allocatable :: sfvsh(:), dsfvsh(:), term_2(:)
   complex (kind=CmplxKind), allocatable :: phi_k1l(:), &
                                    dphi_k2lp(:), phi_k2lp(:), chi_k1l(:)
   complex (kind=CmplxKind), allocatable :: product_1(:)
   
   allocate(phi_k1l(iend), sfvsh(iend), dsfvsh(iend), dphi_k2lp(iend),  &
       term_2(iend), phi_k2lp(iend), chi_k1l(iend), product_1(iend))

   product_1 = CZERO
   rad_val = CZERO
   sfvsh = CZERO
   dsfvsh = CZERO
   phi_k1l = CZERO
   dphi_k2lp = CZERO
   phi_k2lp = CZERO
   chi_k1l = CZERO
   term_2 = CZERO
   mt_coeff = CZERO

   k_2 = lofk(gK2)
   radial_grid => getRmesh(n)
   phi => getRegSolution(spin=is,site=n,atom=ic)

   dphi => getRegSolutionDerivative(spin=is,site=n,atom=ic)


   phi_k1l = phi(1:iend,gK1,gL1)

   dphi_k2lp = dphi(1:iend,gK2, gL2)
   phi_k2lp = phi(1:iend,gK2,gL2)

   if (mt == 0) then
     chi_k1l = qterm(:,gK1,gL)
   endif
   
   rmt = getRadialGridRadius(1, MT=.true.)
   if (mt == 1) then
     mt_coeff = calCGCoeff(gK1, gK2, dir, choice)
     if (choice == 0) then
        do i = 1, iend
          term_2(i) = dphi_k2lp(i) + ((k_2 + 1.0)/radial_grid(i))*phi_k2lp(i)
        enddo
     else if (choice == 1) then
        do i = 1, iend
          term_2(i) = dphi_k2lp(i) - (k_2/radial_grid(i))*phi_k2lp(i)
        enddo
     endif
     do i = 1, iend
        product_1(i) = mt_coeff*phi_k1l(i)*term_2(i)
     enddo
   
   else
     if (choice == 0) then
       sfvsh = calSFsum(n, radial_grid, iend, gL, gK2, dir, -1, qterm)
       dsfvsh = calculateStepFunctionDerivative(iend, radial_grid, sfvsh) 
       do i = 1, iend
          term_2(i) = dphi_k2lp(i)*sfvsh(i) + phi_k2lp(i)*dsfvsh(i) &
            + sfvsh(i)*(k_2 + 1.0)*(phi_k2lp(i)/radial_grid(i))
       enddo
     else if (choice == 1) then
       sfvsh = calSFsum(n, radial_grid, iend, gL, gK2, dir, 1, qterm)
       dsfvsh = calculateStepFunctionDerivative(iend, radial_grid, sfvsh)
       do i = 1, iend
          term_2(i) = dphi_k2lp(i)*sfvsh(i) + phi_k2lp(i)*dsfvsh(i) &
            - sfvsh(i)*k_2*(phi_k2lp(i)/radial_grid(i))
       enddo
     else
       call ErrorHandler('RadialIntegral', 'Incorrect choice of integral', choice)
     endif
    
     do i = 1, iend
       product_1(i) = conjg(chi_k1l(i))*phi_k1l(i)*term_2(i)
     enddo
   endif

!  ---------------------------------------------------------------------
   call calIntegration(iend, radial_grid, product_1, rad_val, 0)
!  ---------------------------------------------------------------------

   if (mt == 1) then
     rad_val = rad_val - 0.5*mt_coeff*conjg(phi_k1l(iend))*phi_k2lp(iend)
   endif
   
   end function RadialIntegral
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSurfaceTerm(iend, n, ic, is, gL1, gL2, gL, gL6, &
                      gK1, gK2, dir, qterm2, qterm1) result(s_term)
!  ===================================================================

   use SSSolverModule, only : getRegSolution
   use RadialGridModule, only : getRmesh
   use IntegrationModule, only : calIntegration

   integer (kind=IntKind), intent(in) :: iend, n, ic, is, gL, &
                                  gL1, gL2, gK1, gK2, dir, gL6
   complex (kind=CmplxKind), intent(in) :: qterm2(iend,kmax_sigma_2,kmax_sigma_2), & 
                                    qterm1(iend, jofk(kmax_sigma))   
   integer (kind=IntKind) :: i, q, l_6, m_6, gL6p, gL6m, prefactor
   real (kind=RealKind) :: cgl6m, cgl6p
   complex (kind=CmplxKind) :: gauntl6m, gauntl6p, s_term
   real (kind=RealKind), pointer :: radial_grid(:)
   complex (kind=CmplxKind) :: s_val
   complex (kind=CmplxKind), pointer :: phi(:,:,:)
   complex (kind=CmplxKind), allocatable :: phi_k1l1(:), phi_k2l2(:), &
                            chi_k1l(:), chi_l6(:), product_2(:)
   
   s_val = CZERO
   s_term = CZERO
   l_6 = lofk(gL6)
   m_6 = mofk(gL6)
   do q = -1, 1
     if (l_6 == m_6 .and. (q == 1 .or. q == 0)) then
        s_val = s_val + CZERO
     else if (l_6 == -m_6 .and. (q == -1 .or. q == 0)) then
        s_val = s_val + CZERO
     else if (l_6 == m_6 + 1 .and. q == 1) then
        s_val = s_val + CZERO
     else if (l_6 == abs(m_6 - 1) .and. q == -1) then
        s_val = s_val + CZERO
     else
        gL6m = (l_6 - 1)**2 + (l_6 - 1) + (m_6 + q) + 1
        cgl6m = getClebschGordanCoefficient(l_6 - 1, m_6 + q, 1, -q, l_6, m_6)
        gauntl6m = gaunt(gK2, gL, gL6m)
        if (q == -1) then 
          s_val = s_val + ep1(dir)*(l_6 + 1)*sqrt(l_6/(2*l_6 + 1.0))*cgl6m*gauntl6m
        else if (q == 0) then
          s_val = s_val + e0(dir)*(l_6 + 1)*sqrt(l_6/(2*l_6 + 1.0))*cgl6m*gauntl6m
        else
          s_val = s_val + em1(dir)*(l_6 + 1)*sqrt(l_6/(2*l_6 + 1.0))*cgl6m*gauntl6m
        endif
     endif
       
     if (l_6 == lmax_sigma) then
        s_val = s_val + CZERO
     else
        gL6p = (l_6 + 1)**2 + (l_6 + 1) + (m_6 + q) + 1
        cgl6p = getClebschGordanCoefficient(l_6 + 1, m_6 + q, 1, -q, l_6, m_6)
        gauntl6p = gaunt(gK2, gL, gL6p)
        if (q == -1) then
          s_val = s_val + ep1(dir)*l_6*sqrt((l_6 + 1)/(2*l_6 + 1.0))*cgl6p*gauntl6p
        else if (q == 0) then
          s_val = s_val + e0(dir)*l_6*sqrt((l_6 + 1)/(2*l_6 + 1.0))*cgl6p*gauntl6p
        else
          s_val = s_val + em1(dir)*l_6*sqrt((l_6 + 1)/(2*l_6 + 1.0))*cgl6p*gauntl6p
        endif
     endif
   enddo

   allocate(phi_k1l1(iend), phi_k2l2(iend), chi_k1l(iend), chi_l6(iend), product_2(iend))

   if (mofk(gL6) < 0) then
     prefactor = (-1)**mofk(gL6)
   else
     prefactor = 1
   endif

   radial_grid => getRmesh(n)
   phi => getRegSolution(spin=is,site=n,atom=ic)
   phi_k1l1 = phi(:,gK1,gL1)
   phi_k2l2 = phi(:,gK2,gL2)
   chi_k1l = qterm2(:,gK1,gL)
   chi_l6 = prefactor*qterm1(:,jofk(gL6))

   do i = 1, iend
      product_2(i) = (s_val/radial_grid(i))*conjg(phi_k1l1(i))*conjg(chi_k1l(i))*phi_k2l2(i)*chi_l6(i)
   enddo

!  ------------------------------------------------------------------
   call calIntegration(iend, radial_grid, product_2, s_term, 0)
!  ------------------------------------------------------------------

   end function calSurfaceTerm
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calculateStepFunctionDerivative(iend, r, B) result(dB)
!  ==================================================================

   integer (kind=IntKind), intent(in) :: iend
   real (kind=RealKind), intent(in) :: r(iend)
   complex (kind=CmplxKind), intent(in) :: B(iend)

   integer (kind=IntKind) :: k
   complex (kind=CmplxKind) :: dB(iend)

   dB(1) = (B(1) - B(2))/(r(1) - r(2))
   do k = 2, iend - 1
      dB(k) = (B(k+1) - B(k-1))/(r(k+1) - r(k-1))
   enddo
   dB(iend) = (B(iend) - B(iend-1))/(r(iend) - r(iend-1))   

   end function calculateStepFunctionDerivative
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calJxFromPhiLrIntegral(n, ic, global_energy, iend, is, qterm1, qterm2, mt)
!  ===================================================================
   use MatrixInverseModule, only : MtxInv_LU
   use WriteMatrixModule, only : writeMatrix
   use SSSolverModule, only : getSineMatrix

   integer (kind=IntKind), intent(in) :: iend, is, mt
   complex (kind=CmplxKind), intent(in) :: global_energy
   complex (kind=CmplxKind), intent(in) :: qterm1(iend,jofk(kmax_sigma)), &
                            qterm2(iend,kmax_sigma_2,kmax_sigma_2)

   integer (kind=IntKind) :: n, ic, dir, calsize
   integer (kind=IntKind) :: ir, gK1, gK2, gL, gL1, gL2, gL6
   complex (kind=CmplxKind) :: kappa, surf
   complex (kind=CmplxKind), pointer :: sine_tmp(:,:), sine_transpose(:,:)
   complex (kind=CmplxKind), allocatable :: Dplus(:,:), Dminus(:,:), temp(:,:), &
                                          sine_mat(:,:), temp2(:,:)
   calsize = master_size
   allocate(sine_mat(calsize, calsize))
   allocate(sine_transpose(calsize, calsize))
   allocate(Dplus(calsize, calsize), Dminus(calsize, calsize), &
         temp(calsize, calsize), temp2(calsize, calsize))

   sine_mat = CZERO
   sine_transpose = CZERO

    
   sine_tmp => getSineMatrix(spin=is,site=n,atom=ic)
   sine_mat = sine_tmp(1:calsize, 1:calsize)
!  -------------------------------------------------------
   call MtxInv_LU(sine_mat, calsize)
!  -------------------------------------------------------
   call zgemm('T', 'n', calsize, calsize, calsize, CONE, sine_mat, &
        calsize, iden, calsize, CZERO, sine_transpose, calsize)
!  -------------------------------------------------------
   do dir = 1, 1
      Dplus = CZERO
      Dminus = CZERO
      temp = CZERO
      temp2 = CZERO
      if (mt == 0) then
         do gL1 = 1, calsize
           do gL2 = 1, calsize
             do gK1 = gL1, gL1
               do gK2 = gL2, gL2
                 do gL = 1, kmax_sigma_2
                   Dplus(gL1, gL2) = Dplus(gL1, gL2) + &
                     RadialIntegral(iend,n,ic,is,gL1,gL2,gL,gK1,gK2,1,0,qterm2,mt)  &
                    -RadialIntegral(iend,n,ic,is,gL1,gL2,gL,gK1,gK2,1,1,qterm2,mt)
                   surf = CZERO
                   do gL6 = 1, kmax_sigma
                     surf = surf + calSurfaceTerm(iend,n,ic,is,gL1,gL2,gL,gL6,gK1,&
                        gK2,1,qterm2,qterm1)
                   enddo
                   Dplus(gL1, gL2) = Dplus(gL1, gL2) + surf
                 enddo
               enddo
             enddo
           enddo
         enddo
      else if (mt == 1) then
         do gL1 = 1, calsize
           do gL2 = 1, calsize
             do gK1 = gL1, gL1
               do gK2 = gL2, gL2
                 Dplus(gL1, gL2) = Dplus(gL1,gL2) + &
                   RadialIntegral(iend,n,ic,is,gL1,gL2,0,gK1,gK2,1,0,qterm2,mt) &
                  -RadialIntegral(iend,n,ic,is,gL1,gL2,0,gK1,gK2,1,1,qterm2,mt)
               enddo
             enddo
           enddo
         enddo
      endif
!     ------------------------------------------------------------------
      call zgemm('n', 'n', calsize, calsize, calsize, -2.0*sqrt(2.0)*SQRTm1, &
           Dplus, calsize, sine_mat, calsize, CZERO, temp, calsize)
!     ------------------------------------------------------------------
      call zgemm('n', 'n', calsize, calsize, calsize, global_energy, &
           sine_transpose, calsize, temp, calsize, CZERO, &
           jspace(1:calsize,1:calsize,n,ic,is,dir), calsize)
!     ------------------------------------------------------------------
   enddo

   end subroutine calJxFromPhiLrIntegral
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calJyzFromJx(n, ic, is)
!  ===================================================================

   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, ic, is
   integer (kind=IntKind) :: dir, gL, gLp, gLpp, lpp, &
                      mpp, lp, mp, l, m, calsize
   
   calsize = master_size
   do gLp = 1, calsize
     mp = mofk(gLp)
     do gL = 1, calsize
       m = mofk(gL)
       jspace(gL, gLp, n, ic, is, 2) =  &
         -SQRTm1*(m - mp)*jspace(gL, gLp, n, ic, is, 1)
     enddo
   enddo

   do gLp = 1, calsize 
     lp = lofk(gLp)
     mp = mofk(gLp)
     do gL = 1, calsize
       l = lofk(gL)
       m = mofk(gL)
       do gLpp = 1, calsize
          lpp = lofk(gLpp)
          mpp = mofk(gLpp)
          jspace(gL, gLp, n, ic,is,3) = jspace(gL, gLp, n, ic, is, 3) + (0.5/SQRTm1)*( &
           (sqrt(lpp*(lpp + 1.0) - mpp*(mpp + 1.0))*kdelta(m,mpp+1) + &
           sqrt(lpp*(lpp + 1.0) - mpp*(mpp - 1.0))*kdelta(m,mpp-1))*kdelta(l,lpp)*jspace(gLpp,gLp,n,ic,is,2) &
          - (sqrt(lp*(lp + 1.0) - mp*(mp + 1.0))*kdelta(mpp,mp+1) + &
           sqrt(lp*(lp + 1.0) - mp*(mp - 1.0))*kdelta(mpp,mp-1))*kdelta(lpp,lp)*jspace(gL,gLpp,n,ic,is,2))
       enddo
       do dir = 1, 3
          jspace2(gL,gLp,n,ic,is,dir) = (-1.0)**lp * jspace(gL,gLp,n,ic,is,dir)
          jspace3(gL,gLp,n,ic,is,dir) = (-1.0)**l * jspace(gL,gLp,n,ic,is,dir)
          jspace4(gL,gLp,n,ic,is,dir) = (-1.0)**(l+lp) * jspace(gL,gLp,n,ic,is,dir) 
       enddo
     enddo
   enddo

   end subroutine calJyzFromJx
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calCurrentMatrix(n, is, eval, pot_type, mode)
!  ===================================================================

   use RadialGridModule, only : getRmesh, getRadialGridRadius
   use StepFunctionModule, only : interpolateStepFunction, getRadialStepFunction

   integer (kind=IntKind), intent(in) :: n, is, pot_type, mode
   complex (kind=CmplxKind), intent(in) :: eval

   integer (kind=IntKind) :: iend, ic, dir
   real (kind=RealKind) :: rmt
   real(kind=RealKind), pointer :: radial_grid(:)
   complex(kind=CmplxKind), allocatable :: sf_term(:,:,:), &
                               sf_single(:,:), sfqsum(:)

   rmt = getRadialGridRadius(n, MT=.true., nr=iend)
   radial_grid => getRmesh(n)

   allocate(sf_term(iend, kmax_sigma_2, kmax_sigma_2), sfqsum(iend), &
     sf_single(iend, jofk(kmax_sigma)))

!  ---------------------------------------------------------------- 
   call interpolateStepFunction(n, iend, radial_grid(1:iend), &
         kmax_sigma_2, kmax_sigma_2, sf_term)
!  ----------------------------------------------------------------
   call getRadialStepFunction(n, 1, iend, radial_grid, &
        lmax_sigma, sf_single)
!  ----------------------------------------------------------------
   
   do ic = 1, num_species(n)
!     -------------------------------------------------------------
      call calJxFromPhiLrIntegral(n,ic,eval,iend, &
           is,sf_single,sf_term,mt=pot_type)
!     -------------------------------------------------------------
      call calJyzFromJx(n, ic, is)
!     -------------------------------------------------------------
   enddo

   if (mode == 3) then
     do ic = 1, num_species(n)
!      ------------------------------------------------------    
       call calJtildeCPA(n, ic, is)
!      ------------------------------------------------------
     enddo
   else if (mode == 4) then
     do ic = 1, num_species(n)
!      ------------------------------------------------------
       call calJtildeSRO(n, ic, is)
!      ------------------------------------------------------
     enddo
   endif

   end subroutine calCurrentMatrix
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calJtildeSRO(n, ic, is)
!  ===================================================================
   
   use SROModule, only : getSROMatrix
   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, ic, is
   integer (kind=IntKind) :: nsize, dir
   complex (kind=CmplxKind), pointer :: Ta(:,:), Tc(:,:), tauc(:,:)
   complex (kind=CmplxKind), allocatable :: taucc(:,:), Tcc(:,:), Tac(:,:)
   complex (kind=CmplxKind), allocatable :: D(:,:), Dt(:,:), D1(:,:), Dt1(:,:)
   complex (kind=CmplxKind), allocatable :: D00(:,:), Dt00(:,:), D100(:,:), Dt100(:,:)
   complex (kind=CmplxKind), allocatable :: tmp1(:,:), tmp2(:,:)
   complex (kind=CmplxKind), allocatable :: buf1(:,:), buf2(:,:), &
                                                 buf3(:,:), buf4(:,:)

   Ta => getSROMatrix('blk-tinv', n, ic, is, nsize)
   Tc => getSROMatrix('blk-tinv', n, 0, is)
   tauc => getSROMatrix('blk-tau', n, 0, is)
  
!  call writeMatrix('Ta', Ta, nsize, nsize)
!  call writeMatrix('Tc', Tc, nsize, nsize)
!  call writeMatrix('tauc', tauc, nsize, nsize)

   allocate(Tac(nsize, nsize), Tcc(nsize, nsize), taucc(nsize, nsize))
   allocate(D(nsize, nsize), Dt(nsize, nsize), D1(nsize, nsize), &
            Dt1(nsize, nsize), tmp1(nsize, nsize), tmp2(nsize, nsize))
   allocate(buf1(kmax_kkr_max, kmax_kkr_max), buf2(kmax_kkr_max, kmax_kkr_max), &
     buf3(kmax_kkr_max, kmax_kkr_max), buf4(kmax_kkr_max, kmax_kkr_max), &
     D00(kmax_kkr_max, kmax_kkr_max), Dt00(kmax_kkr_max, kmax_kkr_max), &
     D100(kmax_kkr_max, kmax_kkr_max), Dt100(kmax_kkr_max, kmax_kkr_max))

   Tac = conjg(Ta)
   Tcc = conjg(Tc)
   taucc = conjg(tauc)

   D = CZERO; Dt = CZERO; D1 = CZERO; Dt1 = CZERO;
   D00 = CZERO; Dt00 = CZERO; D100 = CZERO; Dt100 = CZERO   

   tmp1 = Ta - Tc
   tmp2 = Tac - Tcc
   call computeAprojB('N', nsize, tauc, tmp1, D)
   call computeAprojB('N', nsize, tmp1, tauc, Dt)
   call computeAprojB('N', nsize, taucc, tmp2, D1)
   call computeAprojB('N', nsize, tmp2, taucc, Dt1)

   D00 = D(1:kmax_kkr_max, 1:kmax_kkr_max)
   D100 = D1(1:kmax_kkr_max, 1:kmax_kkr_max)
   Dt00 = Dt(1:kmax_kkr_max, 1:kmax_kkr_max)
   Dt100 = Dt1(1:kmax_kkr_max, 1:kmax_kkr_max)
!  call writeMatrix('D', D, nsize, nsize)
!  call writeMatrix('D1', D1, nsize, nsize)
!  call writeMatrix('Dt', Dt, nsize, nsize)
!  call writeMatrix('Dt1', Dt1, nsize, nsize)

   do dir = 1, 3
     buf1 = CZERO; buf2 = CZERO; buf3 = CZERO; buf4 = CZERO
!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     jspace(:,:,n,ic,is,dir), kmax_kkr_max, D00, &
     kmax_kkr_max, CZERO, buf1, kmax_kkr_max)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     Dt00, kmax_kkr_max, buf1, kmax_kkr_max, &
     CZERO, jtspace(:,:,n,ic,is,dir), kmax_kkr_max)
!    -----------------------------------------------------------------

!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     jspace2(:,:,n,ic,is,dir), kmax_kkr_max, D100, &
     kmax_kkr_max, CZERO,buf2, kmax_kkr_max)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     Dt00, kmax_kkr_max, buf2, kmax_kkr_max, &
     CZERO, jtspace2(:,:,n,ic,is,dir), kmax_kkr_max)
!    -----------------------------------------------------------------

!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     jspace3(:,:,n,ic,is,dir), kmax_kkr_max, D00, &
     kmax_kkr_max, CZERO, buf3, kmax_kkr_max)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     Dt100, kmax_kkr_max, buf3, kmax_kkr_max, &
     CZERO, jtspace3(:,:,n,ic,is,dir), kmax_kkr_max)
!    -----------------------------------------------------------------
    
!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     jspace4(:,:,n,ic,is,dir), kmax_kkr_max, D100, &
     kmax_kkr_max, CZERO, buf4, kmax_kkr_max)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', kmax_kkr_max, kmax_kkr_max, kmax_kkr_max, CONE, &
     Dt100, kmax_kkr_max, buf4, kmax_kkr_max, &
     CZERO, jtspace4(:,:,n,ic,is,dir), kmax_kkr_max)
!    -----------------------------------------------------------------
   enddo

   end subroutine calJtildeSRO
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calJtildeCPA(n, ic, is)
!  ===================================================================
   use SSSolverModule, only : getScatteringMatrix
   use CPAMediumModule, only : getSingleSiteTmat, getCPAMatrix
   use SROModule, only : getSROMatrix
   use MatrixInverseModule, only : MtxInv_LU
   use MatrixModule, only : computeAprojB
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: n, ic, is
   integer (kind=IntKind) :: dsize, i, dir
   complex (kind=CmplxKind), pointer :: t_ctemp(:,:), t_atemp(:,:), tau_ctemp(:,:)
   complex (kind=CmplxKind), allocatable :: t_c(:,:), t_a(:,:), tau_c(:,:) 
   complex (kind=CmplxKind), allocatable :: t_cc(:,:), t_ac(:,:), tau_cc(:,:)
   complex (kind=CmplxKind), allocatable :: Dt(:,:), D(:,:),  &
                                     D1(:,:), Dt1(:,:)
   complex (kind=CmplxKind), allocatable :: temp(:,:), temp2(:,:), temp3(:,:), temp4(:,:)


   t_atemp => getScatteringMatrix('TInv-Matrix', spin=is, site=n, atom=ic)
   t_ctemp => getSingleSiteTmat('TInv-Matrix', spin=is, site=n, atom=0)
   tau_ctemp => getCPAMatrix('Tau',site=n,atom=0)

   dsize = master_size

   allocate(t_c(dsize, dsize), t_a(dsize, dsize), tau_c(dsize, dsize))
   allocate(t_cc(dsize, dsize), t_ac(dsize, dsize), tau_cc(dsize, dsize))

   t_cc = CZERO; tau_cc = CZERO

   t_c = t_ctemp(1:dsize, 1:dsize)
   t_a = t_atemp(1:dsize, 1:dsize)
   tau_c = tau_ctemp(1:dsize, 1:dsize)
   t_ac = conjg(t_a)
   tau_cc = conjg(tau_c)
   t_cc = conjg(t_c)

   allocate(D(dsize, dsize), Dt(dsize, dsize), D1(dsize, dsize), &
        Dt1(dsize, dsize))
   allocate(temp(dsize, dsize), temp2(dsize, dsize), &
        temp3(dsize, dsize), temp4(dsize, dsize))


   D = CZERO
   Dt = CZERO
   D1 = CZERO
   Dt1 = CZERO

!  do i = 1, dsize
!    D(i, i) = CONE
!    Dt(i, i) = CONE
!    D1(i, i) = CONE
!    Dt1(i, i) = CONE
!  enddo

   temp = CZERO; temp2 = CZERO; temp3 = CZERO; temp4 = CZERO

   temp = t_a - t_c
   temp2 = t_ac - t_cc 
 
!  -----------------------------------------------------------------
   call computeAprojB('N', dsize, temp, tau_c, Dt)
!  ----------------------------------------------------------------
   call computeAprojB('N', dsize, tau_c, temp, D)
!  ----------------------------------------------------------------
   call computeAprojB('N', dsize, temp2, tau_cc, Dt1)
   call computeAprojB('N', dsize, tau_cc, temp2, D1)
!  -----------------------------------------------------------------

!  call writeMatrix('D', D, dsize, dsize)
!  call writeMatrix('D1', D1, dsize, dsize)
!  call writeMatrix('Dt', Dt, dsize, dsize)
!  call writeMatrix('Dt1', Dt1, dsize, dsize)
   do dir = 1, 3
     temp = CZERO; temp2 = CZERO; temp3 = CZERO; temp4 = CZERO
!    Calculating Jtilde(E_F + id, E_F + id)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, jspace(:,:,n,ic,is,dir), &
       dsize, D, dsize, CZERO, temp, dsize)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, Dt, dsize, &
       temp, dsize, CZERO, jtspace(:,:,n,ic,is,dir), dsize)
!    -----------------------------------------------------------------

!    Calculating Jtilde(E_F + id, E_F - id) 
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, jspace2(:,:,n,ic,is,dir), &
       dsize, D1, dsize, CZERO, temp2, dsize)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, Dt, dsize, &
       temp2, dsize, CZERO, jtspace2(:,:,n,ic,is,dir), dsize)
!    -----------------------------------------------------------------

!    Calculating Jtilde(E_F - id, E_F + id)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, jspace3(:,:,n,ic,is,dir), &
       dsize, D, dsize, CZERO, temp3, dsize)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, Dt1, dsize, &
       temp3, dsize, CZERO, jtspace3(:,:,n,ic,is,dir), dsize)
!    -----------------------------------------------------------------

!    Calculating Jtilde(E_F - id, E_F - id)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, jspace4(:,:,n,ic,is,dir), &
       dsize, D1, dsize, CZERO, temp4, dsize)
!    -----------------------------------------------------------------
     call zgemm('n', 'n', dsize, dsize, dsize, CONE, Dt1, dsize, &
       temp4, dsize, CZERO, jtspace4(:,:,n,ic,is,dir), dsize)
!    -----------------------------------------------------------------

     if (dir == 1) then
!      call writeMatrix('Jtx', jtspace(:,:,n,1,is,1), dsize, dsize)
!      call writeMatrix('Jtx_2', jtspace2(:,:,n,ic,is,1), dsize, dsize)
!      call writeMatrix('Jtx_3', jtspace3(:,:,n,ic,is,1), dsize, dsize)
!      call writeMatrix('Jtx_4', jtspace4(:,:,n,ic,is,1), dsize, dsize)
     endif
   enddo

   end subroutine calJtildeCPA
!  ===================================================================

!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getJMatrix(n, ic, is, dir, en_type, tilde)  result(Jmat)
!  ===================================================================
  
   integer (kind=IntKind), intent(in) :: n, ic, is, dir, en_type, tilde
   complex (kind=CmplxKind), pointer :: Jmat(:,:)

   if (tilde == 0) then
     if (en_type == 1) then
       Jmat => jspace(:,:,n,ic,is,dir)
     else if (en_type == 2) then
       Jmat => jspace2(:,:,n,ic,is,dir)
     else if (en_type == 3) then
       Jmat => jspace3(:,:,n,ic,is,dir)
     else if (en_type == 4) then
       Jmat => jspace4(:,:,n,ic,is,dir)
     else
       call ErrorHandler('getJMatrix','Incorrect energy type (1-4)', en_type)
     endif
   else if (tilde == 1) then
     if (en_type == 1) then
       Jmat => jtspace(:,:,n,ic,is,dir)
     else if (en_type == 2) then
       Jmat => jtspace2(:,:,n,ic,is,dir)
     else if (en_type == 3) then
       Jmat => jtspace3(:,:,n,ic,is,dir)
     else if (en_type == 4) then
       Jmat => jtspace4(:,:,n,ic,is,dir)
     else
       call ErrorHandler('getJMatrix','Incorrect energy type (1-4)', en_type)
     endif
   else
     call ErrorHandler('getJMatrix','Incorrect value for tilde', tilde)
   endif

   end function getJMatrix
!  ===================================================================
end module CurrentMatrixModule
