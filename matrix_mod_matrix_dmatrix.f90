!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_matrix_dmatrix.f90 ....................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................
!.......................... type/class DMatrix procedures ..........................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................

!...................................................................................................................................
!.......................... Constructors ...........................................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = DMatrix(m,n,is,js,mtype,name), sets MAT % v = 0
!non-optional: m optional: n=m is=1 js=1 mtype='DG' name=''  (inherited from MatrixClass: matrix_metadata_put )
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_zero_d                                (m,n,is,js,mtype,name)   result(MAT)                     !:.:.:
  type(DMatrix)                         :: MAT
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name
  integer                               :: stat

  call matrix_metadata_put(MAT,m,n,is,js,mtype,name)

  allocate(MAT % v(MAT % is:MAT % ie,                     &
                   MAT % js:MAT % je)      ,    stat=stat)

  if(stat > 0) call matrix_error('matrix_construct_zero_d: allocation failure for MAT % v')

  MAT % v = 0.0_dp

 end function       matrix_construct_zero_d
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = DMatrix(C,is=is,js=js,mtype='DG',name='MMAT'); If mtype='DS' the matrix is made symmetric. is,js take default values 1,1.
!TODO: maybe make C allocatable, so that is,js are passed through to the routine (Fortran 2003 feature) 
!WARNING: usage of this constructor should be discouraged. Better use: MAT = Matrix(....), MAT % v = C
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_array2_d                              (C,is,js,mtype,name)     result(MAT)                     !:.:.:
  type(DMatrix)                         :: MAT
  real     (dp)         , intent(in)    :: C(:,:)
  integer     , optional, intent(in)    :: is, js
  character(*), optional, intent(in)    :: mtype,name
  integer                               :: m,n,i,j

  m   = size  (C,1) ; n  = size  (C,2)

  i = 1; if(present(is)) i = is
  j = i; if(present(js)) j = js                                                                                                     ! if js is absent, then js=is

  MAT = matrix_construct_zero_d  (m,n,is=i,js=j,mtype=mtype,name=name)

  MAT % v(i:,j:) = C

  if(MAT % mtype(2:2) == 'S') call MAT % symmetric_set

 end function       matrix_construct_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = DMatrix(r,n)  calls matrix_construct_zero_d, same defaults, MAT % v = r  non-optional: r,m
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_real_d                                (r,m,n,is,js,mtype,name) result(MAT)                     !:.:.:
  type(DMatrix)                         :: MAT
  real     (dp)         , intent(in)    :: r
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name

  MAT     = matrix_construct_zero_d(m,n,is,js,mtype,name)

  MAT % v = r

 end function       matrix_construct_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = DMatrix(r,n)  calls matrix_construct_zero_d, same defaults, MAT % v = r  non-optional: r,m
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_complex_d                             (r,m,n,is,js,mtype,name) result(MAT)                     !:.:.:
  type(DMatrix)                         :: MAT
  complex  (dp)         , intent(in)    :: r
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name

  MAT     = matrix_construct_zero_d(m,n,is,js,mtype,name)

  MAT % v = REAL(r,kind=dp)

 end function       matrix_construct_complex_d
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = DMatrix('uniform',n,mtype='DS')  MAT = DMatrix('gaussian',n,sigma=PI)     non-optional: rtype, m  
!calls matrix_construct_zero_d, Matrix Class: matrix_random_set, matrix_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_random_d                    (rtype,m,n,is,js,mtype,name,sigma) result(MAT)                     !:.:.:
  type(DMatrix)                         :: MAT
  character(*)            , intent(in)  :: rtype  ! random, uniform, gaussian
  integer                 , intent(in)  :: m
  integer      , optional , intent(in)  :: n,is,js
  character(*) , optional , intent(in)  :: mtype,name
  real     (dp), optional , intent(in)  :: sigma
  real     (dp)                         :: s
  character(mtype_len)                  :: mt

  s   = 1.0_dp; if(present(sigma)) s  = sigma
  mt  = 'DG'  ; if(present(mtype)) mt = mtype

  MAT = matrix_construct_zero_d(m,n,is,js,mt,name)

  if( len(rtype) < 1 ) call matrix_error('matrix_construct_random_d: empty rtype')                                                              ! 261125 (MINOR-07)

  select case (rtype(1:1))
  case('r','R', 'u', 'U')
   call matrix_random_set  (MAT)                                                                                                    ! Symmetrization determined by MAT % mtype automatically by matrix_random_set
  case('g','G')
   call matrix_gaussian_set(MAT,sigma=s)                                                                                            ! Symmetrization determined by MAT % mtype automatically by matrix_gaussian_set 
  case('1')  ! a unit matrix!
   if( MAT % m /= MAT % n ) call matrix_error('matrix_construct_random_d: construct unit matrix: m/=n ')
   MAT % v = diagonal(1.0_dp,m)
  case default
   call matrix_error('matrix_construct_random_d: wrong rtype: '//trim(rtype))
  end select

 end function       matrix_construct_random_d
!-----------------------------------------------------------------------------------------------------------------------------------

!...................................................................................................................................
!.......................... Components   ...........................................................................................
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call MAT % symmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_symmetric_set_d                                 (MAT,uplo)                                               !:.:.:
  class(DMatrix)                        :: MAT
  character(*),optional      ,intent(in):: uplo

  if(.not.allocated(MAT % v)) call matrix_error('matrix_symmetric_set_d: MAT % v not allocated')

  if( MAT % m    /= MAT % n ) call matrix_error('matrix_symmetric_set_d: matrix is not a square matrix')

  call symmetric_set( MAT % v , uplo )

  MAT % mtype  = 'DS'

 end subroutine     matrix_symmetric_set_d
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... Operators    ...........................................................................................!:.:.:
!...................................................................................................................................


! 261125 (MINOR-02): removed unused matrix_assignFrom_matrix_d (not in any generic interface; intrinsic assignment is used instead)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_real_d                               (MATB,r)                                                 !:.:.:
  type   (DMatrix), intent(in out)      :: MATB
  real   (dp)    , intent(in)           :: r

  if( .not. allocated(MATB % v )) call matrix_error('matrix_assignFrom_real_d: MATB not allocated')

  MATB % v = r

 end subroutine     matrix_assignFrom_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_complex_d                            (MATB,r)                                                 !:.:.:
  type   (DMatrix), intent(in out)      :: MATB
  complex(dp)     , intent(in)          :: r

  if( .not. allocated(MATB % v )) call matrix_error('matrix_assignFrom_complex_d: MATB not allocated')

  MATB % v     = REAL(r,kind=dp)

 end subroutine     matrix_assignFrom_complex_d
!-----------------------------------------------------------------------------------------------------------------------------------
! matrix name not altered; mtype is reset to 'DG' (261125, BUG-05): call MATB % symmetric_set afterwards if C is symmetric
! WARNING: since C is not allocatable, when MATB  % m == 0, is=js=1 always  (change that in the future?)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_array2_d                             (MATB,C)                                                 !:.:.:
  type   (DMatrix), intent(in out)      :: MATB
  real   (dp)     , intent(in)          :: C(:,:)

  if(MATB  % m == 0) then               !if MATB has not been constructed

   MATB    % m  = size  (C,1) ;  MATB % n  = size  (C,2) 
   MATB    % is = lbound(C,1) ;  MATB % js = lbound(C,2)
   MATB    % ie = ubound(C,1) ;  MATB % je = ubound(C,2)
   MATB    % mtype = 'DG'     ;  MATB % name = ''                ! 261125 (BUG-11): DMatrix, was 'ZG'
   MATB    % v  = C                     ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated)

  else

   if(MATB % m /= size(C,1) .or. MATB % n /= size  (C,2) ) &
        call matrix_error('matrix_assignFrom_array2_d: MATB and C non conformable.')

   if(.not. allocated(MATB % v)) allocate(MATB  % v(MATB % is : MATB % ie, MATB % js : MATB % je))
   MATB                                         % v(MATB % is : MATB % ie, MATB % js : MATB % je)  = C
   MATB % mtype = 'DG'                  ! 261125 (BUG-05): C is not symmetric in general (was: kept 'DS')

  end if

 end subroutine     matrix_assignFrom_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_plus_matrix_d                                     (r,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je          

 end function       real_plus_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_plus_real_d                                     (MATA,r)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je              

 end function       matrix_plus_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_plus_array2_d                                   (MATA,C)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: C(:,:)
  type   (DMatrix)                      :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('matrix_plus_array2_d: MATA and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + C               ! because of MATB = MATA, MATB gets the same is,ie,js,je 

  MATB % mtype = 'DG'

 end function       matrix_plus_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_plus_matrix_d                                   (C,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: C(:,:)
  type   (DMatrix)                      :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('array2_plus_matrix_d: MATA and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + C               ! because of MATB = MATA, MATB gets the same is,ie,js,je     

  MATB % mtype = 'DG'


 end function       array2_plus_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_plus_matrix_d                                   (MATA,MATB)              result(MATC)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA,MATB
  type   (DMatrix)                      :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_plus_matrix_d: MATA and MATB not conformable')

  MATC            = MATA                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'DG'               ! if MATB is not symmetric, the result is not symmetric

  if(MATA % mtype == 'DS' .and. MATB % mtype == 'DS') MATC % mtype = 'DS'

  MATC % v        = MATA % v  + MATB % v                 

 end function       matrix_plus_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_subtract_matrix_d                                 (r,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v =-MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je          

 end function       real_subtract_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_subtract_real_d                                 (MATA,r)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v - r               ! because of MATB = MATA, MATB gets the same is,ie,js,je              

 end function       matrix_subtract_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_subtract_array2_d                               (MATA,C)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)           :: MATA
  real   (dp)     , intent(in)           :: C(:,:)
  type   (DMatrix)                       :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('matrix_subtract_array2_d: MATA and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v - C               ! because of MATB = MATA, MATB gets the same is,ie,js,je                

  MATB % mtype = 'DG'                                                                                                              ! 261125 (BUG-04): as in the '+' versions (was: kept 'DS')

 end function       matrix_subtract_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_subtract_matrix_d                               (C,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: C(:,:)
  type   (DMatrix)                      :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('array2_subtract_matrix_d: MATA and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v =-MATA % v + C               ! because of MATB = MATA, MATB gets the same is,ie,js,je     

  MATB % mtype = 'DG'                                                                                                              ! 261125 (BUG-04): as in the '+' versions (was: kept 'DS')

 end function       array2_subtract_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_subtract_matrix_d                               (MATA,MATB)              result(MATC)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA,MATB
  type   (DMatrix)                      :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_subtract_matrix_d: MATA and MATB not conformable')

  MATC            = MATA                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'DG'               ! if MATB is not Hermitian, the result is not Hermitian

  if(MATA % mtype == 'DS' .and. MATB % mtype == 'DS') MATC % mtype = 'DS'

  MATC % v        = MATA % v  - MATB % v                 

 end function       matrix_subtract_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_minus_matrix_d                           (MATA)                   result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  type   (DMatrix)                      :: MATB

  MATB = MATA

  MATB % v = -MATB % v

 end function       matrix_return_minus_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_mult_matrix_d                                     (r,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v * r                 

 end function       real_mult_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_real_d                                     (MATA,r)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v * r                 

 end function       matrix_mult_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_mult_matrix_d                                  (z,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  complex(dp)     , intent(in)          :: z
  type   ( Matrix)                      :: MATB
  real   (dp)                           :: x,y

  MATB % m     = MATA % m       ; MATB % n    = MATA % n
  MATB % is    = MATA % is      ; MATB % ie   = MATA % ie
  MATB % js    = MATA % js      ; MATB % je   = MATA % je

  MATB % mtype = 'ZG'           ; MATB % name = MATA % name                                                                         ! 261125 (MINOR-06): keep name, was ''

  allocate(MATB % v (MATB % is  : MATB % ie   , MATB % js : MATB % je))

  x            = real(z,kind=dp); y           = aimag(z)

  if      ( x ==        0.0_dp .and. y == 1.0_dp           )then                                                                    ! avoid some multiplications
   MATB % v   =  cmplx(      0.0_dp,      MATA % v, kind=dp)
  else if ( x ==        1.0_dp .and. y == 0.0_dp           )then
   MATB % v   =  cmplx(    MATA % v,   0.0_dp     , kind=dp)
  else if ( x ==        0.0_dp                             )then
   MATB % v   =  cmplx( 0.0_dp     ,  y * MATA % v, kind=dp)
  else if (                          y == 0.0_dp           )then
   MATB % v   =  cmplx(x * MATA % v,      0.0_dp  , kind=dp)
  else
   MATB % v   =  z *   MATA % v
  end if

 end function       complex_mult_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_complex_d                                  (MATA,z)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  complex(dp)     , intent(in)          :: z
  type   ( Matrix)                      :: MATB
  real   (dp)                           :: x,y

  MATB % m     = MATA % m       ; MATB % n    = MATA % n
  MATB % is    = MATA % is      ; MATB % ie   = MATA % ie
  MATB % js    = MATA % js      ; MATB % je   = MATA % je

  MATB % mtype = 'ZG'           ; MATB % name = MATA % name                                                                         ! 261125 (MINOR-06): keep name, was ''

  allocate(MATB % v (MATB % is  : MATB % ie   , MATB % js : MATB % je))

  x            = real(z,kind=dp); y           = aimag(z)

  if      ( x ==        0.0_dp .and. y == 1.0_dp           )then                                                                    ! avoid some multiplications
   MATB % v   =  cmplx(      0.0_dp,      MATA % v, kind=dp)
  else if ( x ==        1.0_dp .and. y == 0.0_dp           )then
   MATB % v   =  cmplx(    MATA % v,   0.0_dp     , kind=dp)
  else if ( x ==        0.0_dp                             )then
   MATB % v   =  cmplx( 0.0_dp     ,  y * MATA % v, kind=dp)
  else if (                          y == 0.0_dp           )then
   MATB % v   =  cmplx(x * MATA % v,      0.0_dp  , kind=dp)
  else
   MATB % v   =  z *   MATA % v
  end if

 end function       matrix_mult_complex_d
!-----------------------------------------------------------------------------------------------------------------------------------
! Non-conformable shapes (inner dimensions differ) abort the program   ! 261125 (CONV-01; was MINOR-03: LMATMUL NaN-filled result)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_mult_array2_d                                   (MATA,C)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: C(:,:)
  type   (DMatrix)                      :: MATB

  
  if( MATA % n /= size(C,1) ) call matrix_error('matrix_mult_array2_d: non-conformable matrices (inner dimensions differ)')   ! 261125 (CONV-01): shape mismatch aborts (was: NaN-filled result)

  MATB % m     = MATA % m      ; MATB % n    = size  (C,2)
  MATB % is    = MATA % is     ; MATB % ie   = MATA % ie
  MATB % js    = lbound(C,2)   ; MATB % je   = ubound(C,2)

  MATB % mtype = 'DG'          ; MATB % name = MATA % name                                                                          ! not Hermitian in general

  allocate(MATB % v (MATB % is : MATB % ie   , MATB % js : MATB % je))

  MATB % v = LMATMUL(MATA % v, C, mtype = MATA % mtype)                 

 end function       matrix_mult_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
! Non-conformable shapes (inner dimensions differ) abort the program   ! 261125 (CONV-01; was MINOR-03: LMATMUL NaN-filled result)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_mult_matrix_d                                   (C,MATA)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: C(:,:)
  type   (DMatrix)                      :: MATB

  if( size(C,2) /= MATA % m ) call matrix_error('array2_mult_matrix_d: non-conformable matrices (inner dimensions differ)')   ! 261125 (CONV-01): shape mismatch aborts (was: NaN-filled result)

  MATB % m     = size  (C,1)   ; MATB % n    = MATA % n 
  MATB % is    = lbound(C,1)   ; MATB % ie   = ubound(C,1)
  MATB % js    = MATA % js     ; MATB % je   = MATA % je

  MATB % mtype = 'DG'          ; MATB % name = MATA % name                ! 261125 (BUG-11): DMatrix, was 'ZG'

  allocate(MATB % v (MATB % is : MATB % ie   , MATB % js : MATB % je))

  MATB % v = LMATMUL(C, MATA % v)                 

 end function       array2_mult_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
! Non-conformable shapes (inner dimensions differ) abort the program   ! 261125 (CONV-01; was MINOR-03: LMATMUL NaN-filled result)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_mult_matrix_d                                   (MATA,MATB)              result(MATC)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA,MATB
  type   (DMatrix)                      :: MATC

  if( MATA % n /= MATB % m ) call matrix_error('matrix_mult_matrix_d: non-conformable matrices (inner dimensions differ)')   ! 261125 (CONV-01): shape mismatch aborts (was: NaN-filled result)

  MATC % m  = MATA % m         ; MATC % n     = MATB % n
  MATC % is = MATA % is        ; MATC % ie    = MATA % ie
  MATC % js = MATB % js        ; MATC % je    = MATB % je

  MATC    % mtype =  'DG'      ; MATC % name  = MATA % name 

  allocate(MATC % v (MATC % is : MATC % ie , MATC % js : MATC % je))

  MATC % v  = LMATMUL(MATA % v,   MATB % v, mtype = MATA % mtype)                 

 end function       matrix_mult_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_divide_real_d                                   (MATA,r)                 result(MATB)                    !:.:.:
  class  (DMatrix), intent(in)          :: MATA
  real   (dp)     , intent(in)          :: r
  type   (DMatrix)                      :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v / r                 

 end function       matrix_divide_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_mult_matrix_sub_d                               (MATA,MATB,MATC)                                         !:.:.:
  class (DMatrix), intent(in)           :: MATA,MATB
  type  (DMatrix), intent(inout)        :: MATC

  if( MATA % n /= MATB % m ) error stop 'matrix_mult_matrix_sub_d: non-conformable matrices'   ! 261125 (CONV-01): inner-dimension mismatch aborts (was: NaN-filled output)

  if( .not. allocated(MATC % v) .or. MATC % m /= MATA % m .or. MATC % n /= MATB % n )then                                         ! 261125 (BUG-08): (re)construct the output
   if( allocated(MATC % v) ) deallocate(MATC % v)                                                                                  !                    (was: segfault if unallocated, silent NaN if wrong shape)
   MATC % m  = MATA % m        ; MATC % n  = MATB % n
   MATC % is = MATA % is       ; MATC % ie = MATA % ie
   MATC % js = MATB % js       ; MATC % je = MATB % je
   allocate(MATC % v(MATC % is : MATC % ie, MATC % js : MATC % je))
  end if

  call  mmmult(MATA % v, MATB % v, MATC % v, mtype= MATA % mtype)

  MATC % mtype = 'DG'                                                                                                              ! 261125 (BUG-08): A.B is not Hermitian/symmetric in general (was: kept old mtype)

 end subroutine     matrix_mult_matrix_sub_d


!...................................................................................................................................
!.......................... Linear Algebra with LAPACK .............................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_inverse                                        (MATA)                   result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MATA
  type (DMatrix)                        :: MATB

  call metadata_copy(MATA,MATB)

  allocate(MATB % v,mold=MATA % v)

  MATB % v = inverse(MATA % v)

 end  function      dmatrix_inverse
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_determinant                                    (MAT)                    result(z)                       !:.:.:
  class(DMatrix), intent(in)            :: MAT
  real   (dp)                           :: z

  z = determinant(MAT % v)

 end  function      dmatrix_determinant
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_lndet                                          (MAT)                    result(z)                       !:.:.:
  class(DMatrix), intent(in)            :: MAT
  real   (dp)                           :: z(2)

  z = lndet(MAT % v)

 end  function      dmatrix_lndet
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_eigenvalues                                    (MAT)                    result(vec)                     !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type ( Vector)                        :: vec

  if( MAT  % m /= MAT  % n) call matrix_error('dmatrix_eigenvalues: MAT not a square matrix') ! 261125 (BUG-07): was m /= m (always false)

  vec     =      Vector( MAT % n, is= MAT % is)

  vec % v = eigenvalues( MAT % v, mtype = MAT % mtype)

 end  function      dmatrix_eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_eigenvectors                                   (MATA,vec)               result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MATA
  type ( Vector), intent(in out)        :: vec
  type ( Matrix)                        :: MATB
  complex(dp)   , allocatable           :: W(:,:)
  integer                               :: n

  if( MATA % m /= MATA % n) call matrix_error('dmatrix_eigenvectors: MATA not a square matrix') ! 261125 (BUG-07): was m /= m (always false)

  n        =      MATA % n ; allocate(W(n,n+1))

  if( vec % n /= n ) vec      =      Vector( n, is= MATA % is)

  call metadata_copy(MATA,MATB);  allocate(MATB % v(MATA % is : MATA % ie, MATA % js : MATA % je))
  
  W        =      eigenvectors(MATA % v, mtype = MATA % mtype)                                                                      ! 261125 (BUG-02): pass mtype -> dsyev for 'DS' (was: always dgeev)
  
  vec  % v =      W(:,1)
  
  MATB % v =      W(:,2:)

  deallocate(W)

 end  function      dmatrix_eigenvectors


!...................................................................................................................................
!.......................... Misc Procedures.........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         dmatrix_antisymmetric_set                              (MAT,uplo)                                               !:.:.:
  class(DMatrix)                        :: MAT
  character(*),optional      ,intent(in):: uplo

  if(.not.allocated(MAT % v)) call matrix_error('dmatrix_antisymmetric_set: MAT % v not allocated')                                   ! 261125 (BUG-06)

  if( MAT % m    /= MAT % n ) call matrix_error('dmatrix_antisymmetric_set: matrix is not a square matrix')                         ! 261125 (BUG-06)

  call antisymmetric_set( MAT % v , uplo)

  MAT % mtype  = 'DG'                                                                                                              ! 261125 (BUG-06): antisymmetric is not symmetric (was: kept 'DS')

 end subroutine     dmatrix_antisymmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_diagonal_get                                   (MAT)                    result(vec)                     !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DVector)                        :: vec
  
  vec     = DVector ( min(MAT % m, MAT % n) , is= MAT % is )                                                                     ! 261125 (BUG-12): min(m,n) elements, built as DVector (was: Vector(MAT % m) + conversion)
  vec % v = diagonal( MAT % v )

 end  function      dmatrix_diagonal_get
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_diagonal_set_from_dmatrix                      (MATA)                   result(MATB)                    !:.:.:
  class (DMatrix), intent(in)           :: MATA
  type  (DMatrix)                       :: MATB

  MATB     = MATA

  MATB % v = diagonalMatrix( MATA % v )

 end  function      dmatrix_diagonal_set_from_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_diagonal_set_from_dvector                      (vec)                    result(MAT)                     !:.:.:
  class (DVector), intent(in)           :: vec
  type  (DMatrix)                       :: MAT

  MAT     = DMatrix ( vec % n , is= vec % is)

  MAT % v = diagonal( vec % v )

 end  function      dmatrix_diagonal_set_from_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_traceless_get                                  (MAT)                    result(MATB)                    !:.:.:
  class (DMatrix), intent(in)           :: MAT
  type  (DMatrix)                       :: MATB

  call metadata_copy(MAT,MATB); allocate(MATB % v, mold = MAT % v)

  MATB % v =  traceless(MAT % v)

 end function       dmatrix_traceless_get
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dmatrix_traceless_set                                  (MAT)                                                    !:.:.:
  class (DMatrix), intent(in out)       :: MAT

  call traceless_set(MAT % v)

 end subroutine     dmatrix_traceless_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_trace2c                                        (MAT)                    result(r)                       !:.:.:
  class (DMatrix), intent(in)           :: MAT                                                                         ! 261125 (MINOR-05): class, was type
  real   (dp)                           :: r

  r = trace2c(MAT % v, mtype= MAT % mtype)

 end  function      dmatrix_trace2c
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_trace2                                         (MAT)                    result(r)                       !:.:.:
  class (DMatrix), intent(in)           :: MAT                                                                         ! 261125 (MINOR-05): class, was type
  real   (dp)                           :: r

  r = trace2(MAT % v, mtype= MAT % mtype)

 end  function      dmatrix_trace2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_trace                                          (MAT)                    result(r)                       !:.:.:
  class (DMatrix), intent(in)           :: MAT                                                                         ! 261125 (MINOR-05): class, was type
  real   (dp)                           :: r

  r = trace(MAT % v)

 end  function      dmatrix_trace
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_norm_d                                          (MAT)                    result(r)                       !:.:.:
  class  (DMatrix), intent(in)          :: MAT                                                                         ! 261125 (MINOR-05): class, was type
  real(dp)                              :: r
  
  r = array2_norm_d(MAT % v)

 end function       matrix_norm_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_dmatrix                                  (MAT)                                                    !:.:.:
  class  (DMatrix)                      :: MAT
  
  call matrix_random_set(MAT)

 end subroutine     random_number_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_dmatrix_gaussian                         (MAT,sigma)                                              !:.:.:
  class  (DMatrix)                      :: MAT
  real   (dp)                           :: sigma
  
  call matrix_gaussian_set(MAT,sigma=sigma)

 end subroutine     random_number_dmatrix_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_transpose_d                              (MATA)                   result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MATA
  type (DMatrix)                        :: MATB

  !Careful m <-> n
  call matrix_metadata_put(MATB,MATA % n,MATA % m, MATA % js, MATA % is, MATA % mtype, MATA % name)

  allocate(MATB % v( MATA % js : MATA % je , MATA % is : MATA % ie ))

  MATB % v = TRANSPOSE(MATA % v)

 end function       matrix_return_transpose_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_dmatrix_complex_return_matrix                  (MATA,MATB)              result(MATC)                    !:.:.:
  class(DMatrix), intent(in)            :: MATA, MATB
  type ( Matrix)                        :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n) &
       call matrix_error('dmatrix_dmatrix_complex_return_matrix: MATA, MATB not conformable.')

  call matrix_metadata_put(MATC,MATA % m,MATA % n, MATA % is, MATA % js, 'ZG' , MATA % name)

  allocate(MATC % v( MATA % is : MATA % ie , MATA % js : MATA % je ))

  MATC % v = CMPLX(MATA % v, MATB % v, kind=dp)

 end function       dmatrix_dmatrix_complex_return_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_is_symmetric                                   (MATA)                   result(r)                       !:.:.:
  class(DMatrix), intent(in)            :: MATA
  real (dp)                             :: r

  r = isSymmetric( MATA % v )

 end function       dmatrix_is_symmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_is_antisymmetric                               (MATA)                   result(r)                       !:.:.:
  class(DMatrix), intent(in)            :: MATA
  real (dp)                             :: r

  r = isAntiSymmetric( MATA % v )

 end function       dmatrix_is_antisymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_is_nan                                         (MAT)                    result(itis)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  logical                               :: itis

  itis = .false.

  if( isNaN(MAT % v) ) itis = .true.

 end  function      dmatrix_is_nan
!-----------------------------------------------------------------------------------------------------------------------------------



!...................................................................................................................................
!.......................... Math and Array Procedures .............................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_abs                                            (MAT)                    result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = abs(MAT % v)

 end  function      dmatrix_abs
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_sin                                             (MAT)                    result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = sin(MAT % v)

 end  function      dmatrix_sin
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_cos                                             (MAT)                    result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = cos(MAT % v)

 end  function      dmatrix_cos
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_exp                                             (MAT)                    result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = exp(MAT % v)

 end  function      dmatrix_exp
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_log                                             (MAT)                    result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = log(MAT % v)

 end  function      dmatrix_log
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_sqrt                                            (MAT)                    result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = sqrt(MAT % v)

 end  function      dmatrix_sqrt
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_power_integer                                  (MAT,n)                  result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  integer       , intent(in)            :: n
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = (MAT % v)**n

 end  function      dmatrix_power_integer
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_power_real                                     (MAT,n)                  result(MATB)                    !:.:.:
  class(DMatrix), intent(in)            :: MAT
  real (dp)     , intent(in)            :: n
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = (MAT % v)**n

 end  function      dmatrix_power_real





!-----------------------------------------------------------------------------------------------------------------------------------
!  Copyright by Konstantinos N. Anagnostopoulos, Physics Department, National Technical University of Athens, 2020
!  konstant@mail.ntua.gr, www.physics.ntua.gr/konstant
!  
!  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
!  published by the Free Software Foundation, version 3 of the License.
!  
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public Liense along with this program. If not, see http://www.gnu.org/licenses
!-----------------------------------------------------------------------------------------------------------------------------------
