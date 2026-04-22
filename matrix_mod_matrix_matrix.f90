!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_matrix_matrix.f90 .....................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................
!.......................... type/class Matrix procedures ...........................................................................!:.:.:
!...................................................................................................................................

!...................................................................................................................................
!.......................... Constructors ...........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
!Basic constructor: allocates MAT % v.  Other constructors call matrix_construct_zero (m,n,is,js,mtype,name)
!Calls: MatrixClass: matrix_metadata_put(MAT,m,n,is,js,mtype,name), inherits defaults: n=m is=1 js=1 mtype='ZG' name=''
!Usage: MAT = Matrix(m,n,is,js,mtype,name), sets MAT % v = 0
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_zero                                  (m,n,is,js,mtype,name)   result(MAT)                     !:.:.:
  type(Matrix)                          :: MAT
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name
  integer                               :: stat

  call matrix_metadata_put(MAT,m,n,is,js,mtype,name)

  allocate(MAT % v(MAT % is:MAT % ie,                     &
                   MAT % js:MAT % je)      ,    stat=stat)

  if(stat > 0) call matrix_error('matrix_construct_zero: allocation failure for MAT % v')

  MAT % v = ZERO

 end function       matrix_construct_zero
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = Matrix(C,is=is,js=js,mtype='ZG',name='MMAT'); If mtype='ZH' the matrix is made hermitian.
!TODO: maybe make C allocatable, so that is,js are passed through to the routine (Fortran 2003 feature) 
!WARNING: usage of this constructor should be discouraged. Better use: MAT = Matrix(....), MAT % v = C
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_array2                                (C,is,js,mtype,name)     result(MAT)                     !:.:.:
  type(Matrix)                          :: MAT
  complex  (dp)         , intent(in)    :: C(:,:)
  integer     , optional, intent(in)    :: is, js
  character(*), optional, intent(in)    :: mtype,name
  integer                               :: m,n,i,j

  m   = size  (C,1) ; n  = size  (C,2)

  i = 1; if(present(is)) i = is
  j = i; if(present(js)) j = js                                                                                                     ! if js is absent, then js=is

  MAT = matrix_construct_zero  (m,n,is=i,js=j,mtype=mtype,name=name)

  MAT % v(i:,j:) = C

  if(MAT % mtype(2:2) == 'H') call MAT % hermitian_set

 end function       matrix_construct_array2
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = Matrix(z,n)  calls matrix_construct_zero, same defaults, MAT % v = c  non-optional: c,m, option 'ZH' sets MAT % v = REAL(c)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_complex                               (c,m,n,is,js,mtype,name) result(MAT)                     !:.:.:
  type(Matrix)                          :: MAT
  complex  (dp)         , intent(in)    :: c
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name
  integer                               :: stat

  MAT     = matrix_construct_zero(m,n,is,js,mtype,name)

  MAT % v = c

  if(MAT % mtype(2:2) == 'H') call MAT % hermitian_set                                                                              ! it makes the matrix effectively real

 end function       matrix_construct_complex
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = Matrix(r,n)  calls matrix_construct_zero, same defaults, MAT % v = r  non-optional: r,m
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_real                                  (r,m,n,is,js,mtype,name) result(MAT)                     !:.:.:
  type(Matrix)                          :: MAT
  real     (dp)         , intent(in)    :: r
  integer               , intent(in)    :: m
  integer     , optional, intent(in)    :: n,is,js
  character(*), optional, intent(in)    :: mtype,name
  integer                               :: stat

  MAT     = matrix_construct_zero(m,n,is,js,mtype,name)

  MAT % v = r

 end function       matrix_construct_real
!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: MAT = Matrix('uniform',n, mtype='ZH')  MAT = Matrix('gaussian',n,sigma=PI) => matrix_construct_random
!calls matrix_construct_zero, Matrix Class: matrix_random_set, matrix_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_construct_random                     (rtype,m,n,is,js,mtype,name, sigma) result(MAT)                     !:.:.:
  type(Matrix)                          :: MAT
  character(*)            , intent(in)  :: rtype  ! random, uniform, gaussian
  integer                 , intent(in)  :: m
  integer      , optional , intent(in)  :: n,is,js
  character(*) , optional , intent(in)  :: mtype,name
  real     (dp), optional , intent(in)  :: sigma
  real     (dp)                         :: s
  character(mtype_len)                  :: mt

  s   = 1.0_dp; if(present(sigma)) s  = sigma
  mt  = 'ZG'  ; if(present(mtype)) mt = mtype

  MAT = matrix_construct_zero(m,n,is,js,mt,name)

  select case (rtype(1:1))
  case('r','R', 'u', 'U')
   call matrix_random_set  (MAT)                                                                                                    ! Symmetrization: determined by MAT % mtype automatically by matrix_random_set        
  case('g','G')
   call matrix_gaussian_set(MAT,sigma=s)                                                                                            ! Symmetrization: determined by MAT % mtype automatically by matrix_gaussian_set 
  case('1')  ! a unit matrix!
   if( MAT % m /= MAT % n ) call matrix_error('matrix_construct_random: construct unit matrix: m/=n ')
   MAT % v = diagonal(ONE,m)
  case default
   call                          matrix_error('matrix_construct_random: wrong rtype: '//trim(rtype))
  end select

 end function       matrix_construct_random
!-----------------------------------------------------------------------------------------------------------------------------------

!...................................................................................................................................
!.......................... Components   ...........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call MAT % hermitian => matrix_hermitian_set  calls: matrix_mod_array: array2_hermitian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_hermitian_set                                   (MAT, uplo)                                              !:.:.:
  class(Matrix)                         :: MAT
  character(*),optional      ,intent(in):: uplo

  if(.not.allocated(MAT % v)) call matrix_error('matrix_hermitian_set: MAT % v not allocated') 

  if( MAT % m    /= MAT % n ) call matrix_error('matrix_hermitian_set: matrix is not a square matrix')

  call hermitian_set( MAT % v , uplo)

  MAT % mtype  = 'ZH'

 end subroutine     matrix_hermitian_set
!-----------------------------------------------------------------------------------------------------------------------------------

!...................................................................................................................................
!.......................... Operators    ...........................................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
! This is not necessary as it is written now: copy data components 'as is' is the default assignment for types. We leave it in case we want to change in the future
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_assignFrom_matrix                               (MATB,MATA)                                              !:.:.:
  type(Matrix), intent(in out)          :: MATB
  type(Matrix), intent(in    )          :: MATA
  character(100)                        :: err

  MATB % m     = MATA % m
  MATB % n     = MATA % n 
  MATB % is    = MATA % is 
  MATB % ie    = MATA % ie 
  MATB % js    = MATA % js  
  MATB % je    = MATA % je  
  MATB % mtype = MATA % mtype  
  MATB % name  = MATA % name

  MATB % v     = MATA % v               ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated)

 end subroutine     matrix_assignFrom_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_dmatrix                              (MATB,MATA)                                              !:.:.:
  type( Matrix), intent(in out)         :: MATB
  type(DMatrix), intent(in    )         :: MATA

  MATB     = Matrix(MATA % m, MATA % n , MATA % is, MATA % js, 'ZG')

  MATB % v =  CMPLX(MATA % v, 0.0_dp, kind = dp)

 end subroutine     matrix_assignFrom_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_real                                 (MATB,r)                                                 !:.:.:
  type   (Matrix), intent(in out)       :: MATB
  real   (dp)    , intent(in)           :: r

  if( .not. allocated(MATB % v )) call matrix_error('matrix_assignFrom_real: MATB not allocated')

  MATB % v = r

 end subroutine     matrix_assignFrom_real
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_complex                              (MATB,r)                                                 !:.:.:
  type   (Matrix), intent(in out)       :: MATB
  complex(dp)    , intent(in)           :: r

  if( .not. allocated(MATB % v )) call matrix_error('matrix_assignFrom_complex: MATB not allocated')

  MATB % v = r

  if( AIMAG(r) /= 0.0_dp) MATB % mtype = 'ZG'

 end subroutine     matrix_assignFrom_complex
!-----------------------------------------------------------------------------------------------------------------------------------
! matrix name and mtype not altered!
! WARNING: since C is not allocatable, when MATB  % m == 0, is=js=1 always  (change that in the future?)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_assignFrom_array2                               (MATB,C)                                                 !:.:.:
  type   (Matrix), intent(in out)       :: MATB
  complex(dp)    , intent(in)           :: C(:,:)

  if(MATB  % m == 0) then               !if MATB has not been constructed
   MATB    % m  = size  (C,1) ;  MATB % n  = size  (C,2) 
   MATB    % is = lbound(C,1) ;  MATB % js = lbound(C,2)
   MATB    % ie = ubound(C,1) ;  MATB % je = ubound(C,2)
   MATB    % mtype = 'ZG'     ;  MATB % name = ''
   MATB    % v  = C                     ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated)
  else
   if(MATB % m /= size(C,1) .or. MATB % n /= size  (C,2) )call matrix_error('matrix_assignFrom_array2: MATB and C non conformable.')
   if(.not. allocated(MATB % v)) allocate(MATB  % v(MATB % is : MATB % ie, MATB % js : MATB % je))
   MATB                                         % v(MATB % is : MATB % ie, MATB % js : MATB % je)  = C
  end if

 end subroutine     matrix_assignFrom_array2
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_plus_matrix                                       (r,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je          

 end function       real_plus_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_plus_real                                       (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je              

 end function       matrix_plus_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_plus_matrix                                    (r,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je              

 end function       complex_plus_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_plus_complex                                    (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je                 

 end function       matrix_plus_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_plus_array2                                     (MATA,C)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: C(:,:)
  type   (Matrix)                       :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('matrix_plus_array2: MATB and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + C               ! because of MATB = MATA, MATB gets the same is,ie,js,je 

  MATB % mtype = 'ZG'

 end function       matrix_plus_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_plus_matrix                                     (C,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: C(:,:)
  type   (Matrix)                       :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('array2_plus_matrix: MATB and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v + C               ! because of MATB = MATA, MATB gets the same is,ie,js,je     

  MATB % mtype = 'ZG'


 end function       array2_plus_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_plus_matrix                                     (MATA,MATB)              result(MATC)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA,MATB
  type   (Matrix)                       :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_plus_matrix: MATA and MATB not conformable')

  MATC            = MATA                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'ZG'               ! if MATB is not Hermitian, the result is not Hermitian

  if(MATA % mtype == 'ZH' .and. MATB % mtype == 'ZH') MATC % mtype = 'ZH'

  MATC % v        = MATA % v  + MATB % v                 

 end function       matrix_plus_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_plus_dmatrix                                    (MATA,MATB)              result(MATC)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  class (DMatrix), intent(in)           :: MATB
  type   (Matrix)                       :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_plus_matrix: MATA and MATB not conformable')

  MATC            = MATA                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'ZG'               

  if(MATA % mtype == 'ZH' .and. MATB % mtype == 'DS') MATC % mtype = 'ZH'

  MATC % v        = MATA % v  + MATB % v                 

 end function       matrix_plus_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_plus_matrix                                    (MATA,MATB)              result(MATC)                    !:.:.:
  class (DMatrix), intent(in)           :: MATA
  class  (Matrix), intent(in)           :: MATB
  type   (Matrix)                       :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_plus_matrix: MATA and MATB not conformable')

  MATC            = MATB                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'ZG'               

  if(MATB % mtype == 'ZH' .and. MATA % mtype == 'DS') MATC % mtype = 'ZH'

  MATC % v        = MATA % v  + MATB % v                 

 end function       dmatrix_plus_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_subtract_matrix                                   (r,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v =-MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je          

 end function       real_subtract_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_subtract_real                                   (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v - r               ! because of MATB = MATA, MATB gets the same is,ie,js,je              

 end function       matrix_subtract_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_subtract_matrix                                (r,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v =-MATA % v + r               ! because of MATB = MATA, MATB gets the same is,ie,js,je              

 end function       complex_subtract_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_subtract_complex                                (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v - r               ! because of MATB = MATA, MATB gets the same is,ie,js,je                 

 end function       matrix_subtract_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_subtract_array2                                 (MATA,C)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: C(:,:)
  type   (Matrix)                       :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('matrix_subtract_array2: MATB and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v - C               ! because of MATB = MATA, MATB gets the same is,ie,js,je                

 end function       matrix_subtract_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_subtract_matrix                                 (C,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: C(:,:)
  type   (Matrix)                       :: MATB

  if(MATA % m /= size(C,1) .or. MATA % n /= size(C,2))call matrix_error('array2_subtract_matrix: MATB and C non conformable.')

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v =-MATA % v + C               ! because of MATB = MATA, MATB gets the same is,ie,js,je     

 end function       array2_subtract_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_subtract_matrix                                 (MATA,MATB)              result(MATC)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA,MATB
  type   (Matrix)                       :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_subtract_matrix: MATA and MATB not conformable')

  MATC            = MATA                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'ZG'               ! if MATB is not Hermitian, the result is not Hermitian

  if(MATA % mtype == 'ZH' .and. MATB % mtype == 'ZH') MATC % mtype = 'ZH'

  MATC % v        = MATA % v  - MATB % v                 

 end function       matrix_subtract_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_subtract_dmatrix                                (MATA,MATB)              result(MATC)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  class (DMatrix), intent(in)           :: MATB
  type   (Matrix)                       :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_subtract_matrix: MATA and MATB not conformable')

  MATC            = MATA                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'ZG'               

  if(MATA % mtype == 'ZH' .and. MATB % mtype == 'DS') MATC % mtype = 'ZH'

  MATC % v        = MATA % v  - MATB % v                 

 end function       matrix_subtract_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_subtract_matrix                                (MATA,MATB)              result(MATC)                    !:.:.:
  class (DMatrix), intent(in)           :: MATA
  class  (Matrix), intent(in)           :: MATB
  type   (Matrix)                       :: MATC

  if(MATA % m /= MATB % m .or. MATA % n /= MATB % n)call matrix_error('matrix_subtract_matrix: MATA and MATB not conformable')

  MATC            = MATB                ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATC % mtype    =  'ZG'               

  if(MATB % mtype == 'ZH' .and. MATA % mtype == 'DS') MATC % mtype = 'ZH'

  MATC % v        = MATA % v  - MATB % v                 

 end function       dmatrix_subtract_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_minus_matrix                             (MATA)                   result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  type   (Matrix)                       :: MATB

  MATB = MATA

  MATB % v = -MATB % v

 end function       matrix_return_minus_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_mult_matrix                                       (r,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v * r                 

 end function       real_mult_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_real                                       (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v * r                 

 end function       matrix_mult_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_mult_matrix                                    (r,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 
  MATB % v = MATA % v * r                 

  if( AIMAG(r) /= ZERO ) MATB % mtype = 'ZG'                                                                                        !not Hermitian anymore

 end function       complex_mult_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_complex                                    (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v * r      

  if( AIMAG(r) /= ZERO ) MATB % mtype = 'ZG'                                                                                        !not Hermitian anymore           

 end function       matrix_mult_complex
!-----------------------------------------------------------------------------------------------------------------------------------
! Errors for incorect matrix shapes are reported by MATMUL
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_array2                                     (MATA,C)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: C(:,:)
  type   (Matrix)                       :: MATB

  
  MATB % m     = MATA % m      ; MATB % n    = size  (C,2)
  MATB % is    = MATA % is     ; MATB % ie   = MATA % ie
  MATB % js    = lbound(C,2)   ; MATB % je   = ubound(C,2)

  MATB % mtype = 'ZG'          ; MATB % name = MATA % name                                                                          ! not Hermitian in general

  allocate(MATB % v (MATB % is : MATB % ie   , MATB % js : MATB % je))

  MATB % v = LMATMUL(MATA % v, C, mtype= MATA % mtype) ! in LMATMUL, the symmetry of the first matrix determines if zhemm will be used instead of zgemm                 

 end function       matrix_mult_array2
!-----------------------------------------------------------------------------------------------------------------------------------
! Errors for incorect matrix shapes are reposrted by MATMUL
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_mult_matrix                                     (C,MATA)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: C(:,:)
  type   (Matrix)                       :: MATB

  MATB % m     = size  (C,1)   ; MATB % n    = MATA % n 
  MATB % is    = lbound(C,1)   ; MATB % ie   = ubound(C,1)
  MATB % js    = MATA % js     ; MATB % je   = MATA % je

  MATB % mtype = 'ZG'          ; MATB % name = MATA % name

  allocate(MATB % v (MATB % is : MATB % ie   , MATB % js : MATB % je))

  MATB % v = LMATMUL(C, MATA % v)                                        ! in LMATMUL, the symmetry of the first matrix determines if zhemm will be used instead of zgemm              

 end function       array2_mult_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! Errors for incorect matrix shapes are reposrted by MATMUL
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_matrix                                     (MATA,MATB)              result(MATC)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA,MATB
  type   (Matrix)                       :: MATC

  MATC % m  = MATA % m         ; MATC % n     = MATB % n
  MATC % is = MATA % is        ; MATC % ie    = MATA % ie
  MATC % js = MATB % js        ; MATC % je    = MATB % je

  MATC    % mtype =  'ZG'      ; MATC % name  = MATA % name 

  allocate(MATC % v (MATC % is : MATC % ie , MATC % js : MATC % je))

  MATC % v  = LMATMUL(MATA % v,   MATB % v, mtype= MATA % mtype)        ! in LMATMUL, the symmetry of the first matrix determines if zhemm will be used instead of zgemm            

 end function       matrix_mult_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! Errors for incorect matrix shapes are reposrted by MATMUL
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_mult_dmatrix                                    (MATA,MATB)              result(MATC)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  class (DMatrix), intent(in)           :: MATB
  type   (Matrix)                       :: MATC
  complex(dp)    , allocatable          :: B(:,:)

  MATC % m  = MATA % m          ; MATC % n     = MATB % n
  MATC % is = MATA % is         ; MATC % ie    = MATA % ie
  MATC % js = MATB % js         ; MATC % je    = MATB % je

  MATC    % mtype =  'ZG'       ; MATC % name  = MATA % name 

  allocate(   B      (MATB % is : MATB % ie , MATB % js : MATB % je))
  allocate(MATC % v  (MATC % is : MATC % ie , MATC % js : MATC % je))

  B         =         MATB % v

  MATC % v  = LMATMUL(MATA % v,      B      , mtype= MATA % mtype)       ! in LMATMUL, the symmetry of the first matrix determines if zhemm will be used instead of zgemm        

 end  function      matrix_mult_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
! Errors for incorect matrix shapes are reposrted by MATMUL
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dmatrix_mult_matrix                                    (MATA,MATB)              result(MATC)                    !:.:.:
  class (DMatrix), intent(in)           :: MATA
  class  (Matrix), intent(in)           :: MATB
  type   (Matrix)                       :: MATC
  complex(dp)    , allocatable          :: A(:,:)

  MATC % m  = MATA % m          ; MATC % n     = MATB % n
  MATC % is = MATA % is         ; MATC % ie    = MATA % ie
  MATC % js = MATB % js         ; MATC % je    = MATB % je

  MATC    % mtype =  'ZG'       ; MATC % name  = MATB % name 

  allocate(   A      (MATA % is : MATA % ie , MATA % js : MATA % je))
  allocate(MATC % v  (MATC % is : MATC % ie , MATC % js : MATC % je))

  A         =         MATA % v

  MATC % v  = LMATMUL(   A    ,   MATB % v)                 

 end  function      dmatrix_mult_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_divide_real                                     (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  real   (dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v / r                 

 end  function      matrix_divide_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_divide_complex                                  (MATA,r)                 result(MATB)                    !:.:.:
  class  (Matrix), intent(in)           :: MATA
  complex(dp)    , intent(in)           :: r
  type   (Matrix)                       :: MATB

  MATB     = MATA                       ! use fortran2003 feature: allocation on assignment for allocatable arrays (MATB % v is reallocated or allocated if not allocated) 

  MATB % v = MATA % v / r      

  if( AIMAG(r) /= ZERO ) MATB % mtype = 'ZG'                                                                                        !not Hermitian anymore           

 end  function      matrix_divide_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_mult_matrix_sub                                 (MATA,MATB,MATC)                                         !:.:.:
  class  (Matrix), intent(in)           :: MATA,MATB
  type   (Matrix), intent(inout)        :: MATC

  call  mmmult(MATA % v, MATB % v, MATC % v, mtype= MATA % mtype)

 end subroutine     matrix_mult_matrix_sub

!...................................................................................................................................
!.......................... Linear Algebra with LAPACK .............................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_inverse                                         (MATA)                   result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MATA
  type ( Matrix)                        :: MATB

  call metadata_copy(MATA,MATB)

  allocate(MATB % v,mold=MATA % v)

  MATB % v = inverse(MATA % v)

 end  function      matrix_inverse
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_determinant                                     (MAT)                    result(z)                       !:.:.:
  class( Matrix), intent(in)            :: MAT
  complex(dp)                           :: z

  z = determinant(MAT % v)

 end  function      matrix_determinant
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_lndet                                           (MAT)                    result(z)                       !:.:.:
  class( Matrix), intent(in)            :: MAT
  complex(dp)                           :: z(2)

  z = lndet(MAT % v)

 end  function      matrix_lndet
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_Pfaffian                                        (MAT)                    result(z)                       !:.:.:
  class( Matrix), intent(in)            :: MAT
  complex(dp)                           :: z

  z =   Pfaffian(MAT % v)

 end  function      matrix_Pfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_lnPfaffian                                      (MAT)                    result(z)                       !:.:.:
  class( Matrix), intent(in)            :: MAT
  complex(dp)                           :: z(2)

  z = lnPfaffian(MAT % v)

 end  function      matrix_lnPfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_eigenvalues                                     (MAT)                    result(vec)                     !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Vector)                        :: vec

  if( MAT  % m /= MAT  % m) call matrix_error('matrix_eigenvalues: MAT not a square matrix')

  vec     =      Vector( MAT % n, is= MAT % is)

  vec % v = eigenvalues( MAT % v, mtype = MAT % mtype)

 end  function      matrix_eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_eigenvectors                                    (MATA,vec)               result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MATA
  type ( Vector), intent(in out)        :: vec
  type ( Matrix)                        :: MATB
  complex(dp)   , allocatable           :: W(:,:)
  integer                               :: n

  if( MATA % m /= MATA % m) call matrix_error('matrix_eigenvectors: MATA not a square matrix')

  n        =      MATA % n ; allocate(W(n,n+1))

  if( vec % n /= n ) vec      =      Vector( n, is= MATA % is)

  MATB     =      MATA
  
  W        =      eigenvectors(MATA % v)
  
  vec  % v =      W(:,1)
  
  MATB % v =      W(:,2:)

  deallocate(W)

 end  function      matrix_eigenvectors


!...................................................................................................................................
!.......................... Misc Procedures.........................................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_symmetric_set                                   (MAT, uplo)                                              !:.:.:
  class(Matrix)                         :: MAT
  character(*),optional      ,intent(in):: uplo

  call symmetric_set( MAT % v , uplo)

 end subroutine     matrix_symmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         matrix_antisymmetric_set                               (MAT, uplo)                                              !:.:.:
  class(Matrix)                         :: MAT
  character(*),optional      ,intent(in):: uplo

  call antisymmetric_set( MAT % v , uplo)

 end subroutine     matrix_antisymmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_diagonal_get                                    (MAT)                    result(vec)                     !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Vector)                        :: vec
  
  vec     =  Vector ( MAT % m , is= MAT % is )
  vec % v = diagonal( MAT % v )

 end  function      matrix_diagonal_get
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_diagonal_set_from_matrix                        (MATA)                   result(MATB)                    !:.:.:
  class ( Matrix), intent(in)           :: MATA
  type  ( Matrix)                       :: MATB

  MATB     = MATA

  MATB % v = diagonalMatrix( MATA % v )

 end  function      matrix_diagonal_set_from_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_diagonal_set_from_vector                        (vec)                    result(MAT)                     !:.:.:
  class ( Vector), intent(in)           :: vec
  type  ( Matrix)                       :: MAT

  MAT     =  Matrix ( vec % n , is= vec % is)

  MAT % v = diagonal( vec % v )

 end  function      matrix_diagonal_set_from_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_traceless_get                                   (MAT)                    result(MATB)                    !:.:.:
  class ( Matrix), intent(in)           :: MAT
  type  ( Matrix)                       :: MATB

  call metadata_copy(MAT,MATB);  allocate(MATB % v, mold = MAT % v)

  MATB % v =  traceless(MAT % v)

 end function       matrix_traceless_get
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_traceless_set                                   (MAT)                                                    !:.:.:
  class ( Matrix), intent(in out)       :: MAT

  call traceless_set(MAT % v)

 end subroutine     matrix_traceless_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_trace2c                                         (MAT)                    result(r)                       !:.:.:
  type  ( Matrix), intent(in)           :: MAT
  complex(dp)                           :: r

  r = trace2c(MAT % v, mtype= MAT % mtype)

 end  function      matrix_trace2c
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_trace2                                          (MAT)                    result(r)                       !:.:.:
  type  ( Matrix), intent(in)           :: MAT
  complex(dp)                           :: r

  r = trace2(MAT % v, mtype= MAT % mtype)

 end  function      matrix_trace2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_trace                                           (MAT)                    result(r)                       !:.:.:
  type  ( Matrix), intent(in)           :: MAT
  complex(dp)                           :: r

  r = trace(MAT % v)

 end  function      matrix_trace
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_norm                                            (MAT)                    result(r)                       !:.:.:
  type   (Matrix), intent(in)           :: MAT
  real(dp)                              :: r
  
  r = array2_norm(MAT % v)

 end function matrix_norm
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_matrix                                   (MAT)                                                    !:.:.:
  type   (Matrix)                       :: MAT
  
  call matrix_random_set(MAT)

 end subroutine     random_number_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_matrix_gaussian                          (MAT,sigma)                                              !:.:.:
  type   (Matrix)                       :: MAT
  real   (dp)                           :: sigma
 
  call matrix_gaussian_set(MAT,sigma=sigma)

 end subroutine     random_number_matrix_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_real_dmatrix                             (MATA)                  result(MATB)                     !:.:.:
  class(Matrix), intent(in)             :: MATA
  type(DMatrix)                         :: MATB

  call matrix_metadata_copy(MATA,MATB)
  
  MATB % mtype = 'DG'

  if(MATA % mtype == 'ZH') MATB % mtype = 'DS'

  allocate(MATB % v( MATA % is : MATA % ie , MATA % js : MATA % je ))

  MATB % v = REAL(MATA % v, kind=dp)

 end function       matrix_return_real_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_imag_dmatrix                             (MATA)                   result(MATB)                    !:.:.:
  class(Matrix), intent(in)             :: MATA
  type(DMatrix)                         :: MATB

  call matrix_metadata_copy(MATA,MATB)
  
  MATB % mtype = 'DG'

  allocate(MATB % v( MATA % is : MATA % ie , MATA % js : MATA % je ))

  MATB % v = AIMAG(MATA % v)

 end function       matrix_return_imag_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_conjg                                    (MATA)                   result(MATB)                    !:.:.:
  class(Matrix), intent(in)             :: MATA
  type (Matrix)                         :: MATB

  call matrix_metadata_copy(MATA,MATB)

  allocate(MATB % v( MATA % is : MATA % ie , MATA % js : MATA % je ))

  MATB % v = CONJG(MATA % v)

 end function       matrix_return_conjg
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_transpose                                (MATA)                   result(MATB)                    !:.:.:
  class(Matrix), intent(in)             :: MATA
  type (Matrix)                         :: MATB

  !Careful m <-> n
  call matrix_metadata_put(MATB,MATA % n,MATA % m, MATA % js, MATA % is, MATA % mtype, MATA % name)

  allocate(MATB % v( MATA % js : MATA % je , MATA % is : MATA % ie ))

  MATB % v = TRANSPOSE(MATA % v)

 end function       matrix_return_transpose
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_return_hermitian                                (MATA)                   result(MATB)                    !:.:.:
  class(Matrix), intent(in)             :: MATA
  type (Matrix)                         :: MATB

  !Careful m <-> n
  call matrix_metadata_put(MATB,MATA % n,MATA % m, MATA % js, MATA % is, MATA % mtype, MATA % name)

  allocate(MATB % v( MATA % js : MATA % je , MATA % is : MATA % ie ))

  MATB % v = TRANSPOSE(CONJG(MATA % v))

 end function       matrix_return_hermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_is_hermitian                                    (MATA)                   result(r)                       !:.:.:
  class(Matrix), intent(in)             :: MATA
  real (dp)                             :: r

  r = isHermitian( MATA % v )

 end function       matrix_is_hermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_is_symmetric                                    (MATA)                   result(r)                       !:.:.:
  class(Matrix), intent(in)             :: MATA
  real (dp)                             :: r

  r = isSymmetric( MATA % v )

 end function       matrix_is_symmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_is_antisymmetric                                (MATA)                   result(r)                       !:.:.:
  class(Matrix), intent(in)             :: MATA
  real (dp)                             :: r

  r = isAntiSymmetric( MATA % v )

 end function       matrix_is_antisymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_is_nan                                          (MAT)                    result(itis)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  logical                               :: itis

  itis = .false.

  if( isNaN(MAT % v) ) itis = .true.

 end  function      matrix_is_nan


!...................................................................................................................................
!.......................... Math and Array Procedures .............................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_abs                                             (MAT)                    result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  type (DMatrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate( MATB % v ( MATB % is : MATB % ie , MATB % js : MATB % je )   )

  if(MAT % mtype == 'ZH')then
   MATB  % mtype =  'DS'
  else
   MATB  % mtype =  'DG'
  end if

  MATB % v = abs(MAT % v)

 end  function      matrix_abs
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_sin                                             (MAT)                    result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = sin(MAT % v)

 end  function      matrix_sin
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_cos                                             (MAT)                    result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = cos(MAT % v)

 end  function      matrix_cos
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_exp                                             (MAT)                    result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % v = exp(MAT % v)

 end  function      matrix_exp
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_log                                             (MAT)                    result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % mtype = 'ZG'

  MATB % v = log(MAT % v)

 end  function      matrix_log
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_sqrt                                            (MAT)                    result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % mtype = 'ZG'

  MATB % v = sqrt(MAT % v)

 end  function      matrix_sqrt
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_power_integer                                   (MAT,n)                  result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  integer       , intent(in)            :: n
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % mtype = 'ZG'

  MATB % v = (MAT % v)**n

 end  function      matrix_power_integer
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_power_real                                      (MAT,n)                  result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  real (dp)     , intent(in)            :: n
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % mtype = 'ZG'

  MATB % v = (MAT % v)**n

 end  function      matrix_power_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      matrix_power_complex                                   (MAT,n)                  result(MATB)                    !:.:.:
  class( Matrix), intent(in)            :: MAT
  complex(dp)   , intent(in)            :: n
  type ( Matrix)                        :: MATB

  call metadata_copy(MAT, MATB);  allocate(MATB % v , mold = MAT % v)

  MATB % mtype = 'ZG'

  MATB % v = (MAT % v)**n

 end  function      matrix_power_complex





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
