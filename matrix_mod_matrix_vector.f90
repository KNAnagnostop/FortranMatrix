!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_matrix_vector.f90 .....................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................
!.......................... type/class Vector procedures ...........................................................................!:.:.:
!...................................................................................................................................

!...................................................................................................................................
!.......................... Constructors ...........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_metadata_put                                    (vec,n,is,name)                                          !:.:.:
  type( Vector),         intent(in out) :: vec
  integer,               intent(in)     :: n
  integer     , optional,intent(in)     :: is
  character(*), optional,intent(in)     :: name
  integer                               :: i
  character(name_len)                   :: nm

  nm = ''; if(present(name)) nm = trim(name)
  i  = 1 ; if(present(is  )) i  = is

  vec % n     = n
  vec % is    = i
  vec % ie    = i + n - 1
  vec % name  = nm

 end subroutine      vector_metadata_put
!-----------------------------------------------------------------------------------------------------------------------------------
!Basic constructor: allocates vec % v.  Other constructors call vector_construct_zero (n,is,name)
!Usage: vec = Vector(n,is,name), sets vec % v = 0
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_construct_zero                                  (n,is,name)              result(vec)                     !:.:.:
  type(Vector)                          :: vec
  integer,                 intent(in)   :: n
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  call vector_metadata_put(vec,n,is,name)

  allocate(vec % v(vec % is:vec % ie))

  vec % v  = ZERO

 end function       vector_construct_zero
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_construct_real                                  (n,r,is,name)            result(vec)                     !:.:.:
  type(Vector)                          :: vec
  integer,                 intent(in)   :: n
  real(dp)    ,            intent(in)   :: r
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  vec     = vector_construct_zero(n,is,name)

  vec % v = r

 end function       vector_construct_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_construct_complex                               (n,z,is,name)            result(vec)                     !:.:.:
  type(Vector)                          :: vec
  integer,                 intent(in)   :: n
  complex(dp) ,            intent(in)   :: z
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  vec     = vector_construct_zero(n,is,name)

  vec % v = z

 end function       vector_construct_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_construct_array1                                (C,is,name)              result(vec)                     !:.:.:
  type(Vector)                          :: vec
  complex(dp) ,            intent(in)   :: C(:)
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  vec     = vector_construct_zero(size(C),is,name)

  vec % v = C

 end function       vector_construct_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 function           vector_construct_random                                (rtype,n,is,name,sigma)  result(vec)                     !:.:.:
  type(Vector)                          :: vec
  character(*)          ,  intent(in)   :: rtype  ! random, uniform, gaussian
  integer               ,  intent(in)   :: n
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name
  real    (dp), optional,  intent(in)   :: sigma
  real    (dp)                          :: s

  s   = 1.0_dp; if(present(sigma)) s  = sigma

  vec     = vector_construct_zero(n,is,name)

  select case (rtype(1:1))
  case('r','R', 'u', 'U')
   call random_number  (vec % v)                                                                                                   
  case('g','G')
   call random_number  (vec % v, sigma=s)                                                                                           
  case default
   call matrix_error   ('vector_construct_random: wrong rtype: '//trim(rtype))
  end select


 end function       vector_construct_random
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... Procedures..............................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 function           vector_sort                                            (veca,by)                result(vecb)                    !:.:.:
  class( Vector),         intent(in)    :: veca
  character(*), optional, intent(in)    :: by
  type ( Vector)                        :: vecb

  vecb = veca

  vecb % v = sort(vecb % v, by=by)

 end  function      vector_sort
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_dot_product                                     (veca,vecb)              result(r)                       !:.:.:
  class( Vector), intent(in)            :: veca, vecb
  complex(dp)                           :: r

  r = dot_product(veca % v,vecb % v)

 end  function      vector_dot_product
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         vector_random_set                                      (vec)                                                    !:.:.:
  class(Vector)         , intent(in out):: vec

  if( .not. allocated(vec % v)) return
  call random_number (vec % v)

 end subroutine     vector_random_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         vector_gaussian_set                                    (vec,sigma)                                              !:.:.:
  class(Vector)         , intent(in out):: vec
  real(dp)    , optional                :: sigma
  real(dp)                              :: s

  s = 1.0_dp;if(present(sigma)) s = sigma

  if( .not. allocated(vec % v)) return
  call random_number (vec % v, sigma=s)

 end subroutine     vector_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
! same as vector_gaussian_set, but sigma not optional: in order to dissambiguate the random_number interface
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_vector_gaussian_set                      (vec,sigma)                                              !:.:.:
  class(Vector)         , intent(in out):: vec
  real(dp)                              :: sigma

  if( .not. allocated(vec % v)) return
  call random_number (vec % v, sigma=sigma)

 end subroutine     random_number_vector_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_metadata_copy_vector                            (veca,vecb)                                              !:.:.:
  class( Vector), intent(in)            :: veca
  class( Vector), intent(in out)        :: vecb

  call vector_metadata_put(vecb, veca % n, veca % is, veca % name)

 end subroutine     vector_metadata_copy_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_metadata_copy_dvector                           (veca,vecb)                                              !:.:.:
  class( Vector), intent(in)            :: veca
  class(DVector), intent(in out)        :: vecb

  call dvector_metadata_put(vecb, veca % n, veca % is, veca % name)

 end subroutine     vector_metadata_copy_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_metadata_copy_vector                           (veca,vecb)                                              !:.:.:
  class(DVector), intent(in)            :: veca
  class( Vector), intent(in out)        :: vecb

  call vector_metadata_put(vecb, veca % n, veca % is, veca % name)

 end subroutine     dvector_metadata_copy_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_metadata_copy_dvector                          (veca,vecb)                                              !:.:.:
  class(DVector), intent(in)            :: veca
  class(DVector), intent(in out)        :: vecb

  call dvector_metadata_put(vecb, veca % n, veca % is, veca % name)

 end subroutine     dvector_metadata_copy_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_norm                                            (vec)                    result(r)                       !:.:.:
  class(Vector), intent(in)             :: vec
  real(dp)                              :: r

  r = norm(vec % v)

 end function       vector_norm
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_return_real_dvector                             (vec)                    result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: vec 
  type(DVector)                         :: vecb

  call vector_metadata_copy_dvector(vec ,vecb)

  allocate( vecb % v(vec  % is : vec  % ie) )

  vecb % v = REAL(vec  % v, kind = dp)
 
 end function       vector_return_real_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_return_imag_dvector                             (vec)                    result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: vec 
  type(DVector)                         :: vecb

  call vector_metadata_copy_dvector(vec ,vecb)
  
  allocate( vecb % v(vec  % is : vec  % ie) )

  vecb % v = AIMAG(vec  % v)

 end function       vector_return_imag_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_return_return_conjg                             (vec)                    result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: vec 
  type (Vector)                         :: vecb

  call vector_metadata_copy_vector(vec ,vecb)
  
  allocate( vecb % v(vec  % is : vec  % ie) )

  vecb % v =CONJG(vec  % v)

 end function       vector_return_return_conjg
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... Operators    ...........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_assignFrom_vector                               (vecb,veca)                                              !:.:.: 
  type( Vector), intent(in out)         :: vecb
  type( Vector), intent(in    )         :: veca

  vecb = veca

 end  subroutine    vector_assignFrom_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_assignFrom_real                                 (vecb,r)                                                 !:.:.: 
  type( Vector), intent(in out)         :: vecb
  real(dp)     , intent(in    )         :: r

  vecb % v = r

 end  subroutine    vector_assignFrom_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_assignFrom_complex                              (vecb,r)                                                 !:.:.: 
  type( Vector), intent(in out)         :: vecb
  complex(dp)  , intent(in    )         :: r

  vecb % v = r

 end  subroutine    vector_assignFrom_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_assignFrom_array1                               (vecb,v)                                                 !:.:.: 
  type( Vector), intent(in out)         :: vecb
  complex(dp)  , intent(in)             :: v(:)

  if(vecb % n == 0)then  ! it has not been constructed
   call vector_metadata_put(vecb,size(v,1))
   vecb % v = v
  else
   if( vecb % n   /=   size(v)   ) then
       vecb % v    =   NaN()
     return
   end if
   if( .not. allocated(vecb % v) ) allocate( vecb % v(vecb % is : vecb % ie) )
   vecb                                           % v(vecb % is : vecb % ie)    = v
  end if
 end subroutine     vector_assignFrom_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_assignFrom_array1_d                             (vecb,v)                                                 !:.:.: 
  type( Vector), intent(in out)         :: vecb
  real   (dp)  , intent(in)             :: v(:)

  if(vecb % n == 0)then  ! it has not been constructed
   call vector_metadata_put(vecb,size(v,1))
   vecb % v = v
  else
   if( vecb % n   /=   size(v)   ) then
       vecb % v    =   NaN()
     return
   end if
   if( .not. allocated(vecb % v) ) allocate( vecb % v(vecb % is : vecb % ie) )
   vecb                                           % v(vecb % is : vecb % ie)    = v
  end if
 end subroutine     vector_assignFrom_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_assignFrom_dvector                              (vecb,veca)                                              !:.:.: 
  type( Vector), intent(in out)         :: vecb
  type(DVector), intent(in)             :: veca

  if(vecb % n == 0)then  ! it has not been constructed
   call dvector_metadata_copy_vector(veca,vecb)
   vecb % v = veca % v
  else
   if( vecb % n   /=   veca % n  ) then
       vecb % v    =   NaN()
    return
   end if
   if( .not. allocated(vecb % v) ) allocate( vecb % v(vecb % is : vecb % ie) )
   vecb                                           % v(vecb % is : vecb % ie)    = veca % v
  end if
 end subroutine     vector_assignFrom_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_plus_vector                                       (r,veca)                 result(vecb)                    !:.:.: 
  class(Vector), intent(in)             :: veca
  real (dp)    , intent(in)             :: r
  type (Vector)                         :: vecb

  vecb     = veca

  vecb % v = vecb % v + r

 end function       real_plus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_plus_real                                       (veca,r)                 result(vecb)                    !:.:.: 
  class(Vector), intent(in)             :: veca
  real (dp)    , intent(in)             :: r
  type (Vector)                         :: vecb

  vecb     = veca

  vecb % v = vecb % v + r

 end function       vector_plus_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_plus_vector                                    (z,veca)                 result(vecb)                    !:.:.: 
  class(Vector), intent(in)             :: veca
  complex(dp)  , intent(in)             :: z
  type (Vector)                         :: vecb

  vecb     = veca

  vecb % v = vecb % v + z

 end function       complex_plus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_plus_complex                                    (veca,z)                 result(vecb)                    !:.:.: 
  class(Vector), intent(in)             :: veca
  complex(dp)  , intent(in)             :: z
  type (Vector)                         :: vecb

  vecb     = veca

  vecb % v = vecb % v + z

 end function       vector_plus_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_plus_array1                                     (veca,v)                 result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: veca
  complex(dp)  , intent(in)             :: v(:)
  type (Vector)                         :: vecb

  vecb     = veca
  
  vecb % v = vecb % v + v

  if( veca % n /= size(v)) vecb % v = NaN()

 end function       vector_plus_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_plus_vector                                     (v,veca)                 result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: veca
  complex(dp)  , intent(in)             :: v(:)
  type (Vector)                         :: vecb

  vecb     = veca
  
  vecb % v = vecb % v + v

  if( veca % n /= size(v)) vecb % v = NaN()

 end function      array1_plus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_plus_array1_d                                   (veca,v)                 result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: veca
  real   (dp)  , intent(in)             :: v(:)
  type (Vector)                         :: vecb

  vecb     = veca
  
  vecb % v = vecb % v + v

  if( veca % n /= size(v)) vecb % v = NaN()

 end function       vector_plus_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_d_plus_vector                                   (v,veca)                 result(vecb)                    !:.:.:
  class(Vector), intent(in)             :: veca
  real   (dp)  , intent(in)             :: v(:)
  type (Vector)                         :: vecb

  vecb     = veca
  
  vecb % v = vecb % v + v

  if( veca % n /= size(v)) vecb % v = NaN()

 end function       array1_d_plus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_plus_vector                                     (veca,vecb)              result(vecc)                    !:.:.:
  class( Vector), intent(in)            :: veca
  class( Vector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc         = veca

  vecc % v     = vecc % v + vecb % v

  if(veca % n /= vecb % n)  vecc % v = NaN()

 end function       vector_plus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_plus_dvector                                    (veca,vecb)              result(vecc)                    !:.:.:
  class( Vector), intent(in)            :: veca
  class(DVector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc         = veca

  vecc % v     = vecc % v + vecb % v

  if(veca % n /= vecb % n)  vecc % v = NaN()

 end function       vector_plus_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_plus_vector                                    (veca,vecb)              result(vecc)                    !:.:.:
  class(DVector), intent(in)            :: veca
  class( Vector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc         = vecb

  vecc % v     = vecc % v + vecb % v

  if(veca % n /= vecb % n)  vecc % v = NaN()

 end function       dvector_plus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_subtract_vector                                   (r,veca)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real(dp)      , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb % v = -vecb % v + r

 end function       real_subtract_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_subtract_real                                   (veca,r)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real(dp)      , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb % v =  vecb % v - r

 end function       vector_subtract_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_subtract_vector                                (r,veca)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb % v = -vecb % v + r

 end function       complex_subtract_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_subtract_complex                                (veca,r)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb % v =  vecb % v - r

 end function       vector_subtract_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_subtract_array1                                 (veca,C)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: C(:)
  type ( Vector)                        :: vecb

  vecb          =          veca

  vecb % v      =          vecb % v - C

  if(vecb % n /= size(C) ) vecb % v = NaN()

 end function       vector_subtract_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_subtract_vector                                 (C,veca)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: C(:)
  type ( Vector)                        :: vecb

  vecb          =          veca

  vecb % v      =         -vecb % v + C

  if(vecb % n /= size(C) ) vecb % v = NaN()

 end function       array1_subtract_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_subtract_array1_d                               (veca,C)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real   (dp)   , intent(in)            :: C(:)
  type ( Vector)                        :: vecb

  vecb          =          veca

  vecb % v      =          vecb % v - C

  if(vecb % n /= size(C) ) vecb % v = NaN()

 end function       vector_subtract_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_d_subtract_vector                               (C,veca)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real   (dp)   , intent(in)            :: C(:)
  type ( Vector)                        :: vecb

  vecb          =          veca

  vecb % v      =         -vecb % v + C

  if(vecb % n /= size(C) ) vecb % v = NaN()

 end function       array1_d_subtract_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_subtract_vector                                 (veca,vecb)              result(vecc)                    !:.:.:
  class( Vector), intent(in)            :: veca
  class( Vector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc          =           veca

  vecc % v      =           veca % v - vecb % v

  if(vecb % n /= veca % n ) vecc % v = NaN()

 end function       vector_subtract_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_subtract_dvector                                (veca,vecb)              result(vecc)                    !:.:.:
  class( Vector), intent(in)            :: veca
  class(DVector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc          =           veca

  vecc % v      =           veca % v - vecb % v

  if(vecb % n /= veca % n ) vecc % v = NaN()

 end  function      vector_subtract_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_subtract_vector                                (veca,vecb)              result(vecc)                    !:.:.:
  class(DVector), intent(in)            :: veca
  class( Vector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc          =           vecb

  vecc % v      =           veca % v - vecb % v

  if(vecb % n /= veca % n ) vecc % v = NaN()

 end  function      dvector_subtract_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_return_minus_vector                             (veca)                   result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  type ( Vector)                        :: vecb

  vecb    = veca

  vecb %v = - veca % v

 end  function      vector_return_minus_vector
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_mult_vector                                       (r,veca)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real (dp)     , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb = vecb % v * r

 end  function      real_mult_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_mult_real                                       (veca,r)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real (dp)     , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb = vecb % v * r

 end  function      vector_mult_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      complex_mult_vector                                    (r,veca)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb = vecb % v * r

 end  function      complex_mult_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_mult_complex                                    (veca,r)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb = vecb % v * r

 end  function      vector_mult_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 function           matrix_mult_vector                                     (MATA,vecb)              result(vecc)                    !:.:.:
  class( Matrix), intent(in)            :: MATA
  class( Vector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc % n      = MATA % m
  vecc % is     = MATA % is
  vecc % ie     = MATA % ie
  vecc % name   = vecb % name

  allocate( vecc % v(vecc % is : vecc % ie) )

  if( MATA % n /= vecb % n)then
   vecc % v     = NaN()
   return
  end if

  if( MATA % mtype(2:2) == 'H') then
   vecc % v     = LMATMUL(MATA % v , vecb % v, type='H')
  else
   vecc % v     = LMATMUL(MATA % v , vecb % v)
  end if

 end  function      matrix_mult_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 function           vector_mult_matrix                                     (vecb,MATA)              result(vecc)                    !:.:.:
  class( Matrix), intent(in)            :: MATA
  class( Vector), intent(in)            :: vecb
  type ( Vector)                        :: vecc

  vecc % n      = MATA % n
  vecc % is     = MATA % js
  vecc % ie     = MATA % je
  vecc % name   = vecb % name

  allocate( vecc % v(vecc % is : vecc % ie) )

  if( MATA % m /= vecb % n)then
   vecc % v     = NaN()
   return
  end if

  if( MATA % mtype(2:2) == 'H') then
   vecc % v     = LMATMUL(CONJG(MATA % v) , vecb % v, type='H')
  else
   vecc % v     = LMATMUL(      MATA % v  , vecb % v, type='T')
  end if

 end  function      vector_mult_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_divide_real                                     (veca,r)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  real   (dp)   , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb = vecb % v / r

 end  function      vector_divide_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_divide_complex                                  (veca,r)                 result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: veca
  complex(dp)   , intent(in)            :: r
  type ( Vector)                        :: vecb

  vecb = veca

  vecb = vecb % v / r

 end  function      vector_divide_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_mult_vector_sub                                (MATA,vecb,vecc)                                          !:.:.:
  class  (Matrix), intent(in)           :: MATA
  class  (Vector), intent(in)           :: vecb
  type   (Vector), intent(inout)        :: vecc

  call  mvmult(MATA % v, vecb % v, vecc % v)

 end subroutine     matrix_mult_vector_sub
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_mult_matrix_sub                                (vecb,MATA,vecc)                                          !:.:.:
  class  (Matrix), intent(in)           :: MATA
  class  (Vector), intent(in)           :: vecb
  type   (Vector), intent(inout)        :: vecc

  call  vmmult(vecb % v, MATA % v, vecc % v)

 end subroutine     vector_mult_matrix_sub


!...................................................................................................................................
!.......................... Math and Array Procedures .............................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_abs                                             (vec)                    result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate( vecb % v ( vecb % is : vecb % ie)   )

  vecb % v = abs(vec % v)

 end  function      vector_abs
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_sin                                             (vec)                    result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = sin(vec % v)

 end  function      vector_sin
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_cos                                             (vec)                    result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = cos(vec % v)

 end  function      vector_cos
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_exp                                             (vec)                    result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = exp(vec % v)

 end  function      vector_exp
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_log                                             (vec)                    result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = log(vec % v)

 end  function      vector_log
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_sqrt                                            (vec)                    result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = sqrt(vec % v)

 end  function      vector_sqrt
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_power_integer                                   (vec,n)                  result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  integer       , intent(in)            :: n
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = (vec % v)**n

 end  function      vector_power_integer
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_power_real                                      (vec,n)                  result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  real (dp)     , intent(in)            :: n
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = (vec % v)**n

 end  function      vector_power_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vector_power_complex                                   (vec,n)                  result(vecb)                    !:.:.:
  class( Vector), intent(in)            :: vec
  complex(dp)   , intent(in)            :: n
  type ( Vector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = (vec % v)**n

 end  function      vector_power_complex


!...................................................................................................................................
!.......................... Utilities ..............................................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
 function           vector_is_nan                                          (vec)                    result(itis)                    !:.:.:
  class( Vector), intent(in)            :: vec
  logical                               :: itis

  itis = .false.

  if( isNaN(vec % v) ) itis = .true.

 end  function      vector_is_nan
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         vector_read                                            (vec,unit)                                               !:.:.:
  class(Vector)         , intent(in out):: vec
  integer     , optional, intent(in)    :: unit
  integer                               :: n,m,i,is,ie,un
  real(dp)                              :: x,y

  un    = f_minput; if(present(unit)) un     = unit

  if( .not. allocated(vec % v) ) return

  n  = vec % n
  is = vec % is; ie = vec % ie

  do i = is, ie
   read (un,*) x, y
   vec % v(i)       = CMPLX(x,y,kind=dp)
  end do

  
 end subroutine     vector_read
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         vector_save                                            (vec,unit,fmt)                                           !:.:.:
  class(Vector)         , intent(in)    :: vec
  integer     , optional, intent(in)    :: unit
  character(*), optional, intent(in)    :: fmt
  character(20)                         :: fmtout
  integer                               :: n,m,i,is,ie,un

  un    = f_mout ;if(present(unit)) un     = unit
  fmtout='G28.17';if(present(fmt )) fmtout = fmt
  
  
  if( .not. allocated(vec % v) ) return

  n  = vec % n
  is = vec % is; ie = vec % ie

  do i = is, ie
   write(un,'(2'//fmtout//')') vec % v(i)
  end do

  
 end subroutine     vector_save
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         vector_print                                           (vec,unit,fmt,form,name,ips,ipe)                         !:.:.:
  class(Vector),           intent(in)   :: vec
  integer,       optional, intent(in)   :: unit
  character(*),  optional, intent(in)   :: fmt
  character(*),  optional, intent(in)   :: form
  integer     ,  optional, intent(in)   :: ips,ipe
  character(*),  optional, intent(in)   :: name
  character(name_len)                   :: namem
  character(25)                         :: fmts,fmtn,myform,eol
  character(math_len)                   :: nleft,nright
  integer                               :: n,m,i,is,ie,un


  if( vec % n == 0              ) call matrix_error('vector_print: vec not initialized: vec % n = 0')
  if( .not. allocated(vec % v)  ) call matrix_error('vector_print: vec % v not allocated')

  un    = f_mout     ; if(present(unit)) un      = unit
  fmtn  = 'G28.17'   ; if(present(fmt )) fmtn    = fmt

  m     = vec % n    ;
  is    = vec % is   ; ie      = vec % ie

  namem = vec % name ; if(present(name)) namem   = name

  if(present(ips) .and. ips >= is) is = ips  ! if asked, print only the ips:ipe, jps:jpe section of the arrays
  if(present(ipe) .and. ipe <= ie) ie = ipe

  myform = 'default'
  if(present(form)) myform = form

  select case (myform(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('d','D')
   write(un,'(A,I6,A,I6,A3,I6,A2)') &
        '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
        '#pvstart: ' // TRIM(namem) // ' =[                    (is a ',m,'  vector) (printing bounds: ',    &
        is,' : ',ie,' )'

   write(fmts,'(A)') '(2'//TRIM(fmtn)//')'
   do i = is, ie
    write(un,fmts) vec % v(i) 
   end do;

   write(un,'(A           )') &
       '#pvend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I')
   write(un,'(A,I12,A,I12,A,I12,A,I12,A,I12)') &
        '#info: vector ' // TRIM(namem) // ' is a ',m,'  vector,  bounds: ',    &
        is,' : ',ie,' :: v(is:ie) :: ', lbound(vec % v),' : ',ubound(vec % v)
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('m','M') ! Mathematica Format
   nleft=Mathematica_nleft; nright=Mathematica_nright ! from matrix_mod_common, a mathematica number in exponential notation must be between nleft and nright
   write(un,'(A)') trim(namem) // '={'

   !--------------------------------------------------------------------------------------------------------------------------------
   do  i = is, ie

    eol = '),'
    if(i == ie) eol = ')'
    write(un,'(A,'//trim(fmtn)//',A,'//trim(fmtn)//',A)') &
          trim(nleft),Real(vec % v(i),kind=dp),trim(nright)//' + I ('//trim(nleft),AIMAG(vec % v(i)),trim(nright)//trim(eol)

   end do !i = is, ie
   !--------------------------------------------------------------------------------------------------------------------------------
   write(un,'(A)') '};'
   !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call matrix_error('vector_print: unknown form for printing vec')
  !---------------------------------------------------------------------------------------------------------------------------------
  end select

 end subroutine vector_print


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
