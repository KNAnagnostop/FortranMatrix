!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_matrix_dvector.f90 ....................................................................!:.:.:
!...................................................................................................................................
!...................................................................................................................................
!.......................... type/class DVector procedures ..........................................................................!:.:.:
!...................................................................................................................................

!...................................................................................................................................
!.......................... Constructors ...........................................................................................!:.:.:
!...................................................................................................................................



!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_metadata_put                                   (vec,n,is,name)                                          !:.:.:
  type(DVector),         intent(in out) :: vec
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

 end subroutine     dvector_metadata_put
!-----------------------------------------------------------------------------------------------------------------------------------
!Basic constructor: allocates vec % v.  Other constructors call dvector_construct_zero (n,is,name)
!Usage: vec = Vector(n,is,name), sets vec % v = 0
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_construct_zero                                 (n,is,name)              result(vec)                     !:.:.:
  type(DVector)                         :: vec
  integer,                 intent(in)   :: n
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  call dvector_metadata_put(vec,n,is,name)

  allocate(vec % v(vec % is:vec % ie))

  vec % v  = ZERO

 end function       dvector_construct_zero
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_construct_real                                 (n,r,is,name)            result(vec)                     !:.:.:
  type(DVector)                         :: vec
  integer,                 intent(in)   :: n
  real(dp)    ,            intent(in)   :: r
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  vec     = dvector_construct_zero(n,is,name)

  vec % v = r

 end function       dvector_construct_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_construct_complex                              (n,z,is,name)            result(vec)                     !:.:.:
  type(DVector)                         :: vec
  integer,                 intent(in)   :: n
  complex(dp) ,            intent(in)   :: z
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  vec     = dvector_construct_zero(n,is,name)

  vec % v = Real(z,kind=dp)

 end function       dvector_construct_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_construct_array1                               (C,is,name)              result(vec)                     !:.:.:
  type(DVector)                         :: vec
  real   (dp) ,            intent(in)   :: C(:)
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name

  vec     = dvector_construct_zero(size(C),is,name)

  vec % v = C

 end function       dvector_construct_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dvector_construct_random                               (rtype,n,is,name,sigma)  result(vec)                     !:.:.:
  type(DVector)                         :: vec
  character(*)          ,  intent(in)   :: rtype  ! random, uniform, gaussian
  integer               ,  intent(in)   :: n
  integer     , optional,  intent(in)   :: is
  character(*), optional,  intent(in)   :: name
  real    (dp), optional,  intent(in)   :: sigma
  real    (dp)                          :: s

  s       = 1.0_dp;  if(present(sigma)) s  = sigma

  vec     = dvector_construct_zero(n,is,name)

  select case (rtype(1:1))
  case('r','R', 'u', 'U')
   call random_number  (vec % v)                                                                                                   
  case('g','G')
   call random_number  (vec % v, sigma=s)                                                                                           
  case default
   call matrix_error   ('dvector_construct_random: wrong rtype: '//trim(rtype))
  end select


 end function       dvector_construct_random
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... Procedures..............................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_maxval                                         (veca,mask)              result(r)                       !:.:.:
  class(DVector), intent(in)            :: veca
  logical       , intent(in), optional  :: mask(:)
  real   (dp)                           :: r

  if(present(mask))then
   r = maxval(veca % v,mask=mask)
  else
   r = maxval(veca % v)
  end if

 end  function      dvector_maxval
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_minval                                         (veca,mask)              result(r)                       !:.:.:
  class(DVector), intent(in)            :: veca
  logical       , intent(in), optional  :: mask(:)
  real   (dp)                           :: r

  if(present(mask))then
   r = minval(veca % v,mask=mask)
  else
   r = minval(veca % v)
  end if

 end  function      dvector_minval
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_dot_product                                    (veca,vecb)              result(r)                       !:.:.:
  class(DVector), intent(in)            :: veca, vecb
  real   (dp)                           :: r

  r = dot_product(veca % v,vecb % v)

 end  function      dvector_dot_product
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         dvector_random_set                                     (vec)                                                    !:.:.:
  class(DVector)        , intent(in out):: vec

  if( .not. allocated(vec % v)) return
  call random_number (vec % v)

 end subroutine     dvector_random_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         dvector_gaussian_set                                   (vec,sigma)                                              !:.:.:
  class(DVector)        , intent(in out):: vec
  real(dp)    , optional                :: sigma
  real(dp)                              :: s

  s = 1.0_dp;if(present(sigma)) s = sigma

  if( .not. allocated(vec % v)) return
  call random_number (vec % v, sigma=s)

 end subroutine     dvector_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
! same as dvector_gaussian_set, but sigma not optional: in order to dissambiguate the random_number interface
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_dvector_gaussian_set                     (vec,sigma)                                              !:.:.:
  class(DVector)        , intent(in out):: vec
  real(dp)                              :: sigma

  if( .not. allocated(vec % v)) return
  call random_number (vec % v, sigma=sigma)

 end subroutine     random_number_dvector_gaussian_set
!-----------------------------------------------------------------------------------------------------------------------------------
!In matrix_mod_matrix_vector.f90 we have the following routines:
!pure subroutine    vector_metadata_copy_vector (veca,vecb)
!                   vector_metadata_copy_dvector(veca,vecb)
!                  dvector_metadata_copy_vector (veca,vecb)
!                  dvector_metadata_copy_dvector(veca,vecb)
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_norm                                           (vec)                    result(r)                       !:.:.:
  class(DVector), intent(in)            :: vec
  real(dp)                              :: r

  r = norm(vec % v)

 end function       dvector_norm
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dvector_sort                                           (veca,by)                result(vecb)                    !:.:.:
  class(DVector),         intent(in)    :: veca
  character(*), optional, intent(in)    :: by
  type (DVector)                        :: vecb

  vecb = veca

  vecb % v = sort(vecb % v, by=by)

 end  function      dvector_sort
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dvector_dvector_complex_return_vector                  (veca,vecb)              result(vecc)                    !:.:.:
  class(DVector), intent(in)            :: veca,vecb
  type ( Vector)                        :: vecc

  if(veca % n /= vecb % n) call matrix_error('dvector_dvector_complex_return_vector: veca, vecb not conformable')

  call metadata_copy(veca, vecc)

  allocate(vecc % v( vecc % is : vecc % ie ) )

  vecc % v = CMPLX ( veca % v  , vecb % v  , kind = dp)

 end  function      dvector_dvector_complex_return_vector



!...................................................................................................................................
!.......................... Operators    ...........................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_assignFrom_vector                              (vecb,veca)                                              !:.:.: 
  type(DVector), intent(in out)         :: vecb
  type( Vector), intent(in    )         :: veca

  if(vecb % n == 0)then  ! it has not been constructed
   call vector_metadata_copy_dvector(veca,vecb)
   vecb % v = veca % v
  else
   if( vecb % n   /=   veca % n  ) then
       vecb % v    =   NaN()
    return
   end if
   if( .not. allocated(vecb % v) ) allocate( vecb % v(vecb % is : vecb % ie) )
   vecb                                           % v(vecb % is : vecb % ie)    = veca % v
  end if

 end  subroutine    dvector_assignFrom_vector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_assignFrom_real                                (vecb,r)                                                 !:.:.: 
  type(DVector), intent(in out)         :: vecb
  real(dp)     , intent(in    )         :: r

  vecb % v = r

 end  subroutine    dvector_assignFrom_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_assignFrom_complex                             (vecb,r)                                                 !:.:.: 
  type(DVector), intent(in out)         :: vecb
  complex(dp)  , intent(in    )         :: r

  vecb % v = r

 end  subroutine    dvector_assignFrom_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_assignFrom_array1                              (vecb,v)                                                 !:.:.: 
  type(DVector), intent(in out)         :: vecb
  complex(dp)  , intent(in)             :: v(:)

  if(vecb % n == 0)then  ! it has not been constructed
   call dvector_metadata_put(vecb,size(v,1))
   vecb % v = v
  else
   if( vecb % n   /=   size(v)   ) return
   if( .not. allocated(vecb % v) ) allocate( vecb % v(vecb % is : vecb % ie) )
   vecb                                           % v(vecb % is : vecb % ie)    = v
  end if
 end subroutine     dvector_assignFrom_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    dvector_assignFrom_array1_d                            (vecb,v)                                                 !:.:.: 
  type(DVector), intent(in out)         :: vecb
  real   (dp)  , intent(in)             :: v(:)

  if(vecb % n == 0)then  ! it has not been constructed
   call dvector_metadata_put(vecb,size(v,1))
   vecb % v = v
  else
   if( vecb % n   /=   size(v)   ) return
   if( .not. allocated(vecb % v) ) allocate( vecb % v(vecb % is : vecb % ie) )
   vecb                                           % v(vecb % is : vecb % ie)    = v
  end if
 end subroutine     dvector_assignFrom_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_plus_dvector                                      (r,veca)                 result(vecb)                    !:.:.: 
  class(DVector),intent(in)             :: veca
  real (dp)    , intent(in)             :: r
  type (DVector)                        :: vecb

  vecb     = veca

  vecb % v = vecb % v + r

 end function       real_plus_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_plus_real                                      (veca,r)                 result(vecb)                    !:.:.: 
  class(DVector),intent(in)             :: veca
  real (dp)    , intent(in)             :: r
  type (DVector)                        :: vecb

  vecb     = veca

  vecb % v = vecb % v + r

 end function       dvector_plus_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_plus_array1_d                                  (veca,v)                 result(vecb)                    !:.:.:
  class(DVector),intent(in)             :: veca
  real   (dp)  , intent(in)             :: v(:)
  type (DVector)                        :: vecb

  vecb     = veca
  
  vecb % v = vecb % v + v

  if( veca % n /= size(v)) vecb % v = NaN()

 end function       dvector_plus_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_d_plus_dvector                                  (v,veca)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real   (dp)  , intent(in)             :: v(:)
  type (DVector)                        :: vecb

  vecb     = veca
  
  vecb % v = vecb % v + v

  if( veca % n /= size(v)) vecb % v = NaN()

 end function       array1_d_plus_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_plus_dvector                                   (veca,vecb)              result(vecc)                    !:.:.:
  class(DVector), intent(in)            :: veca
  class(DVector), intent(in)            :: vecb
  type (DVector)                        :: vecc

  vecc         = veca

  vecc % v     = vecc % v + vecb % v

  if(veca % n /= vecb % n)  vecc % v = NaN()

 end function       dvector_plus_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_subtract_dvector                                  (r,veca)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real(dp)      , intent(in)            :: r
  type (DVector)                        :: vecb

  vecb = veca

  vecb % v = -vecb % v + r

 end function       real_subtract_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_subtract_real                                  (veca,r)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real(dp)      , intent(in)            :: r
  type (DVector)                        :: vecb

  vecb = veca

  vecb % v =  vecb % v - r

 end function       dvector_subtract_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_subtract_array1_d                              (veca,C)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real   (dp)   , intent(in)            :: C(:)
  type (DVector)                        :: vecb

  vecb          =          veca

  vecb % v      =          vecb % v - C

  if(vecb % n /= size(C) ) vecb % v = NaN()

 end function       dvector_subtract_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_d_subtract_dvector                              (C,veca)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real   (dp)   , intent(in)            :: C(:)
  type (DVector)                        :: vecb

  vecb          =          veca

  vecb % v      =         -vecb % v + C

  if(vecb % n /= size(C) ) vecb % v = NaN()

 end function       array1_d_subtract_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_subtract_dvector                               (veca,vecb)              result(vecc)                    !:.:.:
  class(DVector), intent(in)            :: veca
  class(DVector), intent(in)            :: vecb
  type (DVector)                        :: vecc

  vecc          =           veca

  vecc % v      =           veca % v - vecb % v

  if(vecb % n /= veca % n ) vecc % v = NaN()

 end function       dvector_subtract_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_return_minus_dvector                           (veca)                   result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  type (DVector)                        :: vecb

  vecb    = veca

  vecb %v = - veca % v

 end  function      dvector_return_minus_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      real_mult_dvector                                      (r,veca)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real (dp)     , intent(in)            :: r
  type (DVector)                        :: vecb

  vecb = veca

  vecb = vecb % v * r

 end  function      real_mult_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_mult_real                                      (veca,r)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real (dp)     , intent(in)            :: r
  type (DVector)                        :: vecb

  vecb = veca

  vecb = vecb % v * r

 end  function      dvector_mult_real
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dmatrix_mult_dvector                                   (MATA,vecb)              result(vecc)                    !:.:.:
  class(DMatrix), intent(in)            :: MATA
  class(DVector), intent(in)            :: vecb
  type (DVector)                        :: vecc

  vecc % n      = MATA % m
  vecc % is     = MATA % is
  vecc % ie     = MATA % ie
  vecc % name   = vecb % name

  allocate( vecc % v(vecc % is : vecc % ie) )

  if( MATA % n /= vecb % n)then
   vecc % v     = NaN()
   return
  end if

  vecc % v     = LMATMUL(MATA % v , vecb % v)

 end  function      dmatrix_mult_dvector
!-----------------------------------------------------------------------------------------------------------------------------------
 function           dvector_mult_dmatrix                                   (vecb,MATA)              result(vecc)                    !:.:.:
  class(DMatrix), intent(in)            :: MATA
  class(DVector), intent(in)            :: vecb
  type (DVector)                        :: vecc

  vecc % n      = MATA % n
  vecc % is     = MATA % js
  vecc % ie     = MATA % je
  vecc % name   = vecb % name

  allocate( vecc % v(vecc % is : vecc % ie) )

  if( MATA % m /= vecb % n)then
   vecc % v     = NaN()
   return
  end if

  vecc % v     = LMATMUL(MATA % v , vecb % v, type='T')

 end  function      dvector_mult_dmatrix
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_divide_real                                    (veca,r)                 result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: veca
  real   (dp)   , intent(in)            :: r
  type (DVector)                        :: vecb

  vecb = veca

  vecb = vecb % v / r

 end  function      dvector_divide_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    matrix_mult_vector_sub_d                              (MATA,vecb,vecc)                                          !:.:.:
  class (DMatrix), intent(in)           :: MATA
  class (DVector), intent(in)           :: vecb
  type  (DVector), intent(inout)        :: vecc

  call  mvmult(MATA % v, vecb % v, vecc % v)

 end subroutine     matrix_mult_vector_sub_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    vector_mult_matrix_sub_d                              (vecb,MATA,vecc)                                          !:.:.:
  class (DMatrix), intent(in)           :: MATA
  class (DVector), intent(in)           :: vecb
  type  (DVector), intent(inout)        :: vecc

  call  vmmult(vecb % v, MATA % v, vecc % v)

 end subroutine     vector_mult_matrix_sub_d

!...................................................................................................................................
!.......................... Math and Array Procedures .............................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_abs                                            (vec)                    result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = abs(vec % v)

 end  function      dvector_abs
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_sin                                            (vec)                    result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = sin(vec % v)

 end  function      dvector_sin
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_cos                                            (vec)                    result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = cos(vec % v)

 end  function      dvector_cos
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_exp                                            (vec)                    result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = exp(vec % v)

 end  function      dvector_exp
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_log                                            (vec)                    result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = log(vec % v)

 end  function      dvector_log
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_sqrt                                           (vec)                    result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = sqrt(vec % v)

 end  function      dvector_sqrt
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_power_integer                                  (vec,n)                  result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  integer       , intent(in)            :: n
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = (vec % v)**n

 end  function      dvector_power_integer
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      dvector_power_real                                     (vec,n)                  result(vecb)                    !:.:.:
  class(DVector), intent(in)            :: vec
  real (dp)     , intent(in)            :: n
  type (DVector)                        :: vecb

  call metadata_copy(vec, vecb);  allocate(vecb % v , mold = vec % v)

  vecb % v = (vec % v)**n

 end  function      dvector_power_real


!...................................................................................................................................
!.......................... Utilities ..............................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 function           dvector_is_nan                                         (vec)                    result(itis)                    !:.:.:
  class(DVector), intent(in)            :: vec
  logical                               :: itis

  itis = .false.

  if( isNaN(vec % v) ) itis = .true.

 end  function      dvector_is_nan
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         dvector_read                                           (vec,unit)                                               !:.:.:
  class(DVector)        , intent(in out):: vec
  integer     , optional, intent(in)    :: unit
  integer                               :: n,m,i,is,ie,un

  un    = f_minput; if(present(unit)) un     = unit

  if( .not. allocated(vec % v) ) return

  n  = vec % n
  is = vec % is; ie = vec % ie

  do i = is, ie
   read (un,*) vec % v(i)
  end do

  
 end subroutine     dvector_read
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         dvector_save                                           (vec,unit,fmt)                                           !:.:.:
  class(DVector)        , intent(in)    :: vec
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
   write(un,'('//fmtout//')') vec % v(i)
  end do

  
 end subroutine     dvector_save
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         dvector_print                                          (vec,unit,fmt,form,name,ips,ipe)                         !:.:.:
  class(DVector),          intent(in)   :: vec
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
        '#pvstart: ' // TRIM(namem) // ' =[                    (is a ',m,'  dvector) (printing bounds: ',    &
        is,' : ',ie,' )'

   write(fmts,'(A)') '('//TRIM(fmtn)//')'
   do i = is, ie
    write(un,fmts) vec % v(i) 
   end do;

   write(un,'(A           )') &
       '#pvend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I')
   write(un,'(A,I12,A,I12,A,I12,A,I12,A,I12)') &
        '#info: dvector ' // TRIM(namem) // ' is a ',m,'  dvector,  bounds: ',    &
        is,' : ',ie,' :: v(is:ie) :: ', lbound(vec % v),' : ',ubound(vec % v)
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('m','M') ! Mathematica Format
   nleft=Mathematica_nleft; nright=Mathematica_nright ! from matrix_mod_common, a mathematica number in exponential notation must be between nleft and nright
   write(un,'(A)') trim(namem) // '={'

   !--------------------------------------------------------------------------------------------------------------------------------
   do  i = is, ie

    eol  = ','
    if(i == ie) eol = ''
    write(un,'(A,'//trim(fmtn)//',A)') &
          trim(nleft),vec % v(i),trim(nright)//trim(eol)

   end do !i = is, ie
   !--------------------------------------------------------------------------------------------------------------------------------
   write(un,'(A)') '};'
   !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call matrix_error('dvector_print: unknown form for printing vec')
  !---------------------------------------------------------------------------------------------------------------------------------
  end select

 end subroutine     dvector_print


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
