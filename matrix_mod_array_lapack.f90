!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_array_lapack.f90 ......................................................................!:.:.:
!...................................................................................................................................!:.:.:
!...................................................................................................................................

!...................................................................................................................................
!.......................... Matrix Matrix Multiplication ...........................................................................!:.:.:
!...................................................................................................................................
!
! Performance issues:
! https://software.intel.com/en-us/articles/a-simple-example-to-measure-the-performance-of-an-intel-mkl-function
!
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_mm                                (A,B,mtype)              result(C)                       !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)  :: A,B
  character(*) , optional                  , intent(in)  :: mtype
! complex  (dp), dimension(size(A,1),size(B,2))          :: C
  complex  (dp), dimension(:,:), allocatable             :: C ! 260104
  integer                                                :: ma,na,mb,nb,mc,nc
  complex  (dp)                                          :: alpha, beta
  character(1)                                           :: opA, opB, side,uplo
  character(mtype_len)                                   :: mt

  allocate(C(size(A,1),size(B,2)))

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)

  ma = size(A,1) ; na = size(A,2)
  mb = size(B,1) ; nb = size(B,2)
  mc = size(A,1) ; nc = size(B,2)

  C  = ZERO ! 260104

  if( na /= mb ) then
   C = NaN()
   return
  end if

  ! Lapack calculates C = alpha * A . B + beta * C
  !   C    =   A   x   B
  ! m x n    m x k   k x n
  ! call zgemm('N','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) where 'N'  is for normal, 'C' is for Hermitian conjugate
  ! If A is Hermitian, then we can call:
  ! call zhemm('L','U',M,N,  alpha,A,lda,B,ldb,beta,C,ldc) where side= 'L' is for A.B, side= 'R' for B.A and uplo='U'/'L' references the upper/lower part of the A matrix

  alpha = ONE ; beta = ZERO
  opA   = 'N' ; opB  = 'N'  ; side ='L' ; uplo = 'U'

  select case(mt(1:1))
  case('H','h')
    
!  if( ma /= na ) call matrix_error('array2_matmul_lapack_mm: ma /= na')
   if( ma /= na ) then
    C = NaN()
    return
   end if
   if( mc == 0 .or. nc == 0 ) return

   call zhemm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Hermitian               

  case default

   if( mc == 0 .or. nc == 0 .or. na == 0 ) return
    
   call zgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                          

  end select

 end function       array2_matmul_lapack_mm
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_dd                                (A,B,mtype)              result(C)                       !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A,B
  character(*) , optional                  , intent(in)       :: mtype
! real     (dp), dimension(size(A,1),size(B,2))               :: C
  real     (dp), dimension(:,:), allocatable                  :: C ! 260104
  integer                                                     :: ma,na,mb,nb,mc,nc
  real   (dp)                                                 :: alpha, beta
  character(1)                                                :: opA, opB, side,uplo
  character(mtype_len)                                        :: mt

  allocate(C(size(A,1),size(B,2)))

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)

  ma = size(A,1) ; na = size(A,2)
  mb = size(B,1) ; nb = size(B,2)
  mc = size(A,1) ; nc = size(B,2)

  C  = 0.0_dp ! 260104

  !if( na /= mb )  call matrix_error('array2_matmul_lapack_dd: na /= mb')
  if( na /= mb )then
   C = NaN(0.0_dp)
   return
  end if

  ! Lapack calculates C = alpha * A . B + beta * C
  !   C    =   A   x   B
  ! m x n    m x k   k x n
  ! call dgemm('N','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) where 'N'  is for normal, 'T' is for transpose
  ! If A is Hermitian, then we can call:
  ! call dsymm('L','U',M,N,  alpha,A,lda,B,ldb,beta,C,ldc) where side= 'L' is for A.B, side= 'R' for B.A and uplo='U'/'L' references the upper/lower part of the A matrix

  alpha = 1.0_dp ; beta = 0.0_dp
  opA   = 'N' ; opB  = 'N'  ; side ='L' ; uplo = 'U'

  select case(mt(1:1))
  case('S','s')
    
  !if( ma /= na ) call matrix_error('array2_matmul_lapack_dd: ma /= na')
   if( ma /= na )then
    C = NaN(0.0_dp)
    return
   end if
   if( mc == 0 .or. nc == 0 ) return

   call dsymm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Symmetric              

  case default

   if( mc == 0 .or. nc == 0 .or. na == 0 ) return
   
   call dgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                         

  end select

 end function       array2_matmul_lapack_dd
!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^\dagger = array2_matmul_lapack_m(A)    ;   A^\dagger . A = array2_matmul_lapack_m(A,side='H')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_m                                 (A,side)                 result(C)                       !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: A
  character(*) , optional                  , intent(in)       :: side
  complex  (dp), dimension(:,:), allocatable                  :: C
  integer                                                     :: ma, na      
  real   (dp)                                                 :: alpha, beta
  character(1)                                                :: opA, uplo
  integer                                                     :: N,K,i,j

  ma  = size(A,1) ; na = size(A,2)

  opA = 'N';                       ! default is A . A^\dagger
  N   =  ma       ; K  = na

  if(present(side))then
   if(len(side) >= 1)then
    select case(side(1:1))
    case('N','n',' ')                       ! default: A . A^\dagger (A . A^T)
    case('C','c','T','t','H','h','S','s','L','l')

     opA = 'C'
     N   = na      ; K  = ma

    case default                            ! 261124: unrecognised side -> NaN (BUG-05)
     opA = '?'
    end select
   end if
  end if

  ! Lapack calculates C = alpha * A . A^\dagger + beta * C
  ! if opA='N' then   C    =   A   x   A^\dagger          ::          if opA='T'    C    =   A^\dagger    x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m          m x n
  !
  ! call zherk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'

  alpha = 1.0_dp ; beta = 0.0_dp; uplo = 'U'
  allocate(C(N,N)) ; C = ZERO ! 260104
  if( opA == '?' )then                   ! 261124: unrecognised side (BUG-05)
   C = NaN()
   return
  end if
  if( N == 0 .or. K == 0 ) return
  call zherk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                      

  !only the upper triangle part of C is filled:
  do concurrent(i=1:n, j=1:n, i>j)
   C(i,j) = CONJG(C(j,i))
  end do

 end function       array2_matmul_lapack_m

!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^T = array2_matmul_lapack_d(A)    ;   A^T . A = array2_matmul_lapack_d(A,side='S')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_d                                 (A,side)                 result(C)                       !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A
  character(*) , optional                  , intent(in)       :: side
  real     (dp), dimension(:,:), allocatable                  :: C
  integer                                                     :: ma, na      
  real   (dp)                                                 :: alpha, beta
  character(1)                                                :: opA, uplo
  integer                                                     :: N,K,i,j

  ma  = size(A,1) ; na = size(A,2)

  opA = 'N';                       ! default is A . A^T
  N   =  ma       ; K  = na

  if(present(side))then
   if(len(side) >= 1)then
    select case(side(1:1))
    case('N','n',' ')                       ! default: A . A^\dagger (A . A^T)
    case('C','c','T','t','H','h','S','s','L','l')

     opA = 'T'
     N   = na      ; K  = ma

    case default                            ! 261124: unrecognised side -> NaN (BUG-05)
     opA = '?'
    end select
   end if
  end if


  ! Lapack calculates C = alpha * A . A^T + beta * C
  ! if opA='N' then   C    =   A   x   A^T                ::          if opA='T'    C    =     A^T   x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m     m x n
  !
  ! call dsyrk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'


  alpha = 1.0_dp ; beta = 0.0_dp; uplo = 'U'
  allocate(C(N,N)) ; C  = 0.0_dp ! 260104
  if( opA == '?' )then                   ! 261124: unrecognised side (BUG-05)
   C = NaN(0.0_dp)
   return
  end if
  if( N == 0 .or. K == 0 ) return
  call dsyrk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                      
  !only the upper triangle part of C is filled:
  do concurrent(i=1:n, j=1:n, i>j)
   C(i,j) = C(j,i)
  end do


 end function       array2_matmul_lapack_d


!...................................................................................................................................
!.......................... Matrix Vector Multiplication ...........................................................................!:.:.:
!...................................................................................................................................
!
!-----------------------------------------------------------------------------------------------------------------------------------
!type='N'    w = A . v           w_i = A_ij  v_j
!type='T'    w = A^T . v = v . A w_i = v_j   A_ji
!type='H'    w = A . v           w_i = A_ij  v_j    faster, when A = A^\dagger (no check performed)
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_mv                                (A,v,type)               result(w)                       !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: A
  complex  (dp), dimension(:)  , contiguous, intent(in)       :: v
  character(*) , optional                  , intent(in)       :: type
  complex  (dp), dimension(:)  , allocatable                  :: w
  integer                                                     :: ma,na,nv
  complex  (dp)                                               :: alpha, beta
  character(1)                                                :: tp,uplo
  integer                                                     :: incx, incy

  tp = 'N'                                        ! default is w = A.v
  if(present(type))then
   if(len(type) >= 1)then
    tp = type(1:1)
   else
    tp = '?'
   end if
  end if

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1)

  alpha = ONE    ; beta = ZERO
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    allocate(w(ma)); w = ZERO ! 260104

   !if( na /= nv) call matrix_error('array2_matmul_lapack_mv: na /= nv (N)'   )
    if( na /= nv)then
     w = NaN()
     return
    end if
    if( ma == 0 .or. na == 0 ) return

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                

   case('T','C','t','c')

    allocate(w(na)); w = ZERO ! 260104

   !if( ma /= nv) call matrix_error('array2_matmul_lapack_mv: ma /= nv (T,C)' )
    if( ma /= nv)then
     w = NaN()
     return
    end if
    if( ma == 0 .or. na == 0 ) return

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)

   case('H','h')

    uplo = 'U'
    allocate(w(ma)); w = ZERO ! 260104

   !if( ma /= na) call matrix_error('array2_matmul_lapack_mv: ma /= na (H)'   )
   !if( na /= nv) call matrix_error('array2_matmul_lapack_mv: na /= nv (H)'   )
    if( ma /= na .or. na /= nv)then
     w = NaN()
     return
    end if
    if( ma == 0 ) return

    call zhemv(uplo ,ma,alpha,A,ma,v,incx,beta,w,incy)                                 

   case default

    allocate(w(ma)); w = ZERO ! 260104
    w = NaN()
    return
   !call               matrix_error('array2_matmul_lapack_mv: wrong type'     )

  end select


 end function       array2_matmul_lapack_mv

!-----------------------------------------------------------------------------------------------------------------------------------
!type='N'    w = A . v           w_i = A_ij  v_j
!type='T'    w = A^T . v = v . A w_i = v_j   A_ji
!type='S','H' w = A . v          w_i = A_ij  v_j    faster via dsymv, when A = A^T (upper triangle referenced, no check performed)
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_dv                                (A,v,type)               result(w)                       !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A
  real     (dp), dimension(:)  , contiguous, intent(in)       :: v
  character(*) , optional      , intent(in)                   :: type
  real     (dp), dimension(:)  , allocatable                  :: w
  integer                                                     :: ma,na,nv
  real     (dp)                                               :: alpha, beta
  character(1)                                                :: tp, uplo
  integer                                                     :: incx, incy

  tp = 'N'                                        ! default is w = A.v
  if(present(type))then
   if(len(type) >= 1)then
    tp = type(1:1)
   else
    tp = '?'
   end if
  end if

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1)

  alpha = 1.0_dp ; beta = 0.0_dp
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    allocate(w(ma)); w = 0.0_dp ! 260104

   !if( na /= nv) call matrix_error('array2_matmul_lapack_dv: na /= nv')
    if( na /= nv)then
     w = NaN(0.0_dp)
     return
    end if
    if( ma == 0 .or. na == 0 ) return

    call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                              

   case('T','C','t','c')

    allocate(w(na)); w = 0.0_dp ! 260104

   !if( ma /= nv) call matrix_error('array2_matmul_lapack_dv: ma /= nv')
    if( ma /= nv)then
     w = NaN(0.0_dp)
     return
    end if
    if( ma == 0 .or. na == 0 ) return

    call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                

   case('S','s','H','h')

    uplo = 'U'
    allocate(w(ma)); w = 0.0_dp ! 260104

    if( ma /= na .or. na /= nv )then
     w = NaN(0.0_dp)
     return
    end if
    if( ma == 0 ) return

    call dsymv(uplo,ma,alpha,A,ma,v,incx,beta,w,incy)                                       

   case default

   !call               matrix_error('array2_matmul_lapack_dv: wrong type')
    allocate(w(ma)); w = 0.0_dp ! 260104
    w = NaN(0.0_dp)
    return

  end select

 end function       array2_matmul_lapack_dv



!...................................................................................................................................
!.......................... Matrix Matrix Multiplication Subroutines ...............................................................!:.:.:
!...................................................................................................................................
!
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_mm_sub                            (A,B,C,mtype)                                            !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: A,B
  complex  (dp), dimension(:,:), contiguous, intent(inout)    :: C
  character(*) , optional                  , intent(in)       :: mtype
  integer                                                     :: ma,na,mb,nb,mc,nc
  complex  (dp)                                               :: alpha, beta
  character(1)                                                :: opA, opB, side,uplo
  character(mtype_len)                                        :: mt

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)

  ma = size(A,1) ; na = size(A,2)
  mb = size(B,1) ; nb = size(B,2)
  mc = size(C,1) ; nc = size(C,2)

  if( na /= mb .or. mc /= ma .or. nc /= nb) then
   C = NaN()
   return
  end if

  ! Lapack calculates C = alpha * A . B + beta * C
  !   C    =   A   x   B
  ! m x n    m x k   k x n
  ! call zgemm('N','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) where 'N'  is for normal, 'C' is for Hermitian conjugate
  ! If A is Hermitian, then we can call:
  ! call zhemm('L','U',M,N,  alpha,A,lda,B,ldb,beta,C,ldc) where side= 'L' is for A.B, side= 'R' for B.A and uplo='U'/'L' references the upper/lower part of the A matrix

  alpha = ONE ; beta = ZERO
  opA   = 'N' ; opB  = 'N'  ; side ='L' ; uplo = 'U'

  select case(mt(1:1))
  case('H','h')
    
   if( ma /= na ) then
    C = NaN()
    return
   end if
   if( mc == 0 .or. nc == 0 ) return

   call zhemm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Hermitian            

  case default

   if( mc == 0 .or. nc == 0 ) return
   if( na == 0 )then                    ! 261124: empty inner dimension: A.B = 0 (BUG-03)
    C = ZERO
    return
   end if

   call zgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                       

  end select

 end  subroutine    array2_matmul_lapack_mm_sub
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_dd_sub                            (A,B,C,mtype)                                            !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A,B
  real     (dp), dimension(:,:), contiguous, intent(inout)    :: C
  character(*) , optional                  , intent(in)       :: mtype
  integer                                                     :: ma,na,mb,nb,mc,nc
  real   (dp)                                                 :: alpha, beta
  character(1)                                                :: opA, opB, side,uplo
  character(mtype_len)                                        :: mt

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)

  ma = size(A,1) ; na = size(A,2)
  mb = size(B,1) ; nb = size(B,2)
  mc = size(C,1) ; nc = size(C,2)

  if( na /= mb .or. mc /= ma .or. nc /= nb) then
   C = NaN(0.0_dp)
   return
  end if

  ! Lapack calculates C = alpha * A . B + beta * C
  !   C    =   A   x   B
  ! m x n    m x k   k x n
  ! call dgemm('N','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) where 'N'  is for normal, 'T' is for transpose
  ! If A is Hermitian, then we can call:
  ! call dsymm('L','U',M,N,  alpha,A,lda,B,ldb,beta,C,ldc) where side= 'L' is for A.B, side= 'R' for B.A and uplo='U'/'L' references the upper/lower part of the A matrix

  alpha = 1.0_dp ; beta = 0.0_dp
  opA   = 'N' ; opB  = 'N'  ; side ='L' ; uplo = 'U'

  select case(mt(1:1))
  case('S','s')

   if( ma /= na )then
    C = NaN(0.0_dp)
    return
   end if
   if( mc == 0 .or. nc == 0 ) return

   call dsymm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Symmetric                  

  case default

   if( mc == 0 .or. nc == 0 ) return
   if( na == 0 )then                    ! 261124: empty inner dimension: A.B = 0 (BUG-03)
    C = 0.0_dp
    return
   end if
   
   call dgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                             

  end select

 end  subroutine    array2_matmul_lapack_dd_sub
!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^\dagger = array2_matmul_lapack_m_sub(A,C)    ;   A^\dagger . A = array2_matmul_lapack_m_sub(A,C,side='H')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_m_sub                             (A,C,side)                                               !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: A
  complex  (dp), dimension(:,:), contiguous, intent(inout)    :: C
  character(*) , optional                  , intent(in)       :: side
  integer                                                     :: ma,na, mc,nc
  real   (dp)                                                 :: alpha, beta
  character(1)                                                :: opA, uplo
  integer                                                     :: N,K,i,j

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  ma  = size(A,1) ; na = size(A,2)
  mc  = size(C,1) ; nc = size(C,2)

  opA = 'N';                       ! default is A . A^\dagger
  N   =  ma       ; K  = na

  if(present(side))then
   if(len(side) >= 1)then
    select case(side(1:1))
    case('N','n',' ')                       ! default: A . A^\dagger (A . A^T)
    case('C','c','T','t','H','h','S','s','L','l')

     opA = 'C'
     N   = na      ; K  = ma

    case default                            ! 261124: unrecognised side -> NaN (BUG-05)
     opA = '?'
    end select
   end if
  end if

  ! Lapack calculates C = alpha * A . A^\dagger + beta * C
  ! if opA='N' then   C    =   A   x   A^\dagger          ::          if opA='C'    C    =   A^\dagger    x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m          m x n
  !
  ! call zherk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='C'

  alpha = 1.0_dp ; beta = 0.0_dp; uplo = 'U'

  if( opA == '?' )then                   ! 261124: unrecognised side (BUG-05)
   C = NaN()
   return
  end if
  if( mc /= N .or. nc /= N)then
   C = NaN()
   return
  end if
  if( N == 0 ) return
  if( K == 0 )then                     ! 261124: empty inner dimension: A.A^+ = 0 (BUG-03)
   C = ZERO
   return
  end if

  call zherk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                    

  !only the upper triangle part of C is filled:
  do concurrent(i=1:n, j=1:n, i>j)
   C(i,j) = CONJG(C(j,i))
  end do

 end  subroutine    array2_matmul_lapack_m_sub
!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^T = array2_matmul_lapack_d_sub(A,C)    ;   A^T . A = array2_matmul_lapack_d_sub(A,C,side='S')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_d_sub                             (A,C,side)                                               !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A
  real     (dp), dimension(:,:), contiguous, intent(inout)    :: C
  character(*) , optional                  , intent(in)       :: side
  integer                                                     :: ma,na, mc,nc
  real     (dp)                                               :: alpha, beta
  character(1)                                                :: opA, uplo
  integer                                                     :: N,K,i,j

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  ma  = size(A,1) ; na = size(A,2)
  mc  = size(C,1) ; nc = size(C,2)

  opA = 'N';                       ! default is A . A^T
  N   =  ma       ; K  = na

  if(present(side))then
   if(len(side) >= 1)then
    select case(side(1:1))
    case('N','n',' ')                       ! default: A . A^\dagger (A . A^T)
    case('C','c','T','t','H','h','S','s','L','l')

     opA = 'T'
     N   = na      ; K  = ma

    case default                            ! 261124: unrecognised side -> NaN (BUG-05)
     opA = '?'
    end select
   end if
  end if

  ! Lapack calculates C = alpha * A . A^T + beta * C
  ! if opA='N' then   C    =   A   x   A^T                ::          if opA='T'    C    =   A^T          x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m          m x n
  !
  ! call dsyrk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'

  alpha = 1.0_dp ; beta = 0.0_dp; uplo = 'U'

  if( opA == '?' )then                   ! 261124: unrecognised side (BUG-05)
   C = NaN(0.0_dp)
   return
  end if
  if( mc /= N .or. nc /= N)then
   C = NaN(0.0_dp)
   return
  end if
  if( N == 0 ) return
  if( K == 0 )then                     ! 261124: empty inner dimension: A.A^+ = 0 (BUG-03)
   C = 0.0_dp
   return
  end if

  call dsyrk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                    

  !only the upper triangle part of C is filled:
  do concurrent(i=1:n, j=1:n, i>j)
   C(i,j) = C(j,i)
  end do

 end  subroutine    array2_matmul_lapack_d_sub


!...................................................................................................................................
!.......................... Matrix Vector Multiplication Subroutines ...............................................................!:.:.:
!...................................................................................................................................
!
!-----------------------------------------------------------------------------------------------------------------------------------
!type='N'    w = A . v           w_i = A_ij  v_j    (default)
!type='T'    w = A^T . v = v . A w_i = v_j   A_ji
!type='H'    w = A . v           w_i = A_ij  v_j    faster, when A = A^\dagger (no check performed)
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_mv_sub                            (A,v,w,type)                                             !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: A
  complex  (dp), dimension(:)  , contiguous, intent(in)       :: v
  complex  (dp), dimension(:)  , contiguous, intent(inout)    :: w
  character(*) , optional                  , intent(in)       :: type
  integer                                                     :: ma,na,nv,nw
  complex  (dp)                                               :: alpha, beta
  character(1)                                                :: tp,uplo
  integer                                                     :: incx, incy

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)
  tp = 'N'                                        ! default is w = A.v
  if(present(type))then
   if(len(type) >= 1)then
    tp = type(1:1)
   else
    tp = '?'
   end if
  end if

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = ONE    ; beta = ZERO
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    if( na /= nv .or. nw /= ma )then
     w = NaN()
     return
    end if
    if( ma == 0 ) return
    if( na == 0 )then                    ! 261124: empty inner dimension: result = 0 (BUG-03)
     w = ZERO
     return
    end if

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                        

   case('T','C','t','c')

    if( ma /= nv .or. nw /= na )then
     w = NaN()
     return
    end if
    if( na == 0 ) return
    if( ma == 0 )then                    ! 261124: empty inner dimension: result = 0 (BUG-03)
     w = ZERO
     return
    end if

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)

   case('H','h')

    uplo = 'U'

    if( ma /= na .or. na /= nv .or. nw /= ma )then
     w = NaN()
     return
    end if
    if( ma == 0 ) return

    call zhemv(uplo ,ma,alpha,A,ma,v,incx,beta,w,incy)                                         

   case default

    w = NaN()
    return

  end select


 end  subroutine    array2_matmul_lapack_mv_sub
!-----------------------------------------------------------------------------------------------------------------------------------
!type='T'    w = A^T . v = v . A w_i = v_j   A_ji  
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_vm_sub                            (v,A,w)                                                  !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: A
  complex  (dp), dimension(:)  , contiguous, intent(in)       :: v
  complex  (dp), dimension(:)  , contiguous, intent(inout)    :: w
  integer                                                     :: ma,na,nv,nw
  complex  (dp)                                               :: alpha, beta
  character(1)                                                :: tp     
  integer                                                     :: incx, incy

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  tp = 'T'                                        

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = ONE    ; beta = ZERO
  incx  = 1      ; incy = 1

  if( ma /= nv .or. nw /= na )then
   w = NaN()
   return
  end if
  if( na == 0 ) return
  if( ma == 0 )then                    ! 261124: empty inner dimension: result = 0 (BUG-03)
   w = ZERO
   return
  end if

  call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                             


 end  subroutine    array2_matmul_lapack_vm_sub
!-----------------------------------------------------------------------------------------------------------------------------------
!type='N'    w = A . v           w_i = A_ij  v_j
!type='T'    w = A^T . v = v . A w_i = v_j   A_ji
!type='S','H' w = A . v          w_i = A_ij  v_j    faster via dsymv, when A = A^T (upper triangle referenced, no check performed)
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_dv_sub                            (A,v,w,type)                                             !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A
  real     (dp), dimension(:)  , contiguous, intent(in)       :: v
  real     (dp), dimension(:)  , contiguous, intent(inout)    :: w
  character(*) , optional                  , intent(in)       :: type
  integer                                                     :: ma,na,nv,nw
  real     (dp)                                               :: alpha, beta
  character(1)                                                :: tp,uplo
  integer                                                     :: incx, incy

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  tp = 'N'                                        ! default is w = A.v
  if(present(type))then
   if(len(type) >= 1)then
    tp = type(1:1)
   else
    tp = '?'
   end if
  end if

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = 1.0_dp ; beta = 0.0_dp
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    if( na /= nv .or. nw /= ma )then
     w = NaN(0.0_dp)
     return
    end if
    if( ma == 0 ) return
    if( na == 0 )then                    ! 261124: empty inner dimension: result = 0 (BUG-03)
     w = 0.0_dp
     return
    end if

    call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                           

   case('T','C','t','c')

    if( ma /= nv .or. nw /= na )then
     w = NaN(0.0_dp)
     return
    end if
    if( na == 0 ) return
    if( ma == 0 )then                    ! 261124: empty inner dimension: result = 0 (BUG-03)
     w = 0.0_dp
     return
    end if

    call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                               

   case('S','s','H','h')

    uplo = 'U'

    if( ma /= na .or. na /= nv .or. nw /= ma )then
     w = NaN(0.0_dp)
     return
    end if
    if( ma == 0 ) return

    call dsymv(uplo,ma,alpha,A,ma,v,incx,beta,w,incy)                                                             

   case default

    w = NaN(0.0_dp)
    return

  end select

 end  subroutine    array2_matmul_lapack_dv_sub
!-----------------------------------------------------------------------------------------------------------------------------------
!type='T'    w = A^T . v = v . A w_i = v_j   A_ji  
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_vd_sub                            (v,A,w)                                                  !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: A
  real     (dp), dimension(:)  , contiguous, intent(in)       :: v
  real     (dp), dimension(:)  , contiguous, intent(inout)    :: w
  integer                                                     :: ma,na,nv,nw
  real     (dp)                                               :: alpha, beta
  character(1)                                                :: tp     
  integer                                                     :: incx, incy

  ! 261124: no unconditional pre-zeroing: BLAS with beta=0 does not read the output (BUG-03)

  tp = 'T'                                        

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = 1.0_dp ; beta = 0.0_dp
  incx  = 1      ; incy = 1

  if( ma /= nv .or. nw /= na )then
   w = NaN(0.0_dp)
   return
  end if
  if( na == 0 ) return
  if( ma == 0 )then                    ! 261124: empty inner dimension: result = 0 (BUG-03)
   w = 0.0_dp
   return
  end if

  call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                       


 end  subroutine    array2_matmul_lapack_vd_sub



!...................................................................................................................................
!.......................... Matrix Inversion .......................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_inverse                                         (C)                      result(CI)                      !:.:.:
  complex  (dp), dimension(:,:), contiguous, intent(in)       :: C
! complex  (dp), dimension(size(C,1),size(C,2))               :: CI
  complex  (dp), dimension(:,:), allocatable                  :: CI   ! 260104
  integer                                                     :: n, LWORK, info
! integer      , dimension(size(C,1))                         :: ipiv
  integer      , dimension(:)  , allocatable                  :: ipiv ! 260104
  complex(dp)  , dimension(:)  , allocatable                  :: WORK
  complex(dp)  , dimension(1)                                 :: WORKTEST
  character(1024)                                             :: err

  if( size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   allocate(CI(size(C,1),size(C,2))); CI = NaN()
   return
  end if

  n = size(C,1)
  allocate(CI(n,n))
  if( n == 0 ) return
  allocate(ipiv(n))

  !---------------------------------------------------------
  !Use zgetrf to compute LU decomposition: MATIV = P * L * U
  CI = C
  call zgetrf(n,n,CI,n,ipiv,info)                                                                    
  if(info <  0) write(err,*)'array2_inverse:  zgetrf failed to perform LU decomposition. INFO= ',info
  if(info <  0) call matrix_error(err)
  if(info >  0)then ! 261124: singular matrix -> NaN (was: matrix_error)
   CI = NaN()
   return
  end if
  !---------------------------------------------------------
  !Use zgetri to compute inverse:
  LWORK = -1 ! compute optimal LWORK
  call zgetri(n,CI,n,ipiv,WORKTEST,LWORK,info)
  if(info /= 0) write(err,*)'array2_inverse: zgetri failed to determine LWORK. INFO= ',info
  if(info /= 0) call matrix_error(err)
  LWORK = max(1, int(real(WORKTEST(1), kind=dp)) + 10)
  !---------------------------------------------------------
  !Compute inverse:
  ALLOCATE(WORK(LWORK))
  call zgetri(n,CI,n,ipiv,WORK    ,LWORK,info)                                                     
  if(info <  0) write(err,*)'array2_inverse: zgetri failed to compute inverse. INFO= ',info
  if(info <  0) call matrix_error(err)
  if(info >  0)then ! 261124: singular matrix -> NaN (was: matrix_error)
   CI = NaN()
   return
  end if
  !---------------------------------------------------------
  DEALLOCATE(WORK)
  deallocate(ipiv)
 end function array2_inverse
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_inverse_d                                       (C)                      result(CI)                      !:.:.:
  real     (dp), dimension(:,:), contiguous, intent(in)       :: C
! real     (dp), dimension(size(C,1),size(C,2))               :: CI
  real     (dp), dimension(:,:), allocatable                  :: CI ! 260104
  integer                                                     :: n, LWORK, info
! integer      , dimension(size(C,1))                         :: ipiv
  integer      , dimension(:)  , allocatable                  :: ipiv ! 260104
  real   (dp)  , dimension(:)  , allocatable                  :: WORK
  real   (dp)  , dimension(1)                                 :: WORKTEST
  character(1024)                                             :: err

  if( size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   allocate(CI(size(C,1),size(C,2))); CI = NaN(0.0_dp)
   return
  end if

  n = size(C,1)
  allocate(CI(n,n))
  if( n == 0 ) return
  allocate(ipiv(n))

  !---------------------------------------------------------
  !Use dgetrf to compute LU decomposition: MATIV = P * L * U
  CI = C
  call dgetrf(n,n,CI,n,ipiv,info)                                                                 
  if(info <  0) write(err,*)'array2_inverse_d:  dgetrf failed to perform LU decomposition. INFO= ',info
  if(info <  0) call matrix_error(err)
  if(info >  0)then ! 261124: singular matrix -> NaN (was: matrix_error)
   CI = NaN(0.0_dp)
   return
  end if
  !---------------------------------------------------------
  !Use dgetri to compute inverse:
  LWORK = -1 ! compute optimal LWORK
  call dgetri(n,CI,n,ipiv,WORKTEST,LWORK,info)
  if(info /= 0) write(err,*)'array2_inverse_d: dgetri failed to determine LWORK. INFO= ',info
  if(info /= 0) call matrix_error(err)
  LWORK = max(1, int(WORKTEST(1)) + 10)
  !---------------------------------------------------------
  !Compute inverse:
  ALLOCATE(WORK(LWORK))
  call dgetri(n,CI,n,ipiv,WORK    ,LWORK,info)                                                     
  if(info <  0) write(err,*)'array2_inverse_d: dgetri failed to compute inverse. INFO= ',info
  if(info <  0) call matrix_error(err)
  if(info >  0)then ! 261124: singular matrix -> NaN (was: matrix_error)
   CI = NaN(0.0_dp)
   return
  end if
  !---------------------------------------------------------
  DEALLOCATE(WORK)
  deallocate(ipiv)
 end function array2_inverse_d
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... Eigenvalues - Eigenvectors .............................................................................!:.:.:
!...................................................................................................................................
! Notice that for large enough matrices, the order of the eigenvalues for zgeev and dgeev don't come with the same order when job='N'/'V'
! In order to compare, you have to sort them in e.g. real part order and then the difference is pure numerical error. 
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvalues                                     (C,mtype,sortby)         result(eigenval)                !:.:.:
  complex  (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  character(*)                              , optional, intent(in)                     :: mtype, sortby
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec
  character(mtype_len)                                                     :: mt,job

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)
  job = 'N' !only eigenvalues

  select case(mt(1:1))
  case('H','h') !hermitian matrix
   call array2_zheev(C,eigenval,eigenvec,job)
  case default  !general   matrix
   call array2_zgeev(C,eigenval,eigenvec,job)
  end select

! sort eigenvalues, if asked:
   if(present(sortby))  eigenval       = sort(eigenval,by=sortby)

 end function       array2_eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvectors                                    (C,mtype,sortby)         result(evs)                     !:.:.:
  complex  (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  character(*)                              , optional, intent(in)                     :: mtype, sortby
  integer      , dimension(:)               , allocatable                  :: pos
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec, eigenvectmp
  complex  (dp), dimension(:,:)             , allocatable                  :: evs
  character(mtype_len)                                                     :: mt,job
  integer                                                                  :: n,i

  n   = size(C,1)

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)
  job = 'V'     ! eigenvectors too

  select case(mt(1:1))
  case('H','h') ! hermitian matrix
   call array2_zheev(C,eigenval,eigenvec,job)
  case default  ! general   matrix
   call array2_zgeev(C,eigenval,eigenvec,job)
  end select

! sort eigenvalues and eigenvectors, if asked:
  if(present(sortby))then
   allocate   (pos(n),eigenvectmp(n,n))
   eigenval       = sort(eigenval,P=pos,by=sortby)
   eigenvectmp    = eigenvec
   do i = 1, n
    eigenvec(:,i) = eigenvectmp(:,pos(i))
   end do
   deallocate (pos   ,eigenvectmp)
  end if

! store result:
  allocate(evs(n,n+1))

  if( n > 0 )then
   evs(:,1 ) = eigenval         ! The first column has the eigenvalues
   evs(:,2:) = eigenvec         ! Columns 2-(n+1)  has the eigenvectors
  end if

 end function       array2_eigenvectors
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvalues_d                                   (C,mtype,sortby)                result(eigenval)         !:.:.:
  real     (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  character(*)                              , optional, intent(in)                     :: mtype, sortby
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec
  character(mtype_len)                                                     :: mt,job

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)
  job = 'N'     ! only eigenvalues

  select case(mt(1:1))
  case('S','s') ! symmetric matrix
   call array2_dsyev(C,eigenval,eigenvec,job)
  case default  ! general   matrix
   call array2_dgeev(C,eigenval,eigenvec,job)
  end select

! sort eigenvalues, if asked:
   if(present(sortby)) eigenval       = sort(eigenval,by=sortby)

 end function       array2_eigenvalues_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvectors_d                                  (C,mtype,sortby)         result(evs)                     !:.:.:
  real     (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  character(*)                              , optional, intent(in)                     :: mtype, sortby
  integer      , dimension(:)               , allocatable                  :: pos
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec, eigenvectmp
  complex  (dp), dimension(:,:)             , allocatable                  :: evs
  character(mtype_len)                                                     :: mt,job
  integer                                                                  :: n,i

  n   = size(C,1)

  mt = mtype_char(mtype)                ! 261124: accepts 'ZH','DS',... and one-character 'H','S' (BUG-04)
  job = 'V'     ! eigenvectors too

  select case(mt(1:1))
  case('S','s') ! symmetric matrix
   call array2_dsyev(C,eigenval,eigenvec,job)
  case default  ! general   matrix
   call array2_dgeev(C,eigenval,eigenvec,job)
  end select

! sort eigenvalues and eigenvectors, if asked:
  if(present(sortby))then
   allocate   (pos(n),eigenvectmp(n,n))
   eigenval       = sort(eigenval,P=pos,by=sortby)
   eigenvectmp    = eigenvec
   do i = 1, n
    eigenvec(:,i) = eigenvectmp(:,pos(i))
   end do
   deallocate (pos   ,eigenvectmp)
  end if

! store result:
  allocate(evs(n,n+1))

  if( n > 0 )then
   evs(:,1 ) = eigenval        ! The first column has the eigenvalues
   evs(:,2:) = eigenvec        ! Columns 2-(n+1)  has the eigenvectors. eigenvec is real, evs is complex
  end if

 end function       array2_eigenvectors_d
!-----------------------------------------------------------------------------------------------------------------------------------
! call array2_zgeev(C,E,V,job='N'/'V')  'N': eigenvalues only  'V': eigenvalues and right eigenvectors: C . v = λ v
!
! For eigenvalue i:  λ_i = E(i)    v_i = V(:,i)   V = (v_1 ... v_n)    Λ = diag(λ1,...,λn)   then       C . V = V . Λ
!                    MATMUL( C, V(:,i) ) - E(i) * V(:,i) = 0      or   MATMUL(C,V) - MATMUL(V,diagonal(E)) = 0
!
! For eigenvalues:  eigenvectors(1,n)   n = size(C,1)    only allocatable arrays
! For eigenvectors: eigenvectors(n,n)
!
! The upper:lower bounds for eiganvalues/eigenvectors is changed to 1:n, otherwise the result is not correct.
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_zgeev                                           (C,eigenval,eigenvec,job)                                !:.:.:
  complex  (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
! complex  (dp), dimension(size(C,1),size(C,2))                            :: A
  complex  (dp), dimension(:,:)                             , allocatable  :: A ! 260104
  integer                                                                  :: n
  character(1)                                                             :: JOBVL, JOBVR, jch
  integer                                                                  :: info,  LWORK, LDVL,LDVR
  complex  (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: RWORK
  complex  (dp), dimension(:,:)                             , allocatable  :: VL,VR
  complex  (dp), dimension(:)                               , allocatable  :: EV,WORK
  character(1024)                                                          :: err


  n = size(C,1)
  jch = ' ' ; if(len(job) >= 1) jch = job(1:1)

  select case(jch)
  case('N','n')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if( n == 0 ) return
  case('V','v')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if(allocated(eigenvec)) deallocate(eigenvec) ; allocate(eigenvec(n,n))
   if( n == 0 ) return
  case default
   call matrix_error('array2_zgeev: wrong value of job: '//trim(job))
  end select
  if( size(C,1) /= size(C,2) )then ! 261124: non-square input -> NaN (was: matrix_error)
   eigenval = NaN()
   if(jch == 'V' .or. jch == 'v') eigenvec = NaN()
   return
  end if

  allocate(A(n,n)) ; allocate(RWORK(2*n)) ; allocate(EV(n))

  A = C  ! A local copy of the matrix C

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(jch)
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBVL = 'N' ; JOBVR = 'N'; LDVL = 1; LDVR = 1; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,RWORK,info)
   if(info /= 0) write(err,*)                      'array2_zgeev: zgeev LWORK query failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(real(WORKTEST(1), kind=dp)) + 10)
   allocate(WORK(LWORK))
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORK    ,LWORK,RWORK,info)                            
   if(info <  0) write(err,*)                      'array2_zgeev: zgeev failed (N). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    return
   end if
   eigenval = EV
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   JOBVL = 'N' ; JOBVR = 'V'; LDVL = 1; LDVR = n; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,RWORK,info)
   if(info /= 0) write(err,*)                      'array2_zgeev: zgeev LWORK query failed (V). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(real(WORKTEST(1), kind=dp)) + 10)
   allocate(WORK(LWORK))
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORK    ,LWORK,RWORK,info)
   if(info <  0) write(err,*)                      'array2_zgeev: zgeev failed (V). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    eigenvec = NaN()
    return
   end if
   eigenval  = EV
   eigenvec = VR
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(RWORK) ; deallocate(WORK)
  deallocate(VL)    ; deallocate(VR)   ; deallocate(EV)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_zgeev
!-----------------------------------------------------------------------------------------------------------------------------------
!Works as zgeev, but C(:,:) is assumed to be Hermitian (no test performed). 
!
!The ordering of the (real) eigenvalues is by ascending order (not true from zgeev, so don't make any attempt to compare without sorting)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_zheev                                           (C,eigenval,eigenvec,job)                                !:.:.:
  complex  (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
! complex  (dp), dimension(size(C,1),size(C,2))                            :: A
  complex  (dp), dimension(:,:)                             , allocatable  :: A ! 260104
  integer                                                                  :: n
  integer                                                                  :: info,  LWORK
  complex  (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: RWORK
  real     (dp), dimension(:)                               , allocatable  :: EV
  complex  (dp), dimension(:)                               , allocatable  :: WORK
  character(1024)                                                          :: err
  character(1)                                                             :: JOBZ, UPLO, jch


  n = size(C,1)
  jch = ' ' ; if(len(job) >= 1) jch = job(1:1)

  select case(jch)
  case('N','n')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if( n == 0 ) return
  case('V','v')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if(allocated(eigenvec)) deallocate(eigenvec) ; allocate(eigenvec(n,n))
   if( n == 0 ) return
  case default
   call matrix_error('array2_zheev: wrong value of job: '//trim(job))
  end select
  if( size(C,1) /= size(C,2) )then ! 261124: non-square input -> NaN (was: matrix_error)
   eigenval = NaN()
   if(jch == 'V' .or. jch == 'v') eigenvec = NaN()
   return
  end if

  allocate(A(n,n)) ; allocate(RWORK(max(1, 3*n-2))) ; allocate(EV(n))

  A = C  ! A local copy of the matrix C

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(jch)
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBZ='N' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call zheev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,RWORK,info)                                                   
   if(info /= 0) write(err,*)                      'array2_zheev: zheev LWORK query failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(real(WORKTEST(1), kind=dp)) + 10)
   allocate(WORK(LWORK))
   call zheev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,RWORK,info)                                               
   if(info <  0) write(err,*)                      'array2_zheev: zheev failed (N). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    return
   end if
   eigenval = EV
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   JOBZ='V' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call zheev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,RWORK,info)
   if(info /= 0) write(err,*)                      'array2_zheev: zheev LWORK query failed (V). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(real(WORKTEST(1), kind=dp)) + 10)
   allocate(WORK(LWORK))
   call zheev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,RWORK,info)
   if(info <  0) write(err,*)                      'array2_zheev: zheev failed (V). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    eigenvec = NaN()
    return
   end if
   eigenval = EV
   eigenvec = A
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(RWORK) ; deallocate(WORK)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_zheev
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_dgeev                                           (C,eigenval,eigenvec,job)                                !:.:.:
  real     (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
! real     (dp), dimension(size(C,1),size(C,2))                            :: A 
  real     (dp), dimension(:,:)                             , allocatable  :: A ! 260104
  integer                                                                  :: n,j
  character(1)                                                             :: JOBVL, JOBVR, jch
  integer                                                                  :: info,  LWORK, LDVL,LDVR
  real     (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:,:)                             , allocatable  :: VL,VR
  real     (dp), dimension(:)                               , allocatable  :: ReEV,ImEV,WORK
  character(1024)                                                          :: err


  n = size(C,1)
  jch = ' ' ; if(len(job) >= 1) jch = job(1:1)

  select case(jch)
  case('N','n')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if( n == 0 ) return
  case('V','v')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if(allocated(eigenvec)) deallocate(eigenvec) ; allocate(eigenvec(n,n))
   if( n == 0 ) return
  case default
   call matrix_error('array2_dgeev: wrong value of job: '//trim(job))
  end select
  if( size(C,1) /= size(C,2) )then ! 261124: non-square input -> NaN (was: matrix_error)
   eigenval = NaN()
   if(jch == 'V' .or. jch == 'v') eigenvec = NaN()
   return
  end if

  allocate(A(n,n)) ; allocate(ReEV(n)) ; allocate(ImEV(n))

  A = C  ! A local copy of the matrix C

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(jch)
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBVL = 'N' ; JOBVR = 'N'; LDVL = 1; LDVR = 1; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,info)
   if(info /= 0) write(err,*)                      'array2_dgeev: dgeev LWORK query failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(WORKTEST(1)) + 10)
   allocate(WORK(LWORK))
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORK    ,LWORK,info)                            
   if(info <  0) write(err,*)                      'array2_dgeev: dgeev failed (N). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    return
   end if
   eigenval = cmplx(ReEV,ImEV,kind=dp)
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   JOBVL = 'N' ; JOBVR = 'V'; LDVL = 1; LDVR = n; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,info)
   if(info /= 0) write(err,*)                      'array2_dgeev: dgeev LWORK query failed (V). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(WORKTEST(1)) + 10)
   allocate(WORK(LWORK))
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORK    ,LWORK,info)
   if(info <  0) write(err,*)                      'array2_dgeev: dgeev failed (V). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    eigenvec = NaN()
    return
   end if
   eigenval = cmplx(ReEV,ImEV,kind=dp)
   ! If eigenvalues form a complex conjugate pair, dgeev stores ImEV(j) > 0, ImEV(j+1) = -ImEV(j) < 0,
   ! v(j) = VR(:,j) + i*VR(:,j+1), and v(j+1) = VR(:,j) - i*VR(:,j+1).
   j = 1
   do while (j <= n)
    if(ImEV(j) > 0.0_dp .and. j < n)then
     eigenvec(:,j  ) = cmplx(VR(:,j), VR(:,j+1),kind=dp)
     eigenvec(:,j+1) = cmplx(VR(:,j),-VR(:,j+1),kind=dp)
     j = j + 2
    else
     eigenvec(:,j  ) = cmplx(VR(:,j), 0.0_dp,   kind=dp)
     j = j + 1
    end if
   end do
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(WORK)
  deallocate(VL)    ; deallocate(VR)   ; deallocate(ReEV)  ; deallocate(ImEV)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_dgeev
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_dsyev                                           (C,eigenval,eigenvec,job)                                !:.:.:
  real     (dp), dimension(:,:)             , contiguous, intent(in)       :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
! real     (dp), dimension(size(C,1),size(C,2))                            :: A
  real     (dp), dimension(:,:)                             , allocatable  :: A ! 260104
  integer                                                                  :: n
  integer                                                                  :: info,  LWORK
  real     (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: EV
  real     (dp), dimension(:)                               , allocatable  :: WORK
  character(1024)                                                          :: err
  character(1)                                                             :: JOBZ, UPLO, jch


  n = size(C,1)
  jch = ' ' ; if(len(job) >= 1) jch = job(1:1)

  select case(jch)
  case('N','n')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if( n == 0 ) return
  case('V','v')
   if(allocated(eigenval)) deallocate(eigenval) ; allocate(eigenval(n))
   if(allocated(eigenvec)) deallocate(eigenvec) ; allocate(eigenvec(n,n))
   if( n == 0 ) return
  case default
   call matrix_error('array2_dsyev: wrong value of job: '//trim(job))
  end select
  if( size(C,1) /= size(C,2) )then ! 261124: non-square input -> NaN (was: matrix_error)
   eigenval = NaN()
   if(jch == 'V' .or. jch == 'v') eigenvec = NaN()
   return
  end if

  allocate(A(n,n)) ; allocate(EV(n))

  A = C  ! A local copy of the matrix C

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(jch)
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBZ='N' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,info)
   if(info /= 0) write(err,*)                      'array2_dsyev: dsyev LWORK query failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(WORKTEST(1)) + 10)
   allocate(WORK(LWORK))
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,info)                                                      
   if(info <  0) write(err,*)                      'array2_dsyev: dsyev failed (N). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    return
   end if
   eigenval = EV
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   JOBZ='V' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,info)
   if(info /= 0) write(err,*)                      'array2_dsyev: dsyev LWORK query failed (V). INFO= ',info
   if(info /= 0) call matrix_error(err)
   LWORK = max(1, int(WORKTEST(1)) + 10)
   allocate(WORK(LWORK))
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,info)
   if(info <  0) write(err,*)                      'array2_dsyev: dsyev failed (V). INFO= ',info
   if(info <  0) call matrix_error(err)
   if(info >  0)then ! 261124: LAPACK did not converge -> NaN (was: matrix_error)
    eigenval = NaN()
    eigenvec = NaN()
    return
   end if
   eigenval = EV
   eigenvec = A
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(WORK)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_dsyev
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_determinant                                     (C)  result (det)                                        !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp)                                            :: det
! complex(dp), dimension(size(C,1),size(C,2))            :: A
  complex(dp), dimension(:,:), allocatable               :: A     ! 260104
  integer                                                :: n
  integer                                                :: info, isign, i
! integer    , dimension(size(C,1))                      :: ipiv
  integer    , dimension(:)  , allocatable               :: ipiv  ! 260104
  character(1000)                                        :: err

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   det = NaN()
   return
  end if

  n = size(C,1)
  if (n == 0) then
   det = ONE
   return
  end if

  allocate(A   (n,n))
  allocate(ipiv(n))

  A = C

  call zgetrf(n,n,A,n,ipiv,info)

  if (info < 0) then
   write(err,*) 'array2_determinant: zgetrf failed. INFO= ', info
   call matrix_error(trim(err))
  end if
  if (info > 0) then
   det = ZERO
   return
  end if

  isign = 1; det   = ONE

! do concurrent(i  = 1:n)   ! 251202
!   ! $omp parallel do reduction(*:det)
  do i = 1, n
   det  = det * A(i,i)
  end do
!   ! $omp end parallel do

! do concurrent(i  = 1:n)   ! 251202
  do i = 1, n
   if(ipiv(i) /= i) isign = -isign
  end do
  
  det = det * isign

 end function       array2_determinant
!-----------------------------------------------------------------------------------------------------------------------------------
! For n x n matrices with gaussian random numbers, n >~ 192 requires the use of lnDet
!
! det(1) = log(|det(C)|)  det(2) = exp( i arg( det(C) ) )
! Det(C) = exp(det(1))  * det(2)
! For singular C (det(C) = 0): det(1) = (-Inf, 0), det(2) = (0, 0), so exp(det(1)) * det(2) = (0, 0)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_determinant                                 (C)                      result(det)                     !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp), dimension(2)                              :: det
! complex(dp), dimension(size(C,1),size(C,2))            :: A
  complex(dp), dimension(:,:), allocatable               :: A     ! 260104
  integer                                                :: n
  integer                                                :: info, i
! integer    , dimension(size(C,1))                      :: ipiv
  integer    , dimension(:)  , allocatable               :: ipiv  ! 260104
  complex(dp)                                            :: phase
  real   (dp)                                            :: aAbs,lnDetAbs
  character(1000)                                        :: err

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   det = NaN()
   return
  end if

  n = size(C,1)
  if (n == 0) then
   det(1) = ZERO
   det(2) = ONE
   return
  end if

  allocate(A   (n,n))
  allocate(ipiv(n))

  A = C

  call zgetrf(n,n,A,n,ipiv,info)

  if (info < 0) then
   write(err,*) 'array2_log_determinant: zgetrf failed. INFO= ', info
   call matrix_error(trim(err))
  end if
  if (info > 0) then
   det(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
   det(2) = ZERO
   return
  end if

  lnDetAbs  = 0.0_dp
  phase     = ONE

  do i      = 1, n

   aAbs     = ABS(A(i,i))

   if(aAbs <= 0.0_dp) then
    det(1)  = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    det(2)  = ZERO
    return
   end if

   lnDetAbs = lnDetAbs + LOG(aAbs)
   phase    = phase    * (A(i,i)/aAbs)

   if(ipiv(i) /= i )  phase = -phase

  end do

  det(1)    = cmplx(lnDetAbs, 0.0_dp, kind=dp)
  det(2)    = phase

 end function       array2_log_determinant
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_determinant_d                                   (C)                      result(det)                     !:.:.:
  real   (dp), dimension(:,:), contiguous, intent(in)    :: C
  real   (dp)                                            :: det
! real   (dp), dimension(size(C,1),size(C,2))            :: A
  real   (dp), dimension(:,:), allocatable               :: A    ! 260104
  integer                                                :: n
  integer                                                :: info, isign, i
! integer    , dimension(size(C,1))                      :: ipiv
  integer    , dimension(:)  , allocatable               :: ipiv ! 260104
  character(1000)                                        :: err

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   det = NaN(0.0_dp)
   return
  end if

  n = size(C,1)
  if (n == 0) then
   det = 1.0_dp
   return
  end if

  allocate(A   (n,n))
  allocate(ipiv(n))

  A = C

  call dgetrf(n,n,A,n,ipiv,info)

  if (info < 0) then
   write(err,*) 'array2_determinant_d: dgetrf failed. INFO= ', info
   call matrix_error(trim(err))
  end if
  if (info > 0) then
   det = 0.0_dp
   return
  end if

  isign = 1; det   = 1.0_dp

! do concurrent(i  = 1:n)  ! 251202
!  ! $omp parallel do reduction(*:det)
  do i = 1, n
   det  = det * A(i,i)
  end do
!  ! $omp end parallel do

! do concurrent(i  = 1:n)  ! 251202
  do i = 1, n
   if(ipiv(i) /= i) isign = -isign
  end do
  
  det = det * isign

 end function       array2_determinant_d
!-----------------------------------------------------------------------------------------------------------------------------------
! For n x n matrices with gaussian random numbers, n >~ 192 requires the use of lnDet
!
! det(1) = log(|det(C)|)  det(2) = sign(det(C)) (+1 or -1)
! Det(C) = exp(det(1))  * det(2)
! For singular C (det(C) = 0): det(1) = -Inf, det(2) = 0.0_dp, so exp(det(1)) * det(2) = 0.0_dp
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_determinant_d                               (C)                      result(det)                     !:.:.:
  real   (dp), dimension(:,:), contiguous, intent(in)    :: C
  real   (dp), dimension(2)                              :: det
! real   (dp), dimension(size(C,1),size(C,2))            :: A
  real   (dp), dimension(:,:), allocatable               :: A    ! 260104
  integer                                                :: n
  integer                                                :: info, i, isign
! integer    , dimension(size(C,1))                      :: ipiv
  integer    , dimension(:)  , allocatable               :: ipiv ! 260104
  real   (dp)                                            :: aAbs,lnDetAbs
  character(1000)                                        :: err

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   det = NaN(0.0_dp)
   return
  end if

  n = size(C,1)
  if (n == 0) then
   det(1) = 0.0_dp
   det(2) = 1.0_dp
   return
  end if

  allocate(A(n,n))
  allocate(ipiv(n))

  A = C

  call dgetrf(n,n,A,n,ipiv,info)

  if (info < 0) then
   write(err,*) 'array2_log_determinant_d: dgetrf failed. INFO= ', info
   call matrix_error(trim(err))
  end if
  if (info > 0) then
   det(1) = IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF)
   det(2) = 0.0_dp
   return
  end if

  lnDetAbs  = 0.0_dp
  isign     = 1

  do i      = 1, n

   aAbs     = A(i,i)

   if(aAbs == 0.0_dp) then
    det(1)  = IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF)
    det(2)  = 0.0_dp
    return
   end if

   if(aAbs < 0.0_dp) then
    aAbs    = -aAbs
    isign   = -isign
   end if

   lnDetAbs = lnDetAbs + LOG(aAbs)

   if(ipiv(i) /= i )  isign = -isign

  end do

  det(1)    = lnDetAbs
  det(2)    = real(isign, kind=dp)

 end function       array2_log_determinant_d
!-----------------------------------------------------------------------------------------------------------------------------------
! Use PFAPACK as a more efficient algorithm
 function                                                                  array2_pfaffian_pfa(C,alg)   result(pfaffian)            !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp)                                            :: pfaffian
  character(*), optional, intent(in)                                 :: alg
  character(1)                                           :: alg_
  complex(dp), allocatable                               :: C_(:,:)
  integer                                                :: n, info
  character(1000)                                        :: err

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   pfaffian = NaN()
   return
  end if
  n = size(C,1)
  if (MOD(n,2) /= 0) then
   pfaffian = ZERO
   return
  end if
  if (n == 0) then
   pfaffian = ONE
   return
  end if

  alg_ = 'P'
  if(present(alg))then
   if (len(alg) >= 1) alg_ = alg(1:1)
  end if
  
  select case(alg_)
  case ('C','c')
   pfaffian = array2_pfaffian(C)                               ! Catterall algorithm
  case default
   allocate  (C_,source=C)
   call skpfa(C_,pfaffian,info=info)                           ! PFAPACK   algorithm
   if (info /= 0) then
    write(err,*) 'array2_pfaffian_pfa: skpfa failed. INFO= ', info
    call matrix_error(trim(err))
   end if
  end select


 end function       array2_pfaffian_pfa
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_pfaffian_pfa                                (C,alg)                  result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp), dimension(2)                              :: pfaffian
  complex(dp), dimension(2)                              :: lpf
  character(*), optional, intent(in)                                 :: alg
  character(1)                                           :: alg_
  complex(dp), allocatable                               :: C_(:,:)
  integer                                                :: n, info
  character(1000)                                        :: err

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   pfaffian = NaN()
   return
  end if
  n = size(C,1)
  if (MOD(n,2) /= 0) then
   pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
   pfaffian(2) = ZERO
   return
  end if
  if (n == 0) then
   pfaffian(1) = ZERO
   pfaffian(2) = ONE
   return
  end if
       
  alg_ = 'P'
  if(present(alg))then
   if (len(alg) >= 1) alg_ = alg(1:1)
  end if
  
  select case(alg_)
  case ('C','c')
   pfaffian = array2_log_pfaffian(C)                               !Pf(C)= exp(pfaffian(1)) * pfaffian(2)    Catterall algorithm
  case default
   allocate   (C_,source=C)
   call skpf10(C_,lpf,info=info)  ! Pf(C) = lpf(1) * 10**lpf(2)
   if (info /= 0) then
    write(err,*) 'array2_log_pfaffian_pfa: skpf10 failed. INFO= ', info
    call matrix_error(trim(err))
   end if
   if (abs(lpf(1)) == 0.0_dp) then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if
   pfaffian(1) = cmplx(log(abs(lpf(1))) + real(lpf(2), kind=dp) * log(10.0_dp), 0.0_dp, kind=dp) !Pf(C)= exp(pfaffian(1)) * pfaffian(2) PFAPACK algorithm
   pfaffian(2) = exp(IMU * atan2(aimag(lpf(1)),real(lpf(1), kind=dp)))                           !           pfaffian(1) = log(|Pf(C)|)
  end select                                                                                     !           pfaffian(2) = exp( i Arg(Pf(C)) )  (phase of pfaffian e^{i θ} )


 end function array2_log_pfaffian_pfa
!-----------------------------------------------------------------------------------------------------------------------------------
! After C code by Simon Catterall,  Phys. Rev. D 68, 014503 (2003)  https://doi.org/10.1103/PhysRevD.68.014503
! The meaning of rows has been interchanged with columns in the comments because of the translation from C->Fortran for more efficient programming
!
! K. Anagnostopoulos, NTUA Feb 2007
!-----------------------------------------------------------------------------------------------------------------------------------
!Same as in array2_pfaffian2, but more aggressively avoiding outer index in inner loops.
!Checked for accuracy, array2_log_pfaffian is 3 times faster than array2_log_pfaffian2 for the intel compiler and 1200 < n < 2000
!                      array2_pfaffian     is 2 rimes faster than array2_pfaffian2     for the intel compiler and n = 312
!See ~/ikkt/codes/Pfaffian for original code.
!-----------------------------------------------------------------------------------------------------------------------------------
 function                                                                  array2_pfaffian(C)       result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp)                                            :: pfaffian
! complex(dp), dimension(size(C,1),size(C,2))            :: M
  complex(dp), dimension(:,:), allocatable               :: M        ! 260104
  integer                                                :: n
  integer                                                :: i,j,k,jpiv,interchange
  real   (dp)                                            :: pivot
! complex(dp),dimension(size(C,1))                       :: dum,sc
  complex(dp),dimension(:)   , allocatable               :: dum,sc   ! 260104
  complex(dp)                                            :: f      

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   pfaffian = NaN()
   return
  end if

  n = size(C,1)

  if(MOD(n,2) /= 0) then
   pfaffian = ZERO
   return
  end if
  if(n == 0) then
   pfaffian = ONE
   return
  end if

  allocate(M  (n,n))
  allocate(dum(n))
  allocate(sc (n))

  M = C

  interchange=1
      
  !Loop ovel all rows in steps of 2
  do i=1,N-1,2
   !first row i

   !find column whose ith component is biggest to use as pivot
   pivot=ABS(M(i+1,i))
   jpiv=i+1
   do j=i+2,N
    if(ABS(M(j,i))>pivot)then
     pivot=ABS(M(j,i))
     jpiv=j
    endif
   enddo

   if(pivot == 0.0_dp) then
    pfaffian = ZERO
    return
   end if

   !interchange col(i+1) with col(jpiv)
   dum=-M(:,i+1)
   do j=1,N
    M(i+1 ,j)=-M(j,jpiv)
    M(jpiv,j)=dum(j)
   enddo

   !interchange row(i+1) with row(jpiv)
   dum=M(:,i+1)
   do j=1,N
    M(j,i+1 )=M(j,jpiv)
    M(j,jpiv)=dum(j)
   enddo

   if(jpiv /= i+1) interchange = -interchange

   if(abs(M(i+1,i)) == 0.0_dp .or. abs(M(i,i+1)) == 0.0_dp) then
    pfaffian = ZERO
    return
   end if

   !using this,zero progressively elements of row M(j,i),j=i+2,N
   sc(i+2:) = M(i+2:,i)/M(i+1,i) 
   do k=1,N
    M(i+2:,k)=M(i+2:,k)+sc(i+2:)*M(k,i+1)
   enddo

   do j=i+2,N
    !zero out elements along corresponding column M(i,j) too
    M(:,j)=M(:,j)-sc(j)*M(:,i+1)
   enddo!do j=i+2,N

   !next row i+1
   !using this,zero progressively elements of row M(j,i),j=i+2,N
   sc(i+2:)=M(i+2:,i+1)/M(i,i+1)
   do k=1,N
    M(i+2:,k)=M(i+2:,k)+sc(i+2:)*M(k,i)
   end do

   !zero out elements along corresponding column too
   do j=i+2,N
    M(:,j)=M(:,j)-sc(j)*M(:,i)
   enddo!do j=i+2,N
       
  enddo !do i=1,N-1,2,Loop ovel all rows in steps of 2
  
  f=CMPLX(1.0_dp,0.0_dp,kind=dp)                                                                                                    ! 261125 (MINOR-04): was non-standard DCMPLX(1.0D0,0.0D0)
  do i=1,N,2
   f=f*M(i+1,i)
  enddo

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(N/2)
  !we have to change the sign appropriately:
  if(MOD(N/2,2).EQ.1) interchange=-interchange

  pfaffian=f*interchange

 end function       array2_pfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_pfaffian                                    (C)                      result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp), dimension(2)                              :: pfaffian
! complex(dp), dimension(size(C,1),size(C,2))            :: M
  complex(dp), dimension(:,:), allocatable               :: M       ! 260104
  integer                                                :: n
  integer                                                :: i,j,k,jpiv,interchange
  real   (dp)                                            :: pivot
! complex(dp),dimension(size(C,1))                       :: dum,sc
  complex(dp),dimension(:)   , allocatable               :: dum,sc  ! 260104
  real   (dp)                                            :: lnPfaffianAbs,AbsM
  complex(dp)                                            :: phase

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   pfaffian = NaN()
   return
  end if

  n = size(C,1)

  if(MOD(n,2) /= 0) then
   pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
   pfaffian(2) = ZERO
   return
  end if
  if(n == 0) then
   pfaffian(1) = ZERO
   pfaffian(2) = ONE
   return
  end if

  allocate(M  (n,n))
  allocate(dum(n))
  allocate(sc (n))

  M = C

  interchange=1
      
  !Loop ovel all rows in steps of 2
  do i=1,N-1,2
   !first row i

   !find column whose ith component is biggest to use as pivot
   pivot=ABS(M(i+1,i))
   jpiv=i+1
   do j=i+2,N
    if(ABS(M(j,i))>pivot)then
     pivot=ABS(M(j,i))
     jpiv=j
    endif
   enddo

   if(pivot == 0.0_dp) then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if

   !interchange col(i+1) with col(jpiv)
   dum=-M(:,i+1)
   do j=1,N
    M(i+1 ,j)=-M(j,jpiv)
    M(jpiv,j)=dum(j)
   enddo

   !interchange row(i+1) with row(jpiv)
   dum=M(:,i+1)
   do j=1,N
    M(j,i+1 )=M(j,jpiv)
    M(j,jpiv)=dum(j)
   enddo

   if(jpiv /= i+1) interchange = -interchange

   if(abs(M(i+1,i)) == 0.0_dp .or. abs(M(i,i+1)) == 0.0_dp) then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if

   !using this,zero progressively elements of row M(j,i),j=i+2,N
   sc(i+2:) = M(i+2:,i)/M(i+1,i) 
   do k=1,N
    M(i+2:,k)=M(i+2:,k)+sc(i+2:)*M(k,i+1)
   enddo

   do j=i+2,N
    !zero out elements along corresponding column M(i,j) too
    M(:,j)=M(:,j)-sc(j)*M(:,i+1)
   enddo!do j=i+2,N

   !next row i+1
   !using this,zero progressively elements of row M(j,i),j=i+2,N
   sc(i+2:)=M(i+2:,i+1)/M(i,i+1)
   do k=1,N
    M(i+2:,k)=M(i+2:,k)+sc(i+2:)*M(k,i)
   end do

   !zero out elements along corresponding column too
   do j=i+2,N
    M(:,j)=M(:,j)-sc(j)*M(:,i)
   enddo!do j=i+2,N
       
  enddo !do i=1,N-1,2,Loop ovel all rows in steps of 2
  
  !Final results: Compute Pfaffian:
  lnPfaffianAbs = 0.0D0
  phase         = ONE
  
  do i=1,n,2

   AbsM = ABS(M(i+1,i))

   if(AbsM <= 0.0_dp)then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if

   lnPfaffianAbs = lnPfaffianAbs + LOG(AbsM)
   phase         = phase         * (M(i+1,i)/AbsM)

  end do

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(n/2) we have to change the sign appropriately:
  if(MOD(n/2,2).EQ.1) interchange=-interchange

  phase       = phase * interchange

  pfaffian(1) = cmplx(lnPfaffianAbs, 0.0_dp, kind=dp)
  pfaffian(2) = phase

 end function array2_log_pfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
!Slower but easier coding. Keep it for reference. array2_log_pfaffian is approximately 3 times faster than array2_log_pfaffian2 (intel compiler)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_pfaffian2                                       (C)                      result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp)                                            :: pfaffian
! complex(dp), dimension(size(C,1),size(C,2))            :: M
  complex(dp), dimension(:,:), allocatable               :: M    ! 260104
  integer                                                :: n
  integer                                                :: i, j, jpiv, interchange
  real   (dp)                                            :: pivot
! complex(dp),dimension(size(C,1))                       :: dum
  complex(dp), dimension(:)  , allocatable               :: dum  ! 260104
  complex(dp)                                            :: scale,f

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   pfaffian = NaN()
   return
  end if

  n = size(C,1)

  if(MOD(n,2) /= 0) then
   pfaffian = ZERO
   return
  end if
  if(n == 0) then
   pfaffian = ONE
   return
  end if

  allocate(M  (n,n) )
  allocate(dum(n) )

  M = C

  interchange=1
      
  !Loop ovel all rows in steps of 2
  do i=1,n-1,2

   !first row i
   !find column whose ith component is biggest to use as pivot
   pivot=ABS(M(i+1,i))
   jpiv=i+1
   do j=i+2,n
    if(ABS(M(j,i)) > pivot)then
     pivot=ABS(M(j,i))
     jpiv=j
    endif
   enddo

   if(pivot == 0.0_dp) then
    pfaffian = ZERO
    return
   end if

   !interchange col(i+1) with col(jpiv)
   dum = M(i+1,:)
   do j=1,n
    M(i+1 ,j)=M(jpiv,j)
    M(jpiv,j)=dum(j)
   enddo

   !interchange row(i+1) with row(jpiv)
   dum = M(:,i+1)
   do j=1,n
    M(j,i+1 )=M(j,jpiv)
    M(j,jpiv)=dum(j)
   enddo

   if(jpiv /= i+1) interchange = -interchange

   if(abs(M(i+1,i)) == 0.0_dp .or. abs(M(i,i+1)) == 0.0_dp) then
    pfaffian = ZERO
    return
   end if

   !using this,zero progressively elements of row M(j,i),j=i+2,n
   do j=i+2,n
    scale=M(j,i)/M(i+1,i)
    M(j,:)=M(j,:)-scale*M(i+1,:)
    !zero out elements along corresponding column M(i,j) too
    M(:,j)=M(:,j)-scale*M(:,i+1)
   enddo!do j=i+2,n

   !next row i+1
   !using this,zero progressively elements of row M(j,i),j=i+2,n
   do j=i+2,n

    scale=M(j,i+1)/M(i,i+1)
    M(j,:)=M(j,:)-scale*M(i,:)

    !zero out elements along corresponding column too
    M(:,j)=M(:,j)-scale*M(:,i)

   enddo!do j=i+2,n
       
  enddo !do i=1,n-1,2,Loop ovel all rows in steps of 2

  f=CMPLX(1.0_dp,0.0_dp,kind=dp)
  do i=1,n,2
   f=f*M(i+1,i)
  enddo

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(n/2) we have to change the sign appropriately:
  if(MOD(n/2,2).EQ.1) interchange=-interchange

  pfaffian=f*interchange

 end function       array2_pfaffian2
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_pfaffian2                                   (C)                      result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), contiguous, intent(in)    :: C
  complex(dp), dimension(2)                              :: pfaffian
! complex(dp), dimension(size(C,1),size(C,2))            :: M
  complex(dp), dimension(:,:), allocatable               :: M    ! 260104
  integer                                                :: n
  integer                                                :: i, j, jpiv, interchange
  real   (dp)                                            :: pivot
! complex(dp),dimension(size(C,1))                       :: dum
  complex(dp), dimension(:)  , allocatable               :: dum  ! 260104
  complex(dp)                                            :: scale  
  real   (dp)                                            :: lnPfaffianAbs,AbsM
  complex(dp)                                            :: phase

  if(size(C,1)/=size(C,2))then ! 261124: non-square input -> NaN (was: matrix_error)
   pfaffian = NaN()
   return
  end if

  n = size(C,1)

  if(MOD(n,2) /= 0) then
   pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
   pfaffian(2) = ZERO
   return
  end if
  if(n == 0) then
   pfaffian(1) = ZERO
   pfaffian(2) = ONE
   return
  end if

  allocate(M  (n,n))
  allocate(dum(n))

  M = C

  interchange=1
      
  !Loop ovel all rows in steps of 2
  do i=1,n-1,2

   !first row i
   !find column whose ith component is biggest to use as pivot
   pivot=ABS(M(i+1,i))
   jpiv=i+1
   do j=i+2,n
    if(ABS(M(j,i)) > pivot)then
     pivot=ABS(M(j,i))
     jpiv=j
    endif
   enddo

   if(pivot == 0.0_dp) then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if

   !interchange col(i+1) with col(jpiv)
   dum = M(i+1,:)
   do j=1,n
    M(i+1 ,j)=M(jpiv,j)
    M(jpiv,j)=dum(j)
   enddo

   !interchange row(i+1) with row(jpiv)
   dum = M(:,i+1)
   do j=1,n
    M(j,i+1 )=M(j,jpiv)
    M(j,jpiv)=dum(j)
   enddo

   if(jpiv /= i+1) interchange = -interchange

   if(abs(M(i+1,i)) == 0.0_dp .or. abs(M(i,i+1)) == 0.0_dp) then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if

   !using this,zero progressively elements of row M(j,i),j=i+2,n
   do j=i+2,n
    scale=M(j,i)/M(i+1,i)
    M(j,:)=M(j,:)-scale*M(i+1,:)
    !zero out elements along corresponding column M(i,j) too
    M(:,j)=M(:,j)-scale*M(:,i+1)
   enddo!do j=i+2,n

   !next row i+1
   !using this,zero progressively elements of row M(j,i),j=i+2,n
   do j=i+2,n

    scale=M(j,i+1)/M(i,i+1)
    M(j,:)=M(j,:)-scale*M(i,:)

    !zero out elements along corresponding column too
    M(:,j)=M(:,j)-scale*M(:,i)

   enddo!do j=i+2,n
       
  enddo !do i=1,n-1,2,Loop ovel all rows in steps of 2

  !Final results: Compute Pfaffian:
  lnPfaffianAbs = 0.0D0
  phase         = ONE
  
  do i=1,n,2

   AbsM = ABS(M(i+1,i))

   if(AbsM <= 0.0_dp)then
    pfaffian(1) = cmplx(IEEE_VALUE(0.0_dp, IEEE_NEGATIVE_INF), 0.0_dp, kind=dp)
    pfaffian(2) = ZERO
    return
   end if

   lnPfaffianAbs = lnPfaffianAbs + LOG(AbsM)
   phase         = phase         * (M(i+1,i)/AbsM)

  end do

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(n/2) we have to change the sign appropriately:
  if(MOD(n/2,2).EQ.1) interchange=-interchange

  phase       = phase * interchange

  pfaffian(1) = cmplx(lnPfaffianAbs, 0.0_dp, kind=dp)
  pfaffian(2) = phase
  
 end function       array2_log_pfaffian2
 

!-----------------------------------------------------------------------------------------------------------------------------------
!...................................................................................................................................!:.:.:
!-----------------------------------------------------------------------------------------------------------------------------------



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
