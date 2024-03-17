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
  complex  (dp), dimension(:,:), intent(in)       :: A,B
  character(*) , optional      , intent(in)       :: mtype
  complex  (dp), dimension(size(A,1),size(B,2))   :: C
  integer                                         :: ma,na,mb,nb,mc,nc
  complex  (dp)                                   :: alpha, beta
  character(1)                                    :: opA, opB, side,uplo
  character(mtype_len)                            :: mt
! logical      , save                             :: have_I_not_warned_256 = .true.

  mt = 'ZG'; if(present(mtype)) mt = mtype

  ma = size(A,1) ; na = size(A,2)
  mb = size(B,1) ; nb = size(B,2)
  mc = size(A,1) ; nc = size(B,2)

!  if(have_I_not_warned_256 .and. (MOD(ma,256) == 0 .or. MOD(mb,256) == 0 .or. MOD(mc,256) == 0) )then
   !From here: https://software.intel.com/en-us/articles/a-simple-example-to-measure-the-performance-of-an-intel-mkl-function
!   have_I_not_warned_256 = .false.
!   write(f_mout,'(A)') '# WARNING: array2_matmul_lapack_mm: One of the leading dimensions is a multiple of 256. Performance issues?'
!  end if

!  if( na /= mb ) call matrix_error('array2_matmul_lapack_mm: na /= mb')

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

  select case(mt(2:2))
  case('H','h')
    
!  if( ma /= na ) call matrix_error('array2_matmul_lapack_mm: ma /= na')
   if( ma /= na ) then
    C = NaN()
    return
   end if

   call zhemm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Hermitian                                               !:.:.:

  case default
    
   call zgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                                                         !:.:.:

  end select

 end function       array2_matmul_lapack_mm
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_dd                                (A,B,mtype)              result(C)                       !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A,B
  character(*) , optional      , intent(in)       :: mtype
  real     (dp), dimension(size(A,1),size(B,2))   :: C
  integer                                         :: ma,na,mb,nb,mc,nc
  real   (dp)                                     :: alpha, beta
  character(1)                                    :: opA, opB, side,uplo
  character(mtype_len)                            :: mt
! logical      , save                             :: have_I_not_warned_256 = .true.

  mt = 'DG'; if(present(mtype)) mt = mtype

  ma = size(A,1) ; na = size(A,2)
  mb = size(B,1) ; nb = size(B,2)
  mc = size(A,1) ; nc = size(B,2)

  !if( na /= mb )  call matrix_error('array2_matmul_lapack_dd: na /= mb')
  if( na /= mb )then
   C = NaN()
   return
  end if

  !if(have_I_not_warned_256 .and. (MOD(ma,256) == 0 .or. MOD(mb,256) == 0 .or. MOD(mc,256) == 0) )then
  !From here: https://software.intel.com/en-us/articles/a-simple-example-to-measure-the-performance-of-an-intel-mkl-function
  ! have_I_not_warned_256 = .false.
  ! write(f_mout,'(A)') '# WARNING: array2_matmul_lapack_dd: One of the leading dimensions is a multiple of 256. Performance issues?'
  !end if

  ! Lapack calculates C = alpha * A . B + beta * C
  !   C    =   A   x   B
  ! m x n    m x k   k x n
  ! call dgemm('N','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) where 'N'  is for normal, 'T' is for transpose
  ! If A is Hermitian, then we can call:
  ! call dsymm('L','U',M,N,  alpha,A,lda,B,ldb,beta,C,ldc) where side= 'L' is for A.B, side= 'R' for B.A and uplo='U'/'L' references the upper/lower part of the A matrix

  alpha = 1.0_dp ; beta = 0.0_dp
  opA   = 'N' ; opB  = 'N'  ; side ='L' ; uplo = 'U'

  select case(mt(2:2))
  case('S','s')
    
  !if( ma /= na ) call matrix_error('array2_matmul_lapack_dd: ma /= na')
   if( ma /= na )then
    C = NaN()
    return
   end if

   call dsymm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Symmetric                                               !:.:.:

  case default
   
   call dgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                                                         !:.:.:

  end select

 end function       array2_matmul_lapack_dd
!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^\dagger = array2_matmul_lapack_m(A)    ;   A^\dagger . A = array2_matmul_lapack_m(A,side='H')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_m                                 (A,side)                 result(C)                       !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: A
  character(*) , optional      , intent(in)       :: side
  complex  (dp), dimension(:,:), allocatable      :: C
  integer                                         :: ma,na, mc,nc
  complex(dp)                                     :: alpha, beta
  character(1)                                    :: opA, uplo
  integer                                         :: N,K,i,j

  ma  = size(A,1) ; na = size(A,2)

  opA = 'N';                       ! default is A . A^\dagger
  N   =  ma       ; K  = na

  if(present(side))then
   select case(side(1:1))
   case('C','c','T','t','H','h','S','s','L','l')

    opA = 'C'
    N   = na      ; K  = ma

   end select
  end if

  ! Lapack calculates C = alpha * A . A^\dagger + beta * C
  ! if opA='N' then   C    =   A   x   A^\dagger          ::          if opA='T'    C    =   A^\dagger    x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m          m x n
  !
  ! call zherk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'

  alpha = ONE ; beta = ZERO; uplo = 'U'
  allocate(C(N,N))
  call zherk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                                                    !:.:.:

  !only the upper triangle part of C is filled:
  do concurrent(i=1:n, j=1:n, i>j)
   C(i,j) = CONJG(C(j,i))
  end do

 end function       array2_matmul_lapack_m

!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^T = array2_matmul_lapack_m(A)    ;   A^T . A = array2_matmul_lapack_m(A,side='S')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_d                                 (A,side)                 result(C)                       !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A
  character(*) , optional      , intent(in)       :: side
  real     (dp), dimension(:,:), allocatable      :: C
  integer                                         :: ma,na, mc,nc
  real   (dp)                                     :: alpha, beta
  character(1)                                    :: opA, uplo
  integer                                         :: N,K,i,j

  ma  = size(A,1) ; na = size(A,2)

  opA = 'N';                       ! default is A . A^T
  N   =  ma       ; K  = na

  if(present(side))then
   select case(side(1:1))
   case('C','c','T','t','H','h','S','s','L','l')

    opA = 'T'
    N   = na      ; K  = ma

   end select
  end if


  ! Lapack calculates C = alpha * A . A^T + beta * C
  ! if opA='N' then   C    =   A   x   A^T                ::          if opA='T'    C    =     A^T   x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m     m x n
  !
  ! call dsyrk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'


  alpha = 1.0_dp ; beta = 0.0_dp; uplo = 'U'
  allocate(C(N,N))
  call dsyrk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                                                    !:.:.:
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
!type='T'    w = v . A^T         w_i = v_j   A_ji
!type='H'    w = A . v           w_i = A_ij  v_j    faster, when A = A^\dagger (no check performed)
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_mv                                (A,v,type)               result(w)                       !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: A
  complex  (dp), dimension(:)  , intent(in)       :: v
  character(*) , optional      , intent(in)       :: type
  complex  (dp), dimension(:)  , allocatable      :: w
  integer                                         :: ma,na,nv
  complex  (dp)                                   :: alpha, beta
  character(1)                                    :: tp,uplo
  integer                                         :: incx, incy

  tp = 'N'                                        ! default is w = A.v
  if(present(type)) tp  = type(1:1)

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1)

  alpha = ONE    ; beta = ZERO
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    allocate(w(ma))

   !if( na /= nv) call matrix_error('array2_matmul_lapack_mv: na /= nv (N)'   )
    if( na /= nv)then
     w = NaN()
     return
    end if

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                              !:.:.:

   case('T','C','t','c')

    allocate(w(na))

   !if( ma /= nv) call matrix_error('array2_matmul_lapack_mv: ma /= nv (T,C)' )
    if( ma /= nv)then
     w = NaN()
     return
    end if

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)

   case('H','h')

    uplo = 'U'
    allocate(w(ma))

   !if( ma /= na) call matrix_error('array2_matmul_lapack_mv: ma /= na (H)'   )
   !if( na /= nv) call matrix_error('array2_matmul_lapack_mv: na /= nv (H)'   )
    if( ma /= na .or. na /= nv)then
     w = NaN()
     return
    end if

    call zhemv(uplo ,ma,alpha,A,ma,v,incx,beta,w,incy)                                                                              !:.:.:

   case default

    allocate(w(ma))
    w = NaN()
    return
   !call               matrix_error('array2_matmul_lapack_mv: wrong type'     )

  end select


 end function       array2_matmul_lapack_mv

!-----------------------------------------------------------------------------------------------------------------------------------
!type='N'    w = A . v           w_i = A_ij  v_j
!type='T'    w = v . A^T         w_i = v_j   A_ji
!type='S'    w = v . A^T         w_i = v_j   A_ji
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_lapack_dv                                (A,v,type)               result(w)                       !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A
  real     (dp), dimension(:)  , intent(in)       :: v
  character(*) , optional      , intent(in)       :: type
  real     (dp), dimension(:)  , allocatable      :: w
  integer                                         :: ma,na,nv
  real     (dp)                                   :: alpha, beta
  character(1)                                    :: tp
  integer                                         :: incx, incy

  tp = 'N'                                        ! default is w = A.v
  if(present(type)) tp  = type(1:1)
  if(tp == 'S'    ) tp  = 'T'                     ! allow S for Symmetric (zgemv accepts only 'T')

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1)

  alpha = 1.0_dp ; beta = 0.0_dp
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    allocate(w(ma))

   !if( na /= nv) call matrix_error('array2_matmul_lapack_dv: na /= nv')
    if( na /= nv)then
     w = NaN()
     return
    end if

   case('T','C','t','c')

    allocate(w(na))

   !if( ma /= nv) call matrix_error('array2_matmul_lapack_dv: ma /= nv')
    if( ma /= nv)then
     w = NaN()
     return
    end if

   case default

   !call               matrix_error('array2_matmul_lapack_dv: wrong type')
    allocate(w(ma))
    w = NaN()
    return

  end select

  call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                                !:.:.:

 end function       array2_matmul_lapack_dv



!...................................................................................................................................
!.......................... Matrix Matrix Multiplication Subroutines ...............................................................!:.:.:
!...................................................................................................................................
!
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_mm_sub                            (A,B,C,mtype)                                            !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: A,B
  complex  (dp), dimension(:,:), intent(inout)    :: C
  character(*) , optional      , intent(in)       :: mtype
  integer                                         :: ma,na,mb,nb,mc,nc
  complex  (dp)                                   :: alpha, beta
  character(1)                                    :: opA, opB, side,uplo
  character(mtype_len)                            :: mt

  mt = 'ZG'; if(present(mtype)) mt = mtype

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

  select case(mt(2:2))
  case('H','h')
    
   if( ma /= na ) then
    C = NaN()
    return
   end if

   call zhemm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Hermitian                                               !:.:.:

  case default

   call zgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                                                         !:.:.:

  end select

 end  subroutine    array2_matmul_lapack_mm_sub
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_dd_sub                            (A,B,C,mtype)                                            !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A,B
  real     (dp), dimension(:,:), intent(inout)    :: C
  character(*) , optional      , intent(in)       :: mtype
  integer                                         :: ma,na,mb,nb,mc,nc
  real   (dp)                                     :: alpha, beta
  character(1)                                    :: opA, opB, side,uplo
  character(mtype_len)                            :: mt
! logical      , save                             :: have_I_not_warned_256 = .true.

  mt = 'DG'; if(present(mtype)) mt = mtype

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
  ! call dgemm('N','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc) where 'N'  is for normal, 'T' is for transpose
  ! If A is Hermitian, then we can call:
  ! call dsymm('L','U',M,N,  alpha,A,lda,B,ldb,beta,C,ldc) where side= 'L' is for A.B, side= 'R' for B.A and uplo='U'/'L' references the upper/lower part of the A matrix

  alpha = 1.0_dp ; beta = 0.0_dp
  opA   = 'N' ; opB  = 'N'  ; side ='L' ; uplo = 'U'

  select case(mt(2:2))
  case('S','s')

   if( ma /= na )then
    C = NaN()
    return
   end if

   call dsymm(side,uplo,mc,nc   ,alpha,A,ma,B,mb,beta,C,mc)  ! A is assumed Symmetric                                               !:.:.:

  case default
   
   call dgemm(opA ,opB ,mc,nc,na,alpha,A,ma,B,mb,beta,C,mc)                                                                         !:.:.:

  end select

 end  subroutine    array2_matmul_lapack_dd_sub
!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^\dagger = array2_matmul_lapack_m(A)    ;   A^\dagger . A = array2_matmul_lapack_m(A,side='H')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_m_sub                             (A,C,side)                                               !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: A
  complex  (dp), dimension(:,:), intent(inout)    :: C
  character(*) , optional      , intent(in)       :: side
  integer                                         :: ma,na, mc,nc
  complex(dp)                                     :: alpha, beta
  character(1)                                    :: opA, uplo
  integer                                         :: N,K,i,j

  ma  = size(A,1) ; na = size(A,2)
  mc  = size(C,1) ; nc = size(C,2)

  opA = 'N';                       ! default is A . A^\dagger
  N   =  ma       ; K  = na

  if(present(side))then
   select case(side(1:1))
   case('C','c','T','t','H','h','S','s','L','l')

    opA = 'C'
    N   = na      ; K  = ma

   end select
  end if

  ! Lapack calculates C = alpha * A . A^\dagger + beta * C
  ! if opA='N' then   C    =   A   x   A^\dagger          ::          if opA='T'    C    =   A^\dagger    x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m          m x n
  !
  ! call zherk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'

  alpha = ONE ; beta = ZERO; uplo = 'U'

  if( mc /= N .or. nc /= N)then
   C = NaN()
   return
  end if

  call zherk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                                                    !:.:.:

  !only the upper triangle part of C is filled:
  do concurrent(i=1:n, j=1:n, i>j)
   C(i,j) = CONJG(C(j,i))
  end do

 end  subroutine    array2_matmul_lapack_m_sub
!-----------------------------------------------------------------------------------------------------------------------------------
! A . A^T = array2_matmul_lapack_m(A)    ;   A^T . A = array2_matmul_lapack_m(A,side='S')
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_matmul_lapack_d_sub                             (A,C,side)                                               !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A
  real     (dp), dimension(:,:), intent(inout)    :: C
  character(*) , optional      , intent(in)       :: side
  integer                                         :: ma,na, mc,nc
  real     (dp)                                   :: alpha, beta
  character(1)                                    :: opA, uplo
  integer                                         :: N,K,i,j

  ma  = size(A,1) ; na = size(A,2)
  mc  = size(C,1) ; nc = size(C,2)

  opA = 'N';                       ! default is A . A^\dagger
  N   =  ma       ; K  = na

  if(present(side))then
   select case(side(1:1))
   case('C','c','T','t','H','h','S','s','L','l')

    opA = 'C'
    N   = na      ; K  = ma

   end select
  end if

  ! Lapack calculates C = alpha * A . A^\dagger + beta * C
  ! if opA='N' then   C    =   A   x   A^\dagger          ::          if opA='T'    C    =   A^\dagger    x   A    
  !                 m x m    m x n   n x m                                        n x n      n x m          m x n
  !
  ! call dsyrk('U','N',N,K,alpha,A,lda,beta,C,ldc)  N=m, K=n for opA='N' and N=n, K=m for opA='T'

  alpha = 1.0_dp ; beta = 0.0_dp; uplo = 'U'

  if( mc /= N .or. nc /= N)then
   C = NaN()
   return
  end if

  call dsyrk(uplo, opA ,N,K,alpha,A,ma,beta,C,N)                                                                                    !:.:.:

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
!type='T'    w = v . A^T         w_i = v_j   A_ji
!type='H'    w = A . v           w_i = A_ij  v_j    faster, when A = A^\dagger (no check performed)
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_mv_sub                            (A,v,w,type)                                             !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: A
  complex  (dp), dimension(:)  , intent(in)       :: v
  complex  (dp), dimension(:)  , intent(inout)    :: w
  character(*) , optional      , intent(in)       :: type
  integer                                         :: ma,na,nv,nw
  complex  (dp)                                   :: alpha, beta
  character(1)                                    :: tp,uplo
  integer                                         :: incx, incy

  tp = 'N'                                        ! default is w = A.v
  if(present(type)) tp  = type(1:1)

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

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                              !:.:.:

   case('T','C','t','c')

    if( ma /= nv .or. nw /= na )then
     w = NaN()
     return
    end if

    call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)

   case('H','h')

    uplo = 'U'

    if( ma /= na .or. na /= nv .or. nw /= ma )then
     w = NaN()
     return
    end if

    call zhemv(uplo ,ma,alpha,A,ma,v,incx,beta,w,incy)                                                                              !:.:.:

   case default

    w = NaN()
    return

  end select


 end  subroutine    array2_matmul_lapack_mv_sub
!-----------------------------------------------------------------------------------------------------------------------------------
!type='T'    w = v . A^T         w_i = v_j   A_ji  
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_vm_sub                            (v,A,w)                                                  !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: A
  complex  (dp), dimension(:)  , intent(in)       :: v
  complex  (dp), dimension(:)  , intent(inout)    :: w
  integer                                         :: ma,na,nv,nw
  complex  (dp)                                   :: alpha, beta
  character(1)                                    :: tp,uplo
  integer                                         :: incx, incy

  tp = 'T'                                        

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = ONE    ; beta = ZERO
  incx  = 1      ; incy = 1

  if( ma /= nv .or. nw /= na )then
   w = NaN()
   return
  end if

  call zgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                                !:.:.:


 end  subroutine    array2_matmul_lapack_vm_sub
!-----------------------------------------------------------------------------------------------------------------------------------
!type='N'    w = A . v           w_i = A_ij  v_j
!type='T'    w = v . A^T         w_i = v_j   A_ji
!type='S'    w = v . A^T         w_i = v_j   A_ji
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_dv_sub                            (A,v,w,type)                                             !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A
  real     (dp), dimension(:)  , intent(in)       :: v
  real     (dp), dimension(:)  , intent(inout)    :: w
  character(*) , optional      , intent(in)       :: type
  integer                                         :: ma,na,nv,nw
  real     (dp)                                   :: alpha, beta
  character(1)                                    :: tp,uplo
  integer                                         :: incx, incy

  tp = 'N'                                        ! default is w = A.v
  if(present(type)) tp  = type(1:1)
  if(tp == 'S'    ) tp  = 'T'                     ! allow S for Symmetric (zgemv accepts only 'T')

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = 1.0_dp ; beta = 0.0_dp
  incx  = 1      ; incy = 1

  select case (tp)
   case('N','n')

    if( na /= nv .or. nw /= ma )then
     w = NaN()
     return
    end if

   case('T','C','t','c')

    if( ma /= nv .or. nw /= na )then
     w = NaN()
     return
    end if

   case('H','h')

    uplo = 'U'

    if( ma /= na .or. na /= nv .or. nw /= ma )then
     w = NaN()
     return
    end if

   case default

    w = NaN()
    return

  end select

  call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                                !:.:.:

 end  subroutine    array2_matmul_lapack_dv_sub
!-----------------------------------------------------------------------------------------------------------------------------------
!type='T'    w = v . A^T         w_i = v_j   A_ji  
!-----------------------------------------------------------------------------------------------------------------------------------
pure subroutine     array2_matmul_lapack_vd_sub                            (v,A,w)                                                  !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: A
  real     (dp), dimension(:)  , intent(in)       :: v
  real     (dp), dimension(:)  , intent(inout)    :: w
  integer                                         :: ma,na,nv,nw
  real     (dp)                                   :: alpha, beta
  character(1)                                    :: tp,uplo
  integer                                         :: incx, incy

  tp = 'T'                                        

  ma = size(A,1) ;  na  = size(A,2)
  nv = size(v,1) ;  nw  = size(w,1)

  alpha = 1.0_dp ; beta = 0.0_dp
  incx  = 1      ; incy = 1

  if( ma /= nv .or. nw /= na )then
   w = NaN()
   return
  end if

  call dgemv(tp,ma,na,alpha,A,ma,v,incx,beta,w,incy)                                                                                !:.:.:


 end  subroutine    array2_matmul_lapack_vd_sub



!...................................................................................................................................
!.......................... Matrix Inversion .......................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_inverse                                         (C)                      result(CI)                      !:.:.:
  complex  (dp), dimension(:,:), intent(in)       :: C
  complex  (dp), dimension(size(C,1),size(C,2))   :: CI
  integer                                         :: n, LWORK, info
  integer      , dimension(size(C,1))             :: ipiv
  complex(dp)  , dimension(:)       ,allocatable  :: WORK
  complex(dp)  , dimension(1)                     :: WORKTEST
  character(1024)                                 :: err

  if( size(C,1)/=size(C,2)) call matrix_error('array2_inverse: Input matrix not a square matrix. m /= n')

  n = size(C,1)
  !---------------------------------------------------------
  !Use zgetrf to compute LU decomposition: MATIV = P * L * U
  CI = C
  call zgetrf(n,n,CI,n,ipiv,info)                                                                                                   !:.:.:
  if(info /= 0) write(err,*)'array2_inverse:  zgetrf failed to perform LU decomposition. INFO= ',info
  if(info /= 0) call matrix_error(err)
  !---------------------------------------------------------
  !Use zgetri to compute inverse:
  LWORK = -1 ! compute optimal LWORK
  call zgetri(n,CI,n,ipiv,WORKTEST,LWORK,info)
  if(info /= 0) write(err,*)'array2_inverse: zgetri failed to determine LWORK. INFO= ',info
  if(info /= 0) call matrix_error(err)
  LWORK = WORKTEST(1)
  !---------------------------------------------------------
  !Compute inverse:
  ALLOCATE(WORK(LWORK))
  call zgetri(n,CI,n,ipiv,WORK    ,LWORK,info)                                                                                      !:.:.:
  if(info /= 0) write(err,*)'array2_inverse: zgetri failed to compute inverse. INFO= ',info
  if(info /= 0) call matrix_error(err)
  !---------------------------------------------------------
  DEALLOCATE(WORK)

 end function array2_inverse
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_inverse_d                                       (C)                      result(CI)                      !:.:.:
  real     (dp), dimension(:,:), intent(in)       :: C
  real     (dp), dimension(size(C,1),size(C,2))   :: CI
  integer                                         :: n, LWORK, info
  integer      , dimension(size(C,1))             :: ipiv
  real   (dp)  , dimension(:)       ,allocatable  :: WORK
  real   (dp)  , dimension(1)                     :: WORKTEST
  character(1024)                                 :: err

  if( size(C,1)/=size(C,2)) call matrix_error('array2_inverse: Input matrix not a square matrix. m /= n')

  n = size(C,1)
  !---------------------------------------------------------
  !Use dgetrf to compute LU decomposition: MATIV = P * L * U
  CI = C
  call dgetrf(n,n,CI,n,ipiv,info)                                                                                                   !:.:.:
  if(info /= 0) write(err,*)'array2_inverse_d:  dgetrf failed to perform LU decomposition. INFO= ',info
  if(info /= 0) call matrix_error(err)
  !---------------------------------------------------------
  !Use dgetri to compute inverse:
  LWORK = -1 ! compute optimal LWORK
  call dgetri(n,CI,n,ipiv,WORKTEST,LWORK,info)
  if(info /= 0) write(err,*)'array2_inverse_d: dgetri failed to determine LWORK. INFO= ',info
  if(info /= 0) call matrix_error(err)
  LWORK = WORKTEST(1) + 10
  !---------------------------------------------------------
  !Compute inverse:
  ALLOCATE(WORK(LWORK))
  call dgetri(n,CI,n,ipiv,WORK    ,LWORK,info)                                                                                      !:.:.:
  if(info /= 0) write(err,*)'array2_inverse: dgetri failed to compute inverse. INFO= ',info
  if(info /= 0) call matrix_error(err)
  !---------------------------------------------------------
  DEALLOCATE(WORK)

 end function array2_inverse_d
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... Eigenvalues - Eigenvectors .............................................................................!:.:.:
!...................................................................................................................................
! Notice that for large enough matrices, the order of the eigenvalues for zgeev and dgeev don't come with the same order when job='N'/'V'
! In order to compare, you have to sort them in e.g. real part order and then the difference is pure numerical error. 
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvalues                                     (C,mtype)                result(eigenval)                !:.:.:
  complex  (dp), dimension(:,:)             , intent(in)                   :: C
  character(*)                              , optional                     :: mtype
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec
  character(mtype_len)                                                     :: mt,job

  mt  = 'ZG'; if(present(mtype)) mt = mtype
  job = 'N' !only eigenvalues

  select case(mt(2:2))
  case('H','h') !hermitian matrix
   call array2_zheev(C,eigenval,eigenvec,job)
  case default  !general   matrix
   call array2_zgeev(C,eigenval,eigenvec,job)
  end select

 end function       array2_eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvectors                                    (C,mtype,sortby)         result(evs)                     !:.:.:
  complex  (dp), dimension(:,:)             , intent(in)                   :: C
  character(*)                              , optional                     :: mtype, sortby
  integer      , dimension(:)               , allocatable                  :: pos
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec, eigenvectmp
  complex  (dp), dimension(:,:)             , allocatable                  :: evs
  character(mtype_len)                                                     :: mt,job
  integer                                                                  :: n,i

  n   = size(C,1)

  mt  = 'ZG'; if(present(mtype)) mt = mtype
  job = 'V'     ! eigenvectors too

  select case(mt(2:2))
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

  evs(:,1 ) = eigenval         ! The first column has the eigenvalues
  evs(:,2:) = eigenvec         ! Columns 2-(n+1)  has the eigenvectors

 end function       array2_eigenvectors
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvalues_d                                   (C,mtype)                result(eigenval)                !:.:.:
  real     (dp), dimension(:,:)             , intent(in)                   :: C
  character(*)                              , optional                     :: mtype
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec
  character(mtype_len)                                                     :: mt,job

  mt  = 'DG'; if(present(mtype)) mt = mtype
  job = 'N'     ! only eigenvalues

  select case(mt(2:2))
  case('S','s') ! symmetric matrix
   call array2_dsyev(C,eigenval,eigenvec,job)
  case default  ! general   matrix
   call array2_dgeev(C,eigenval,eigenvec,job)
  end select

 end function       array2_eigenvalues_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_eigenvectors_d                                  (C,mtype,sortby)         result(evs)                     !:.:.:
  real     (dp), dimension(:,:)             , intent(in)                   :: C
  character(*)                              , optional                     :: mtype, sortby
  integer      , dimension(:)               , allocatable                  :: pos
  complex  (dp), dimension(:)               , allocatable                  :: eigenval
  complex  (dp), dimension(:,:)             , allocatable                  :: eigenvec, eigenvectmp
  complex  (dp), dimension(:,:)             , allocatable                  :: evs
  character(mtype_len)                                                     :: mt,job
  integer                                                                  :: n,i

  n   = size(C,1)

  mt  = 'DG'; if(present(mtype)) mt = mtype
  job = 'V'     ! eigenvectors too

  select case(mt(2:2))
  case('S','s') ! hermitian matrix
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

  evs(:,1 ) = eigenval        ! The first column has the eigenvalues
  evs(:,2:) = eigenvec        ! Columns 2-(n+1)  has the eigenvectors. eigenvec is real, evs is complex

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
  complex  (dp), dimension(:,:)             , intent(in)                   :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
  complex  (dp), dimension(size(C,1),size(C,2))                            :: A
  integer                                                                  :: n
  character(1)                                                             :: JOBVL, JOBVR
  integer                                                                  :: info,  LWORK, LDVL,LDVR
  complex  (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: RWORK
  complex  (dp), dimension(:,:)                             , allocatable  :: VL,VR
  complex  (dp), dimension(:)                               , allocatable  :: EV,WORK
  character(1024)                                                          :: err
  
  n =size(C,1)

  if(size(C,1) /= size(C,2)   )  call matrix_error('array2_zgeev: Input matrix must be a square matrix:  m /= n')

  if(allocated(eigenval)      )  deallocate(eigenval) ; allocate(eigenval(n))

  A = C  ! A local copy of the matrix C

  allocate(RWORK(2*n)); allocate(EV(n))

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(job(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBVL = 'N' ; JOBVR = 'N'; LDVL = 1; LDVR = 1; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,RWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVL,WORK    ,LWORK,RWORK,info)                                                       !:.:.:
   if(info /= 0) write(err,*)                      'array2_zgeev: zgeev failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = EV
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   if(allocated(eigenvec)) deallocate(eigenvec); allocate(eigenvec(n,n))
   JOBVL = 'N' ; JOBVR = 'V'; LDVL = 1; LDVR = n; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,RWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call zgeev(JOBVL,JOBVR,n,A,n,EV,VL,LDVL,VR,LDVR,WORK    ,LWORK,RWORK,info)
   if(info /= 0) write(err,*)                      'array2_zgeev: zgeev failed (V). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval  = EV
   eigenvec = VR
  !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call                               matrix_error('array2_zgeev: wrong value of job: '//trim(job))
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
  complex  (dp), dimension(:,:)             , intent(in)                   :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
  complex  (dp), dimension(size(C,1),size(C,2))                            :: A
  integer                                                                  :: n
  integer                                                                  :: info,  LWORK
  complex  (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: RWORK
  real     (dp), dimension(:)                               , allocatable  :: EV
  complex  (dp), dimension(:)                               , allocatable  :: WORK
  character(1024)                                                          :: err
  character(1)                                                             :: JOBZ, UPLO

  n = size(C,1)

  if( size(C,1) /= size(C,2)  )  call matrix_error('array2_zheev: Input matrix must be a square matrix:  m /= n')

  if(allocated(eigenval)      )  deallocate(eigenval) ; allocate(eigenval(n))

  A = C  ! A local copy of the matrix C

  allocate(RWORK(3*n-2)); allocate(EV(n))

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(job(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBZ='N' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call zheev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,RWORK,info)                                                   
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call zheev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,RWORK,info)                                                                         !:.:.:
   if(info /= 0) write(err,*)                      'array2_zheev: zheev failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = EV
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   if(allocated(eigenvec)) deallocate(eigenvec); allocate(eigenvec(n,n))
   JOBZ='V' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call zheev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,RWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call zheev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,RWORK,info)
   if(info /= 0) write(err,*)                      'array2_zheev: zheev failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = EV
   eigenvec = A
  !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call                               matrix_error('array2_zheev: wrong value of job: '//trim(job))
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(RWORK) ; deallocate(WORK)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_zheev
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_dgeev                                           (C,eigenval,eigenvec,job)                                !:.:.:
  real     (dp), dimension(:,:)             , intent(in)                   :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
  real     (dp), dimension(size(C,1),size(C,2))                            :: A
  integer                                                                  :: n,j
  character(1)                                                             :: JOBVL, JOBVR
  integer                                                                  :: info,  LWORK, LDVL,LDVR
  real     (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: RWORK
  real     (dp), dimension(:,:)                             , allocatable  :: VL,VR
  real     (dp), dimension(:)                               , allocatable  :: ReEV,ImEV,WORK
  character(1024)                                                          :: err
  
  n =size(C,1)

  if(size(C,1) /= size(C,2)   )  call matrix_error('array2_dgeev: Input matrix must be a square matrix:  m /= n')

  if(allocated(eigenval)      )  deallocate(eigenval) ; allocate(eigenval(n))

  A = C  ! A local copy of the matrix C

  allocate(ReEV(n)); allocate(ImEV(n)); 

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(job(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBVL = 'N' ; JOBVR = 'N'; LDVL = 1; LDVR = 1; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVL,WORK    ,LWORK,info)                                                      !:.:.:
   if(info /= 0) write(err,*)                      'array2_dgeev: dgeev failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = cmplx(ReEV,ImEV,kind=dp)
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   if(allocated(eigenvec)) deallocate(eigenvec); allocate(eigenvec(n,n))
   JOBVL = 'N' ; JOBVR = 'V'; LDVL = 1; LDVR = n; LWORK = -1
   allocate(VL(LDVL,n)); allocate(VR(LDVR,n))
   !Query optimal LWORK:
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORKTEST,LWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call dgeev(JOBVL,JOBVR,n,A,n,ReEV,ImEV,VL,LDVL,VR,LDVR,WORK    ,LWORK,info)
   if(info /= 0) write(err,*)                      'array2_dgeev: dgeev failed (V). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = cmplx(ReEV,ImEV,kind=dp)
   eigenvec = VR
   do concurrent (j=1:n, ImEV(j) > 0.0_dp)               ! If eigenvalues are complex conjugate pairs, then eigenvectors must be computed by assuming that the positive imaginary part comes first, then comes
    eigenvec(:,j  ) = cmplx(VR(:,j), VR(:,j+1),kind=dp)  ! the negative imaginary part (hence the mask ImEV(j) > 0.0_dp) and then applying the rule in the documentation of dgeev:
    eigenvec(:,j+1) = cmplx(VR(:,j),-VR(:,j+1),kind=dp)  !"If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR. If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, then
   end do                                                ! v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1)"
  !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call                               matrix_error('array2_dgeev: wrong value of job: '//trim(job))
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(WORK)
  deallocate(VL)    ; deallocate(VR)   ; deallocate(ReEV)  ; deallocate(ImEV)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_dgeev
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_dsyev                                           (C,eigenval,eigenvec,job)                                !:.:.:
  real     (dp), dimension(:,:)             , intent(in)                   :: C
  complex  (dp), dimension(:)               , intent(in out), allocatable  :: eigenval                                              ! if lbound is not one, the result was not correct. Make local allocation
  complex  (dp), dimension(:,:)             , intent(in out), allocatable  :: eigenvec
  character(*)                              , intent(in)                   :: job
  real     (dp), dimension(size(C,1),size(C,2))                            :: A
  integer                                                                  :: n
  integer                                                                  :: info,  LWORK
  real     (dp)                                                            :: WORKTEST(1)
  real     (dp), dimension(:)                               , allocatable  :: EV
  real     (dp), dimension(:)                               , allocatable  :: WORK
  character(1024)                                                          :: err
  character(1)                                                             :: JOBZ, UPLO

  n = size(C,1)

  if( size(C,1) /= size(C,2)  )  call matrix_error('array2_dsyev: Input matrix must be a square matrix:  m /= n')

  if(allocated(eigenval)      )  deallocate(eigenval) ; allocate(eigenval(n))

  allocate(EV(n))

  A = C  ! A local copy of the matrix C

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(job(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case('N','n')
   JOBZ='N' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,info)                                                                               !:.:.:
   if(info /= 0) write(err,*)                      'array2_dsyev: dsyev failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = EV
  !---------------------------------------------------------------------------------------------------------------------------------
  case('V','v')
   if(allocated(eigenvec)) deallocate(eigenvec); allocate(eigenvec(n,n))
   JOBZ='V' ; UPLO='U'; LWORK = -1
   !Query optimal LWORK:
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORKTEST,LWORK,info)
   LWORK = WORKTEST(1) + 10
   allocate(WORK(LWORK))
   call dsyev(JOBZ,UPLO,n,A,n,EV,WORK    ,LWORK,info)
   if(info /= 0) write(err,*)                      'array2_dsyev: dsyev failed (N). INFO= ',info
   if(info /= 0) call matrix_error(err)
   eigenval = EV
   eigenvec = A
  !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call                               matrix_error('array2_dsyev: wrong value of job: '//trim(job))
  !---------------------------------------------------------------------------------------------------------------------------------
  end select
  !---------------------------------------------------------------------------------------------------------------------------------
  deallocate(WORK)
  !---------------------------------------------------------------------------------------------------------------------------------
 end subroutine     array2_dsyev
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_determinant                                     (C)  result (det)                                        !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp)                                     :: det
  complex(dp), dimension(size(C,1),size(C,2))     :: A
  integer                                         :: n
  integer                                         :: info, isign, i
  integer    , dimension(size(C,1))               :: ipiv
  character(1000)                                 :: err

  if(size(C,1)/=size(C,2)) call matrix_error('array2_determinant: Matrix is not a square matrix')

  n = size(C,1)

  A = C

  call zgetrf(n,n,A,n,ipiv,info)                                                                                                    !:.:.:

  if( info /= 0          ) write(err,*)      'array2_determinant: zgetrf failed. INFO= ',info
  if( info /= 0          ) call matrix_error(err)

  isign = 1; det   = ONE

  do concurrent(i  = 1:n)
   det  = det * A(i,i)
  end do

  do concurrent(i  = 1:n)
   if(ipiv(i) /= i) isign = -isign
  end do
  
  det = det * isign

 end function       array2_determinant
!-----------------------------------------------------------------------------------------------------------------------------------
! For n x n matrices with gaussian random numbers, n >~ 192 requires the use of lnDet
!
! det(1) = log(|det(C)|)  det(2) = exp( i arg( det(C) ) )
! Det(C) = exp(det(1))  * det(2)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_determinant                                 (C)                      result(det)                     !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp), dimension(2)                       :: det
  complex(dp), dimension(size(C,1),size(C,2))     :: A
  integer                                         :: n
  integer                                         :: info, i
  integer    , dimension(size(C,1))               :: ipiv
  complex(dp)                                     :: phase
  real   (dp)                                     :: aAbs,lnDetAbs
  character(1000)                                 :: err

  if(size(C,1)/=size(C,2)) call matrix_error('array2_log_determinant: Matrix is not a square matrix')

  n = size(C,1)

  A = C

  call zgetrf(n,n,A,n,ipiv,info)                                                                                                    !:.:.:

  if( info /= 0          ) write(err,*)      'array2_log_determinant: zgetrf failed. INFO= ',info
  if( info /= 0          ) call matrix_error(err)

  lnDetAbs  = 0.0_dp
  phase     = ONE

  do i      = 1, n

   aAbs     = ABS(A(i,i))

   if(aAbs <= 0.0_dp) then
    det(1)  = ZERO
    det(2)  = ONE
    return
   end if

   lnDetAbs = lnDetAbs + LOG(aAbs)
   phase    = phase    * (A(i,i)/aAbs)

   if(ipiv(i) /= i )  phase = -phase

  end do

  det(1)    = lnDetAbs
  det(2)    = phase

 end function       array2_log_determinant
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_determinant_d                                   (C)                      result(det)                     !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C
  real   (dp)                                     :: det
  real   (dp), dimension(size(C,1),size(C,2))     :: A
  integer                                         :: n
  integer                                         :: info, isign, i
  integer    , dimension(size(C,1))               :: ipiv
  character(1000)                                 :: err

  if(size(C,1)/=size(C,2)) call matrix_error('array2_determinant_d: Matrix is not a square matrix')

  n = size(C,1)

  A = C

  call dgetrf(n,n,A,n,ipiv,info)                                                                                                    !:.:.:

  if( info /= 0          ) write(err,*)      'array2_determinant_d: dgetrf failed. INFO= ',info
  if( info /= 0          ) call matrix_error(err)

  isign = 1; det   = ONE

  do concurrent(i  = 1:n)
   det  = det * A(i,i)
  end do

  do concurrent(i  = 1:n)
   if(ipiv(i) /= i) isign = -isign
  end do
  
  det = det * isign

 end function       array2_determinant_d
!-----------------------------------------------------------------------------------------------------------------------------------
! For n x n matrices with gaussian random numbers, n >~ 192 requires the use of lnDet
!
! det(1) = log(|det(C)|)  det(2) = exp( i arg( det(C) ) )
! Det(C) = exp(det(1))  * det(2)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_determinant_d                               (C)                      result(det)                     !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C
  real   (dp), dimension(2)                       :: det
  real   (dp), dimension(size(C,1),size(C,2))     :: A
  integer                                         :: n
  integer                                         :: info, i, isign
  integer    , dimension(size(C,1))               :: ipiv
  real   (dp)                                     :: aAbs,lnDetAbs
  character(1000)                                 :: err

  if(size(C,1)/=size(C,2)) call matrix_error('array2_log_determinant_d: Matrix is not a square matrix')

  n = size(C,1)

  A = C

  call dgetrf(n,n,A,n,ipiv,info)                                                                                                    !:.:.:

  if( info /= 0          ) write(err,*)      'array2_log_determinant_d: dgetrf failed. INFO= ',info
  if( info /= 0          ) call matrix_error(err)

  lnDetAbs  = 0.0_dp
  isign     = 1

  do i      = 1, n

   aAbs     = A(i,i)

   if(aAbs <=  0 )then
    aAbs    = -aAbs
    isign   = -isign
   end if

   if(aAbs == 0.0_dp) then
    det(1)  = 0.0_dp
    det(2)  = 1.0_dp
    return
   end if

   lnDetAbs = lnDetAbs + LOG(aAbs)

   if(ipiv(i) /= i )  isign = -isign

  end do

  det(1)    = lnDetAbs
  det(2)    = isign

 end function       array2_log_determinant_d
!-----------------------------------------------------------------------------------------------------------------------------------
! After C code by Simon Catterall,  Phys. Rev. D 68, 014503 (2003) 
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
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp)                                     :: pfaffian
  complex(dp), dimension(size(C,1),size(C,2))     :: M
  integer                                         :: n
  integer                                         :: i,j,k,jpiv,interchange
  real   (dp)                                     :: pivot
  complex(dp),dimension(size(C,1))                :: dum,sc
  complex(dp)                                     :: scale,f

  
  if(size(C,1)/=size(C,2)) call matrix_error('array2_pfaffian: Matrix is not a square matrix')

  n = size(C,1)

  if(MOD(n,2) /= 0       ) call matrix_error('array2_pfaffian: Matrix is not even dimensional')

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
  
  f=DCMPLX(1.0D0,0.0D0)
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
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp), dimension(2)                       :: pfaffian
  complex(dp), dimension(size(C,1),size(C,2))     :: M
  integer                                         :: n
  integer                                         :: i,j,k,jpiv,interchange
  real   (dp)                                     :: pivot
  complex(dp),dimension(size(C,1))                :: dum,sc
  complex(dp)                                     :: scale,f
  real   (dp)                                     :: lnPfaffianAbs,AbsM
  complex(dp)                                     :: phase

  
  if(size(C,1)/=size(C,2)) call matrix_error('array2_pfaffian: Matrix is not a square matrix')

  n = size(C,1)

  if(MOD(n,2) /= 0       ) call matrix_error('array2_pfaffian: Matrix is not even dimensional')

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
    pfaffian(1) = ZERO
    pfaffian(2) = ONE
    return
   end if

   lnPfaffianAbs = lnPfaffianAbs + LOG(AbsM)
   phase         = phase         * (M(i+1,i)/AbsM)

  end do

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(n/2) we have to change the sign appropriately:
  if(MOD(n/2,2).EQ.1) interchange=-interchange

  phase       = phase * interchange

  pfaffian(1) = lnPfaffianAbs
  pfaffian(2) = phase

 end function array2_log_pfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
!Slower but easier coding. Keep it for reference. array2_log_pfaffian is approximately 3 times faster than array2_log_pfaffian2 (intel compiler)
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_pfaffian2                                       (C)                      result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp)                                     :: pfaffian
  complex(dp), dimension(size(C,1),size(C,2))     :: M
  integer                                         :: n
  integer                                         :: i,j,k,jpiv,interchange
  real   (dp)                                     :: pivot
  complex(dp),dimension(size(C,1))                :: dum
  complex(dp)                                     :: scale,f

  
  if(size(C,1)/=size(C,2)) call matrix_error('array2_pfaffian: Matrix is not a square matrix')

  n = size(C,1)

  if(MOD(n,2) /= 0       ) call matrix_error('array2_pfaffian: Matrix is not even dimensional')

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

  f=CMPLX(1.0D0,0.0D0)
  do i=1,n,2
   f=f*M(i+1,i)
  enddo

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(n/2) we have to change the sign appropriately:
  if(MOD(n/2,2).EQ.1) interchange=-interchange

  pfaffian=f*interchange

 end function       array2_pfaffian2
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array2_log_pfaffian2                                   (C)                      result(pfaffian)                !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp), dimension(2)                       :: pfaffian
  complex(dp), dimension(size(C,1),size(C,2))     :: M
  integer                                         :: n
  integer                                         :: i,j,k,jpiv,interchange
  real   (dp)                                     :: pivot
  complex(dp),dimension(size(C,1))                :: dum
  complex(dp)                                     :: scale,f
  real   (dp)                                     :: lnPfaffianAbs,AbsM
  complex(dp)                                     :: phase

  
  if(size(C,1)/=size(C,2)) call matrix_error('array2_pfaffian: Matrix is not a square matrix')

  n = size(C,1)

  if(MOD(n,2) /= 0       ) call matrix_error('array2_pfaffian: Matrix is not even dimensional')

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
    pfaffian(1) = ZERO
    pfaffian(2) = ONE
    return
   end if

   lnPfaffianAbs = lnPfaffianAbs + LOG(AbsM)
   phase         = phase         * (M(i+1,i)/AbsM)

  end do

  !since we interchanged rows with columns and Pf(A^T)=(-1)^(n/2) we have to change the sign appropriately:
  if(MOD(n/2,2).EQ.1) interchange=-interchange

  phase       = phase * interchange

  pfaffian(1) = lnPfaffianAbs
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
