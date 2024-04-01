!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_array.f90 .............................................................................!:.:.:
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
module              matrix_mod_array                                                                                                !:.:.:
 use matrix_mod_common
 implicit none
 save
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          operator(.mm.)                                                                                                  !:.:.:
  module procedure :: array2_matmul_array2, array2_matmul_array2_d, array2_d_matmul_array2, array2_d_matmul_array2_d
  module procedure :: array2_matmul_array1     , array1_matmul_array2
  module procedure :: array2_matmul_array1_d   , array1_d_matmul_array2
  module procedure :: array2_d_matmul_array1   , array1_matmul_array2_d
  module procedure :: array2_d_matmul_array1_d , array1_d_matmul_array2_d
 end interface      operator(.mm.)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          random_number                                                                                                   !:.:.:
  module procedure :: random_number_real_scalar_gaussian, random_number_complex_scalar_gaussian, random_number_complex_scalar
  module procedure :: random_number_array1,               random_number_array1_gaussian,         random_number_array1_gaussian_d
  module procedure :: random_number_array2,               random_number_array2_gaussian,         random_number_array2_gaussian_d
  module procedure :: random_number_array3,               random_number_array3_gaussian,         random_number_array3_gaussian_d
 end interface      random_number
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          mmmult                                                                                                          !:.:.:
  module procedure :: array2_matmul_lapack_mm_sub,array2_matmul_lapack_dd_sub,array2_matmul_lapack_m_sub,array2_matmul_lapack_d_sub
 end interface      mmmult
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          mvmult                                                                                                          !:.:.:
  module procedure :: array2_matmul_lapack_mv_sub,array2_matmul_lapack_dv_sub
 end interface      mvmult
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          vmmult                                                                                                          !:.:.:
  module procedure :: array2_matmul_lapack_vm_sub,array2_matmul_lapack_vd_sub
 end interface      vmmult
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          lmatmul                                                                                                         !:.:.:
  module procedure :: array2_matmul_lapack_mm,array2_matmul_lapack_dd, array2_matmul_lapack_m, array2_matmul_lapack_d
  module procedure :: array2_matmul_lapack_mv,array2_matmul_lapack_dv
 end interface      lmatmul
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          eigenvalues                                                                                                     !:.:.:
  module procedure :: array2_eigenvalues,     array2_eigenvalues_d
 end interface      eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          eigenvectors                                                                                                    !:.:.:
  module procedure :: array2_eigenvectors,    array2_eigenvectors_d
 end interface      eigenvectors
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          determinant                                                                                                     !:.:.:
  module procedure :: array2_determinant,     array2_determinant_d
 end interface      determinant
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          lndet                                                                                                           !:.:.:
  module procedure :: array2_log_determinant, array2_log_determinant_d
 end interface      lndet
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          pfaffian                                                                                                        !:.:.:
  module procedure :: array2_pfaffian
 end interface      pfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          lnPfaffian                                                                                                      !:.:.:
  module procedure :: array2_log_pfaffian
 end interface      lnPfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          inverse                                                                                                         !:.:.:
  module procedure :: array2_inverse,         array2_inverse_d
 end interface      inverse
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          norm                                                                                                            !:.:.:
  module procedure :: array3_norm, array3_norm_d
  module procedure :: array2_norm, array2_norm_d
  module procedure :: array1_norm, array1_norm_d
 end interface      norm
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          hermitian                                                                                                       !:.:.:
  module procedure :: array2_hermitian_get
 end interface      hermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          hermitian_set                                                                                                   !:.:.:
  module procedure :: array2_hermitian_set
 end interface      hermitian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          symmetric                                                                                                       !:.:.:
  module procedure :: array2_symmetric_get_d,array2_symmetric_get
 end interface      symmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          symmetric_set                                                                                                   !:.:.:
  module procedure :: array2_symmetric_set_d,array2_symmetric_set
 end interface      symmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          antisymmetric                                                                                                   !:.:.:
  module procedure :: array2_antisymmetric_get_d,array2_antisymmetric_get
 end interface      antisymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          antisymmetric_set                                                                                               !:.:.:
  module procedure :: array2_antisymmetric_set_d,array2_antisymmetric_set
 end interface      antisymmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          diagonal                                                                                                        !:.:.:
  module procedure :: array2_diagonal_get  ,array2_diagonal_set  ,array2_diagonal_set_from_complex ! get means obtaining the diagonal as vector, set means to contruct a diagonal matrix
  module procedure :: array2_diagonal_get_d,array2_diagonal_set_d,array2_diagonal_set_from_real_d
 end interface      diagonal
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          diagonalMatrix                                                                                                  !:.:.:
  module procedure :: array2_diagonal_set_from_complex_matrix    ,array2_diagonal_set_from_real_matrix_d
 end interface      diagonalMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          identitymatrix                                                                                                  !:.:.:
  module procedure :: array2_diagonal_set_identity_complex_matrix 
 end interface      identitymatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          cidentitymatrix                                                                                                  !:.:.:
  module procedure :: array2_diagonal_set_identity_complex_matrix 
 end interface      cidentitymatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          didentitymatrix                                                                                                  !:.:.:
  module procedure :: array2_diagonal_set_identity_real_matrix 
 end interface      didentitymatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          paulimatrix                                                                                                     !:.:.:
  module procedure :: array2_PauliMatrix
 end interface paulimatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          trace                                                                                                           !:.:.:
  module procedure :: array2_trace,array2_trace_d
 end interface      trace
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          trace2                                                                                                          !:.:.:
  module procedure :: array2_trace2,array2_trace2_d,array2_trace2_2matrices,array2_trace2_d_2matrices
 end interface      trace2
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          trace2c                                                                                                         !:.:.:
  module procedure :: array2_trace2_connected,array2_trace2_connected_d
 end interface      trace2c
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          traceless                                                                                                       !:.:.:
  module procedure :: array2_traceless_get, array2_traceless_get_d
 end interface      traceless
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          traceless_set                                                                                                   !:.:.:
  module procedure :: array2_traceless_set, array2_traceless_set_d
 end interface      traceless_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isHermitian                                                                                                     !:.:.:
  module procedure :: array2_is_Hermitian
 end interface      isHermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isSymmetric                                                                                                     !:.:.:
  module procedure :: array2_is_Symmetric,array2_is_Symmetric_d
 end interface      isSymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isAntiSymmetric                                                                                                 !:.:.:
  module procedure :: array2_is_AntiSymmetric,array2_is_AntiSymmetric_d
 end interface      isAntiSymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          sort                                                                                                            !:.:.:
  module procedure :: array1_sort,array1_sort_d
 end interface      sort
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          random_sort                                                                                                     !:.:.:
  module procedure :: random_sort_fun_i, random_sort_fun_d, random_sort_fun_c, random_sort_fun_char
  module procedure :: random_list_n    
 end interface      random_sort
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          random_sort_array                                                                                               !:.:.:
  module procedure :: random_sort_i    , random_sort_d    , random_sort_c
 end interface      random_sort_array
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          print                                                                                                           !:.:.:
  module procedure :: array2_print,array2_print_d
  module procedure :: array1_print,array1_print_d
 end interface      print
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          printna                                                                                                         !:.:.:
  module procedure :: array2_print_nonallocatable,array2_print_nonallocatable_d
  module procedure :: array1_print_nonallocatable,array1_print_nonallocatable_d
 end interface      printna
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          save                                                                                                            !:.:.:
  module procedure :: array3_save_matrix, array3_save_matrix_d
  module procedure :: array2_save_matrix, array2_save_matrix_d
  module procedure :: array1_save_matrix, array1_save_matrix_d
 end interface      save
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          read                                                                                                            !:.:.:
  module procedure :: array1_read_matrix, array1_read_matrix_d
  module procedure :: array2_read_matrix, array2_read_matrix_d
  module procedure :: array3_read_matrix, array3_read_matrix_d
 end interface      read
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isNaN                                                                                                           !:.:.:
  module procedure :: array3_is_nan, array3_is_nan_d
  module procedure :: array2_is_nan, array2_is_nan_d
  module procedure :: array1_is_nan, array1_is_nan_d
 end interface isNaN
!-----------------------------------------------------------------------------------------------------------------------------------
! Lapack and BLAS interfaces:  Notice that complex(8), real(8) are fixed to Lapack/BLAS requirements. Don't use parameter like dp.
!-----------------------------------------------------------------------------------------------------------------------------------
 interface
!-----------------------------------------------------------------------------------------------------------------------------------
! BLAS
  pure subroutine   zgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
   complex  (8), intent(in)             :: alpha, beta
   character(1), intent(in)             :: transa,transb
   integer     , intent(in)             :: m, n, k, lda, ldb, ldc
   complex  (8), intent(in)             :: A(lda,*), B(ldb,*)         ! Careful here, don't use C(:,:), use exactly as in LAPACK
   complex  (8), intent(inout)          :: C(ldc,*)
  end  subroutine   zgemm
  pure subroutine   zhemm(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc)
   complex  (8), intent(in)             :: alpha, beta
   character(1), intent(in)             :: side,uplo
   integer     , intent(in)             :: m,n,lda,ldb,ldc
   complex  (8), intent(in)             :: A(lda,*), B(ldb,*)
   complex  (8), intent(inout)          :: C(ldc,*)
  end  subroutine   zhemm
  pure subroutine   zgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy)
   complex  (8), intent(in)             :: alpha,beta
   character(1), intent(in)             :: trans
   integer     , intent(in)             :: m,n,lda,incx,incy
   complex  (8), intent(in)             :: A(lda,*), X(*)
   complex  (8), intent(inout)          :: Y(*)
  end  subroutine   zgemv
  pure subroutine   zhemv(uplo,n,alpha,A,lda,X,incx,beta,Y,incy)
   complex  (8), intent(in)             :: alpha,beta
   character(1), intent(in)             :: uplo
   integer     , intent(in)             :: n,lda,incx,incy
   complex  (8), intent(in)             :: A(lda,*),X(*)
   complex  (8), intent(inout)          :: Y(*)
  end  subroutine   zhemv
  pure subroutine   zherk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc)
   complex  (8), intent(in)             :: alpha,beta
   character(1), intent(in)             :: uplo,trans
   integer     , intent(in)             :: n,k,lda,ldc
   complex  (8), intent(in)             :: A(lda,*)
   complex  (8), intent(inout)          :: C(ldc,*)
  end  subroutine   zherk
  pure subroutine   dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
   real     (8), intent(in)             :: alpha, beta
   character(1), intent(in)             :: transa,transb
   integer     , intent(in)             :: m, n, k, lda, ldb, ldc
   real     (8), intent(in)             :: A(lda,*), B(ldb,*)
   real     (8), intent(inout)          :: C(ldc,*)
  end  subroutine   dgemm
  pure subroutine   dsymm(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc)
   real     (8), intent(in)             :: alpha, beta
   character(1), intent(in)             :: side,uplo
   integer     , intent(in)             :: m,n,lda,ldb,ldc
   real     (8), intent(in)             :: A(lda,*), B(ldb,*)
   real     (8), intent(inout)          :: C(ldc,*)
  end  subroutine   dsymm
  pure subroutine   dsymv(uplo,n,alpha,A,lda,X,incx,beta,Y,incy)
   real     (8), intent(in)             :: alpha,beta
   character(1), intent(in)             :: uplo
   integer     , intent(in)             :: n,lda,incx,incy
   real     (8), intent(in)             :: A(lda,*),X(*)
   real     (8), intent(inout)          :: Y(*)
  end  subroutine   dsymv
  pure subroutine   dgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy)
   real     (8), intent(in)             :: alpha,beta
   character(1), intent(in)             :: trans
   integer     , intent(in)             :: m,n,lda,incx,incy
   real     (8), intent(in)             :: A(lda,*), X(*)
   real     (8), intent(inout)          :: Y(*)
  end  subroutine   dgemv
  pure subroutine   dsyrk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc)
   real     (8), intent(in)             :: alpha,beta
   character(1), intent(in)             :: uplo,trans
   integer     , intent(in)             :: n,k,lda,ldc
   real     (8), intent(in)             :: A(lda,*)
   real     (8), intent(inout)          :: C(ldc,*)
  end  subroutine   dsyrk
! end BLAS
!-----------------------------------------------------------------------------------------------------------------------------------
! Lapack
  pure subroutine   zgeev(jobvl,jobvr,n,A,lda,W,VL,ldvl,VR,ldvr,WORK,lwork,RWORK,info)
   character(1), intent(in)             :: jobvl,jobvr
   integer     , intent(in)             :: n,lda,ldvl,ldvr,lwork
   complex  (8), intent(inout)          :: A(lda,*),W(*),VL(ldvl,*),VR(ldvr,*),WORK(*)
   real     (8), intent(inout)          :: RWORK(*)
   integer     , intent(inout)          :: info
  end  subroutine   zgeev
  pure subroutine   zheev(jobz,uplo,n,A,lda,W,WORK,lwork,RWORK,info)
   character(1), intent(in)             :: jobz,uplo
   integer     , intent(in)             :: n,lda,lwork
   complex  (8), intent(inout)          :: A(lda,*),WORK(*)
   real     (8), intent(inout)          :: W(*),RWORK(*)
   integer     , intent(inout)          :: info
  end  subroutine   zheev
  pure subroutine   zgetrf(m,n,A,lda,IPIV,info)
   integer     , intent(in)             :: m,n,lda
   complex  (8), intent(inout)          :: A(lda,*)
   integer     , intent(inout)          :: IPIV(*)
   integer     , intent(inout)          :: info
  end  subroutine   zgetrf
  pure subroutine   zgetri(n,A,lda,IPIV,WORK,lwork,info)
   integer     , intent(in)             :: n,lda,lwork
   integer     , intent(in)             :: IPIV(*)
   complex  (8), intent(inout)          :: A(lda,*)
   complex  (8), intent(inout)          :: WORK(*)
   integer     , intent(inout)          :: info
  end  subroutine   zgetri
  pure subroutine   dgeev(jobvl,jobvr,n,A,lda,WR,WI,VL,ldvl,VR,ldvr,WORK,lwork,info)
   character(1), intent(in)             :: jobvl,jobvr
   integer     , intent(in)             :: n,lda,ldvl,ldvr,lwork
   real     (8), intent(inout)          :: A(lda,*),WR(*),WI(*),VL(ldvl,*),VR(ldvr,*),WORK(*)
   integer     , intent(inout)          :: info
  end  subroutine   dgeev
  pure subroutine   dsyev(jobz,uplo,n,A,lda,W,WORK,lwork,info)
   character(1), intent(in)             :: jobz,uplo
   integer     , intent(in)             :: n,lda,lwork
   real     (8), intent(inout)          :: A(lda,*),W(*),WORK(*)
   integer     , intent(inout)          :: info
  end  subroutine   dsyev
  pure subroutine   dgetrf(m,n,A,lda,IPIV,info)
   integer     , intent(in)             :: m,n,lda
   real     (8), intent(inout)          :: A(lda,*)
   integer     , intent(inout)          :: IPIV(*)
   integer     , intent(inout)          :: info
  end  subroutine   dgetrf
  pure subroutine   dgetri(n,A,lda,IPIV,WORK,lwork,info)
   integer     , intent(in)             :: n,lda,lwork
   integer     , intent(in)             :: IPIV(*)
   real     (8), intent(inout)          :: A(lda,*)
   real     (8), intent(inout)          :: WORK(*)
   integer     , intent(inout)          :: info
  end  subroutine   dgetri
! end Lapack
!-----------------------------------------------------------------------------------------------------------------------------------
 end interface
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------


!...................................................................................................................................
!.......................... array2 matrix functions ................................................................................!:.:.:
!...................................................................................................................................
 include 'matrix_mod_array_lapack.f90'

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace2_connected                                (C,mtype)                result(t)                       !:.:.:
  complex  (dp), dimension(:,:),intent(in)        :: C
  character(*) , optional      ,intent(in)        :: mtype
  complex  (dp)                                   :: t
  complex  (dp), dimension(size(C,1),size(C,2))   :: Cc
  character(mtype_len)                            :: type
  integer                                         :: i,j

  type = 'ZG'

  if(present(mtype))  type = mtype
  
  Cc = C; call array2_traceless_set(Cc)
  
  t  = array2_trace2(Cc,type)

  end function      array2_trace2_connected
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace2_connected_d                              (C,mtype)                result(t)                       !:.:.:
  real     (dp), dimension(:,:),intent(in)        :: C
  character(*) , optional      ,intent(in)        :: mtype
  real     (dp)                                   :: t
  real     (dp), dimension(size(C,1),size(C,2))   :: Cc
  character(mtype_len)                            :: type
  integer                                         :: i,j

  type = 'DG'

  if(present(mtype))  type = mtype
  
  Cc = C; call array2_traceless_set_d(Cc)
  
  t  = array2_trace2_d(Cc,type)

  end function      array2_trace2_connected_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace2                                          (C,mtype)                result(t)                       !:.:.:
  complex  (dp), dimension(:,:),intent(in)        :: C
  character(*) , optional      ,intent(in)        :: mtype
  complex  (dp)                                   :: t
  character(mtype_len)                            :: type
  integer                                         :: i,j

  type = 'G'

  t    = ZERO

  if(present(mtype))  type = mtype(2:2)

  select case (type)
   case('H') ! faster computation, if C is hermitian
    do concurrent( i = 1:size(C,1), j = 1:size(C,2) )
     t = t + C(i,j) * CONJG(C(i,j))
    end do
   case default
    do concurrent( i = 1:size(C,1), j = 1:size(C,2) )
     t = t + C(i,j) * C(j,i)
    end do
   end select

 end function       array2_trace2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace2_2matrices                                (C1,C2,mtype)            result(t)                       !:.:.:
  complex  (dp), dimension(:,:),intent(in)        :: C1,C2
  character(*) , optional      ,intent(in)        :: mtype
  complex  (dp)                                   :: t
  character(mtype_len)                            :: type
  integer                                         :: i,j,m1,m2,n1,n2

  type = 'G'

  t    = ZERO

  if(present(mtype))  type = mtype(2:2)

  select case (type)
   case('H') ! faster computation, if C2 is hermitian
    do concurrent( i = 1:size(C1,1), j = 1:size(C1,2) )
     t = t + C1(i,j) * CONJG(C2(i,j))
    end do
   case default
    do concurrent( i = 1:size(C1,1), j = 1:size(C1,2) )
     t = t + C1(i,j) * C2(j,i)
    end do
   end select

  m1 = size(C1,1) ; m2 = size(C2,1)
  n1 = size(C1,2) ; n2 = size(C2,2)
  if (n1 /= m1 ) t= nan() 
  if (n2 /= m2 ) t= nan()
  if (n1 /= n2 ) t= nan() 

 end function       array2_trace2_2matrices
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace2_d                                        (C,mtype)                result(t)                       !:.:.:
  real     (dp), dimension(:,:),intent(in)        :: C
  character(*) , optional      ,intent(in)        :: mtype
  real     (dp)                                   :: t
  character(mtype_len)                            :: type
  integer                                         :: i,j

  type = 'G'

  t    = 0.0_dp

  if(present(mtype))  type = mtype(2:2)

  select case (type)
   case('S') ! faster computation, if C is hermitian
    do concurrent( i = 1:size(C,1), j = 1:size(C,2) )
     t = t + C(i,j) * C(i,j)
    end do
   case default
    do concurrent( i = 1:size(C,1), j = 1:size(C,2) )
     t = t + C(i,j) * C(j,i)
    end do
   end select

  end function      array2_trace2_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace2_d_2matrices                             (C1,C2,mtype)             result(t)                       !:.:.:
  real     (dp), dimension(:,:),intent(in)        :: C1,C2
  character(*) , optional      ,intent(in)        :: mtype
  real     (dp)                                   :: t
  character(mtype_len)                            :: type
  integer                                         :: i,j,m1,m2,n1,n2

  type = 'G'

  t    = 0.0_dp

  if(present(mtype))  type = mtype(2:2)

  select case (type)
   case('S') ! faster computation, if C2 is hermitian
    do concurrent( i = 1:size(C1,1), j = 1:size(C1,2) )
     t = t + C1(i,j) * C2(i,j)
    end do
   case default
    do concurrent( i = 1:size(C1,1), j = 1:size(C1,2) )
     t = t + C1(i,j) * C2(j,i)
    end do
   end select

  m1 = size(C1,1) ; m2 = size(C2,1)
  n1 = size(C1,2) ; n2 = size(C2,2)
  if (n1 /= m1 ) t= nan() 
  if (n2 /= m2 ) t= nan()
  if (n1 /= n2 ) t= nan() 

 end function       array2_trace2_d_2matrices
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_traceless_get                                   (C)                      result(B)                       !:.:.:
  complex(dp), dimension(:,:),    intent(in)      :: C
  complex(dp), dimension(size(C,1),size(C,1))     :: B
  complex(dp)                                     :: t
  integer                                         :: i

  t = - array2_trace(C) / size(C,1)

  do concurrent ( i = 1:size(C,1) )
   B(i,i) = C(i,i) + t
  end do

 end function       array2_traceless_get
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_traceless_get_d                                 (C)                      result(B)                       !:.:.:
  real   (dp), dimension(:,:),    intent(in)      :: C
  real   (dp), dimension(size(C,1),size(C,1))     :: B
  real   (dp)                                     :: t
  integer                                         :: i

  t = - array2_trace_d(C) / size(C,1)

  do concurrent ( i = 1:size(C,1) )
   B(i,i) = C(i,i) + t
  end do

 end function       array2_traceless_get_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_traceless_set                                   (C)                                                      !:.:.:
  complex(dp), dimension(:,:),intent(in out)      :: C
  complex(dp)                                     :: t
  integer                                         :: i

  t = - array2_trace(C) / size(C,1)

  do concurrent ( i = 1:size(C,1) )
   C(i,i) = C(i,i) + t
  end do

 end subroutine     array2_traceless_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_traceless_set_d                                 (C)                                                      !:.:.:
  real   (dp), dimension(:,:),intent(in out)      :: C
  real   (dp)                                     :: t
  integer                                         :: i

  t = - array2_trace_d(C) / size(C,1)

  do concurrent ( i = 1:size(C,1) )
   C(i,i) = C(i,i) + t
  end do

 end subroutine     array2_traceless_set_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace                                           (C)                      result(t)                       !:.:.:
  complex(dp), dimension(:,:),intent(in)          :: C
  complex(dp)                                     :: t

  t = SUM(array2_diagonal_get(C))

 end function       array2_trace
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_trace_d                                         (C)                      result(t)                       !:.:.:
  real   (dp), dimension(:,:),intent(in)          :: C
  real   (dp)                                     :: t

  t = SUM(array2_diagonal_get_d(C))

 end function       array2_trace_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_get                                    (C)                      result(d)                       !:.:.:
  complex(dp), dimension(:,:),intent(in)          :: C
  complex(dp), dimension(size(C,1))               :: d
  integer                                         :: i

  do concurrent( i = 1:size(d) )
   d(i) = C(i,i)
  end do

 end function       array2_diagonal_get
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_get_d                                  (C)                      result(d)                       !:.:.:
  real   (dp), dimension(:,:),intent(in)          :: C
  real   (dp), dimension(size(C,1))               :: d
  integer                                         :: i

  do concurrent( i = 1:size(d) )
   d(i) = C(i,i)
  end do

 end function       array2_diagonal_get_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set                                    (d)                      result(C)                       !:.:.:
  complex(dp), dimension(:), intent(in)           :: d
  complex(dp), dimension(size(d),size(d))         :: C
  integer                                         :: i

  C = ZERO
  do concurrent( i = 1:size(d) )
   C(i,i) = d(i)
  end do

 end function       array2_diagonal_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_d                                  (d)                      result(C)                       !:.:.:
  real   (dp), dimension(:), intent(in)           :: d
  real   (dp), dimension(size(d),size(d))         :: C
  integer                                         :: i

  C = 0.0_dp
  do concurrent( i = 1:size(d) )
   C(i,i) = d(i)
  end do

 end function       array2_diagonal_set_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_from_real                          (r,n)                    result(C)                       !:.:.:
  real   (dp), intent(in)                         :: r
  integer    , intent(in)                         :: n
  complex(dp), dimension(n,n)                     :: C
  integer                                         :: i

  C = ZERO

  do concurrent( i = 1:n )
   C(i,i) = r
  end do

 end function       array2_diagonal_set_from_real
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_from_real_d                        (r,n)                    result(C)                       !:.:.:
  real   (dp), intent(in)                         :: r
  integer    , intent(in)                         :: n
  real   (dp), dimension(n,n)                     :: C
  integer                                         :: i

  C = 0.0_dp

  do concurrent( i = 1:n )
   C(i,i) = r
  end do

 end function       array2_diagonal_set_from_real_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_from_complex                       (z,n)                    result(C)                       !:.:.:
  complex(dp), intent(in)                         :: z
  integer    , intent(in)                         :: n
  complex(dp), dimension(n,n)                     :: C
  integer                                         :: i

  C = ZERO

  do concurrent( i = 1:n )
   C(i,i) = z
  end do

 end function       array2_diagonal_set_from_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_from_complex_matrix                (C)                      result(D)                       !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  complex(dp), dimension(size(C,1),size(C,1))     :: D
  integer                                         :: i

  D = ZERO

  do concurrent( i = 1:size(C,1) )
   D(i,i) = C(i,i)
  end do

 end function       array2_diagonal_set_from_complex_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_from_real_matrix_d                 (C)                      result(D)                       !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C
  real   (dp), dimension(size(C,1),size(C,1))     :: D
  integer                                         :: i

  D = 0.0_dp

  do concurrent( i = 1:size(C,1) )
   D(i,i) = C(i,i)
  end do

 end function       array2_diagonal_set_from_real_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_identity_complex_matrix            (n)                      result(C)                       !:.:.:
  integer    , intent   (in)                      :: n
  complex(dp), dimension(n,n)                     :: C
  integer                                         :: i

  C = ZERO

  do concurrent( i = 1:n )
   C(i,i) = 1.0_dp
  end do


 end function       array2_diagonal_set_identity_complex_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_diagonal_set_identity_real_matrix               (n)                      result(C)                       !:.:.:
  integer    , intent   (in)                      :: n
  real   (dp), dimension(n,n)                     :: C
  integer                                         :: i
  C = 0.0_dp

  do concurrent( i = 1:n )
   C(i,i) = 1.0_dp
  end do


 end function       array2_diagonal_set_identity_real_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_Hermitian                                    (C)                      result(r)                       !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  real   (dp)                                     :: r
  integer                                         :: i

  r = norm( C - CONJG(TRANSPOSE(C)))

 end function       array2_is_Hermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_Symmetric                                    (C)                      result(r)                       !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  real   (dp)                                     :: r

  r = norm( C - TRANSPOSE(C))

 end function       array2_is_Symmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_Symmetric_d                                  (C)                      result(r)                       !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C
  real   (dp)                                     :: r

  r = norm( C - TRANSPOSE(C))

 end function       array2_is_Symmetric_d

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_AntiSymmetric                                (C)                      result(r)                       !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C
  real   (dp)                                     :: r

  r = norm( C + TRANSPOSE(C))

 end function       array2_is_AntiSymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_AntiSymmetric_d                              (C)                      result(r)                       !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C
  real   (dp)                                     :: r

  r = norm( C + TRANSPOSE(C))

 end function       array2_is_AntiSymmetric_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_PauliMatrix                                     (n)                      result(C)                       !:.:.:
  integer    , intent(in)                         :: n
  complex(dp), dimension(2,2)                     :: C


  select case(n)
   case(0)
    C(1,1)  =  ONE   ; C(1,2) =  ZERO
    C(2,1)  =  ZERO  ; C(2,2) =  ONE
   case(1)
    C(1,1)  =  ZERO  ; C(1,2) =  ONE
    C(2,1)  =  ONE   ; C(2,2) =  ZERO
   case(2)
    C(1,1)  =  ZERO  ; C(1,2) = -IMU
    C(2,1)  =  IMU   ; C(2,2) =  ZERO
   case(3)
    C(1,1)  =  ONE   ; C(1,2) =  ZERO
    C(2,1)  =  ZERO  ; C(2,2) = -ONE
   case(4)
    C(1,1)  =  ONE   ; C(1,2) =  ZERO
    C(2,1)  =  ZERO  ; C(2,2) =  ONE
   case default
    C(:,:)  =  nan()
  end select

 end function       array2_PauliMatrix



!...................................................................................................................................
!.......................... array2 procedures.......................................................................................!:.:.:
!...................................................................................................................................


!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array3_norm                                            (C)                      result(r)                       !:.:.:
  complex(dp), dimension(:,:,:),intent(in)        :: C
  real   (dp)                                     :: r

  r = SUM(ABS(C)) / SIZE(C)

 end function       array3_norm
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array3_norm_d                                          (C)                      result(r)                       !:.:.:
  real   (dp), dimension(:,:,:),intent(in)        :: C
  real   (dp)                                     :: r

  r = SUM(ABS(C)) / SIZE(C)

 end function       array3_norm_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_norm                                            (C)                      result(r)                       !:.:.:
  complex(dp), dimension(:,:),intent(in)          :: C
  real   (dp)                                     :: r

  r = SUM(ABS(C)) / SIZE(C)

 end function       array2_norm
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_norm_d                                          (C)                      result(r)                       !:.:.:
  real   (dp), dimension(:,:),intent(in)          :: C
  real   (dp)                                     :: r

  r = SUM(ABS(C)) / SIZE(C)

 end function       array2_norm_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_norm                                            (C)                      result(r)                       !:.:.:
  complex(dp), dimension(:)  ,intent(in)          :: C
  real   (dp)                                     :: r

  r = SUM(ABS(C)) / SIZE(C)

 end function       array1_norm
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_norm_d                                          (C)                      result(r)                       !:.:.:
  real   (dp), dimension(:)  ,intent(in)          :: C
  real   (dp)                                     :: r

  r = SUM(ABS(C)) / SIZE(C)

 end function       array1_norm_d
!-----------------------------------------------------------------------------------------------------------------------------------
!The matrix is assumed a square matrix:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_hermitian_set                                   (C,uplo)                                                 !:.:.:
  complex(dp), dimension(:,:),intent(in out)      :: C
  character(*),optional      ,intent(in)          :: uplo
  integer                                         :: i,j
  character(1)                                    :: ul

  ul = 't'; if(present(uplo)) ul = uplo(1:1)

  select case (ul)
  case('u','U') ! conserve the upper part of the matrix i<j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i>j)
    C(i,j) = CONJG(C(j,i))
   end do
  case('l','L') ! conserve the lower part of the matrix i>j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i<j)
    C(i,j) = CONJG(C(j,i))
   end do
  case default
   C  = 0.5_dp * ( C + CONJG(TRANSPOSE(C)) )        ! careful: gaussian random number generation for Matrix/DMatrix depends on this relation: sigma -> sigma*sqrt(2) for off diagonal for X -> (X1+X2)/2
  end select

  do concurrent   (i=1:size(C,1))
   C(i,i) = real(C(i,i), kind=dp)
  end do

 end subroutine     array2_hermitian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_hermitian_get                                   (C)                      result(CH)                      !:.:.:
  complex(dp), dimension(:,:),intent(in)          :: C
  complex(dp), dimension(size(C,2),size(C,1))     :: CH
  integer                                         :: i,j


  CH = CONJG(TRANSPOSE(C))


 end function       array2_hermitian_get
!-----------------------------------------------------------------------------------------------------------------------------------
!The matrix is assumed a square matrix:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_symmetric_set                                   (C,uplo)                                                 !:.:.:
  complex(dp), dimension(:,:),intent(in out)      :: C
  character(*),optional      ,intent(in)          :: uplo
  integer                                         :: i,j
  character(1)                                    :: ul

  ul = 't'; if(present(uplo)) ul = uplo(1:1)

  select case (ul)
  case('u','U') ! conserve the upper part of the matrix i<j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i>j)
    C(i,j) =       C(j,i)
   end do
  case('l','L') ! conserve the lower part of the matrix i>j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i<j)
    C(i,j) =       C(j,i)
   end do
  case default
   C  = 0.5_dp * ( C +      TRANSPOSE(C) )        ! careful: gaussian random number generation for Matrix/DMatrix depends on this relation: sigma -> sigma*sqrt(2) for off diagonal for X -> (X1+X2)/2
  end select

 end subroutine     array2_symmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_symmetric_get                                   (C)                      result(CS)                      !:.:.:
  complex(dp), dimension(:,:),intent(in)          :: C
  complex(dp), dimension(size(C,2),size(C,1))     :: CS

  CS = TRANSPOSE(C) 

 end function       array2_symmetric_get
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_symmetric_set_d                                 (C,uplo)                                                 !:.:.:
  real   (dp), dimension(:,:),intent(in out)      :: C
  character(*),optional      ,intent(in)          :: uplo
  integer                                         :: i,j
  character(1)                                    :: ul

  ul = 't'; if(present(uplo)) ul = uplo(1:1)

  select case (ul)
  case('u','U') ! conserve the upper part of the matrix i<j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i>j)
    C(i,j) =       C(j,i)
   end do
  case('l','L') ! conserve the lower part of the matrix i>j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i<j)
    C(i,j) =       C(j,i)
   end do
  case default
   C  = 0.5_dp * ( C +      TRANSPOSE(C) )        ! careful: gaussian random number generation for Matrix/DMatrix depends on this relation: sigma -> sigma*sqrt(2) for off diagonal for X -> (X1+X2)/2
  end select

 end subroutine     array2_symmetric_set_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_symmetric_get_d                                 (C)                      result(CS)                      !:.:.:
  real   (dp), dimension(:,:),intent(in)          :: C
  real   (dp), dimension(size(C,2),size(C,1))     :: CS

  CS = TRANSPOSE(C) 

 end function       array2_symmetric_get_d

!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_antisymmetric_set                               (C,uplo)                                                 !:.:.:
  complex(dp), dimension(:,:),intent(in out)      :: C
  character(*),optional      ,intent(in)          :: uplo
  integer                                         :: i,j
  character(1)                                    :: ul

  ul = 't'; if(present(uplo)) ul = uplo(1:1)

  select case (ul)
  case('u','U') ! conserve the upper part of the matrix i<j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i>j)
    C(i,j) =      -C(j,i)
   end do
  case('l','L') ! conserve the lower part of the matrix i>j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i<j)
    C(i,j) =      -C(j,i)
   end do
  case default
   C  = 0.5_dp * ( C -      TRANSPOSE(C) )        ! careful: gaussian random number generation for Matrix/DMatrix depends on this relation: sigma -> sigma*sqrt(2) for off diagonal for X -> (X1+X2)/2
  end select

  do concurrent   (i=1:size(C,1))
   C(i,i) = ZERO
  end do

 end subroutine     array2_antisymmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_antisymmetric_get                               (C)                      result(CS)                      !:.:.:
  complex(dp), dimension(:,:),intent(in)          :: C
  complex(dp), dimension(size(C,2),size(C,1))     :: CS

  CS =-TRANSPOSE(C) 

 end function       array2_antisymmetric_get
!-----------------------------------------------------------------------------------------------------------------------------------
 pure subroutine    array2_antisymmetric_set_d                             (C,uplo)                                                 !:.:.:
  real   (dp), dimension(:,:),intent(in out)      :: C
  character(*),optional      ,intent(in)          :: uplo
  integer                                         :: i,j
  character(1)                                    :: ul

  ul = 't'; if(present(uplo)) ul = uplo(1:1)

  select case (ul)
  case('u','U') ! conserve the upper part of the matrix i<j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i>j)
    C(i,j) =      -C(j,i)
   end do
  case('l','L') ! conserve the lower part of the matrix i>j
   do concurrent   (j = 1:size(C,2), i = 1:size(C,1), i<j)
    C(i,j) =      -C(j,i)
   end do
  case default
   C  = 0.5_dp * ( C -      TRANSPOSE(C) )        ! careful: gaussian random number generation for Matrix/DMatrix depends on this relation: sigma -> sigma*sqrt(2) for off diagonal for X -> (X1+X2)/2
  end select

  do concurrent   (i=1:size(C,1))
   C(i,i) = 0.0_dp
  end do

 end subroutine     array2_antisymmetric_set_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_antisymmetric_get_d                             (C)                      result(CS)                      !:.:.:
  real   (dp), dimension(:,:),intent(in)          :: C
  real   (dp), dimension(size(C,2),size(C,1))     :: CS

  CS =-TRANSPOSE(C) 

 end function       array2_antisymmetric_get_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_gauss_set                                       (C,sigma)                                                !:.:.:
  complex(dp), dimension(:,:),intent(out)         :: C
  real   (dp), optional      ,intent(in)          :: sigma
  real   (dp), dimension(size(C,1),size(C,2))     :: x,y
  real   (dp)                                     :: s
  
  s = 1.0_dp; if(present(sigma)) s = sigma
  call random_number(x); call random_number(y)
  C = s * sqrt(-2.0_dp * log(x)) * CMPLX(cos(TWOPI*y),sin(TWOPI*y),kind=dp)
  
 end subroutine     array2_gauss_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_gauss_set_d                                     (C,sigma)                                                !:.:.:
  real   (dp), dimension(:,:),intent(out)         :: C
  real   (dp), optional      ,intent(in)          :: sigma
  real   (dp), dimension(size(C,1),size(C,2))     :: x,y
  real   (dp)                                     :: s

  s = 1.0_dp; if(present(sigma)) s = sigma
  call random_number(x); call random_number(y)
  C = s * sqrt(-2.0_dp * log(x)) * cos(TWOPI*y)

 end subroutine     array2_gauss_set_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_random_set                                      (C)                                                      !:.:.:
  complex(dp), dimension(:,:)                     :: C
  real   (dp), dimension(size(C,1),size(C,2))     :: x,y
  
  call random_number(x); call random_number(y)
  C = CMPLX(x,y,kind=dp)

 end subroutine     array2_random_set
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_random_set_d                                    (C)                                                      !:.:.:
  real   (dp), dimension(:,:)                     :: C

  call random_number(C); 

 end subroutine     array2_random_set_d
!-----------------------------------------------------------------------------------------------------------------------------------
! Sort using options: reverse: from largest to smallest   forward: from smallest to largest
! MR/MF   modulus   reverse/modulus   forward
! RR/RF   real part reverse/real part forward
! IR/IF   imag part reverse/imag part forward
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array1_sort                                            (C,P,by)                 result(D)                       !:.:.:
  complex(dp), dimension(:), intent(in)           :: C
  integer    , dimension(:), optional             :: P
  character(*), optional   , intent(in)           :: by
  complex(dp), dimension(size(C,1))               :: D
  character(2)                                    :: b
  integer                                         :: i,j
  complex(dp)                                     :: c_i

  b = 'MF';                                               ! default value 
  if(present(by)) b = by(1:2)
  D = C

  if(present(P))then

   p = [(i, i=1,size(C,1))]

   select case (b)
   case('mr','MR')
    call array1_reversequicksortZbyModulus_pos (D,P,1,size(D,1)) ! reverse sorting by modulus:   from largest  to smallest
   case('mf','MF')
    call array1_quicksortZbyModulus_pos        (D,P,1,size(D,1)) ! forward sorting by modulus:   from smallest to largest
   case('rr','RR')
    call array1_reversequicksortZbyRealPart_pos(D,P,1,size(D,1)) ! reverse sorting by real part: from largest  to smallest
   case('rf','RF')
    call array1_quicksortZbyRealPart_pos       (D,P,1,size(D,1)) ! reverse sorting by real part: from smallest to largest
   case('ir','IR')
    call array1_reversequicksortZbyImagPart_pos(D,P,1,size(D,1)) ! reverse sorting by imag part: from largest  to smallest
   case('if','IF')
    call array1_quicksortZbyImagPart_pos       (D,P,1,size(D,1)) ! reverse sorting by imag part: from smallest to largest
   case default
    call matrix_error('aray1_sort: wrong "by" option: '//b)
   end select

  else

   select case (b)
   case('mr','MR')
    call array1_reversequicksortZbyModulus     (D,  1,size(D,1)) ! reverse sorting by modulus:   from largest  to smallest
   case('mf','MF')
    call array1_quicksortZbyModulus            (D,  1,size(D,1)) ! forward sorting by modulus:   from smallest to largest
   case('rr','RR')
    call array1_reversequicksortZbyRealPart    (D,  1,size(D,1)) ! reverse sorting by real part: from largest  to smallest
   case('rf','RF')
    call array1_quicksortZbyRealPart           (D,  1,size(D,1)) ! reverse sorting by real part: from smallest to largest
   case('ir','IR')
    call array1_reversequicksortZbyImagPart    (D,  1,size(D,1)) ! reverse sorting by imag part: from largest  to smallest
   case('if','IF')
    call array1_quicksortZbyImagPart           (D,  1,size(D,1)) ! reverse sorting by imag part: from smallest to largest
   case default
    call matrix_error('aray1_sort: wrong "by" option: '//b)
   end select



  end if
  !DEBUG: careful, if degeneracy is present, then not all indices are included in P(i), some are repeated
  !if(present(P))then !we record the position of each element in D in the original array C
  ! P            = -huge(i)
  ! iloop:  do i = 1, size(C,1)
  !  c_i         = D(i)
  !  jloop: do j = 1, size(C,1)
  !   if( c_i   == C(j)) then 
  !    P   (i)   =   j
  !    exit  jloop
  !   end if
  !  end  do jloop
  ! end   do iloop
  ! !test if error:
  ! do i         = 1,  size(C,1)
  !  if( P(i)   == -huge(i)) call matrix_error('sort: array1_sort: failed to ')
  ! end do
  !end    if ! if(present(P))

 end function       array1_sort
!-----------------------------------------------------------------------------------------------------------------------------------
! Sort using options: reverse: from largest to smallest   forward: from smallest to largest
! RR/RF   real     value reverse/real part forward
! AR/AF   absolute value reverse/imag part forward
!-----------------------------------------------------------------------------------------------------------------------------------
 function           array1_sort_d                                          (C,P,by)                 result(D)                       !:.:.:
  real   (dp), dimension(:), intent(in)           :: C
  integer    , dimension(:), optional             :: P
  character(*), optional   , intent(in)           :: by
  real   (dp), dimension(size(C,1))               :: D
  character(2)                                    :: b
  integer                                         :: i,j
  real   (dp)                                     :: c_i

  b = 'RF';                                               ! default value 
  if(present(by)) b = by(1:2)
  D = C

  if(present(P))then

   p = [(i, i=1,size(C,1))]

   select case (b)
   case('vr','VR','rr','RR')
    call array1_reversequicksortDbyValue_d_pos  (D,P,1,size(D,1)) ! reverse sorting by real     value: from largest  to smallest
   case('vf','VF','rf','RF')
    call array1_quicksortDbyValue_d_pos         (D,P,1,size(D,1)) ! forward sorting by real     value: from smallest to largest
   case('mr','mR','ar','AR')
    call array1_reversequicksortDbyModulus_d_pos(D,P,1,size(D,1)) ! reverse sorting by absolute value: from largest  to smallest
   case('mf','mF','af','AF')
    call array1_quicksortDbyModulus_d_pos       (D,P,1,size(D,1)) ! reverse sorting by absolute value: from smallest to largest
   case default
    call matrix_error('aray1_sort_d: wrong "by" option: '//b)
   end select

  else

   select case (b)
   case('vr','VR','rr','RR')
    call array1_reversequicksortDbyValue_d      (D,  1,size(D,1)) ! reverse sorting by real     value: from largest  to smallest
   case('vf','VF','rf','RF')
    call array1_quicksortDbyValue_d             (D,  1,size(D,1)) ! forward sorting by real     value: from smallest to largest
   case('mr','mR','ar','AR')
    call array1_reversequicksortDbyModulus_d    (D,  1,size(D,1)) ! reverse sorting by absolute value: from largest  to smallest
   case('mf','mF','af','AF')
    call array1_quicksortDbyModulus_d           (D,  1,size(D,1)) ! reverse sorting by absolute value: from smallest to largest
   case default
    call matrix_error('aray1_sort_d: wrong "by" option: '//b)
   end select

  end if

 end function       array1_sort_d
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortZbyModulus                           (A,first,last)                                           !:.:.:
  complex(dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) < ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) < ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortZbyModulus(a, first, i-1 )
  if (j+1   < last)  call array1_quicksortZbyModulus(a, j+1  , last)

 end subroutine     array1_quicksortZbyModulus
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortZbyModulus                    (A,first,last)                                           !:.:.:
  complex(dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) > ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) > ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortZbyModulus(a, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortZbyModulus(a, j+1  , last)

 end subroutine     array1_reversequicksortZbyModulus
!-----------------------------------------------------------------------------------------------------------------------------------
!Simple sort routine, suitable for small arrays
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_sortZbyModulus                                  (C)                                                      !:.:.:
  complex(dp), dimension(:)                       :: C
  complex(dp)                                     :: a
  integer                                         :: first, last
  integer                                         :: i,j,n

  n = size(C,1)
  do j=2, n
   a=C(j)
   do i=j-1,1,-1
    if (ABS(C(i))>=ABS(a)) goto 10
    C(i+1)=C(i)
   end do
   i=0
10 C(i+1)=a
  end do

 end subroutine     array1_sortZbyModulus
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortZbyRealPart                          (A,first,last)                                           !:.:.:
  complex(dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (REAL(a(i), kind=dp) < REAL(x   , kind=dp)   )
    i=i+1
   end do
   do while (REAL(x   , kind=dp) < REAL(a(j), kind=dp)   )
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortZbyRealPart(a, first, i-1 )
  if (j+1   < last)  call array1_quicksortZbyRealPart(a, j+1  , last)

 end subroutine     array1_quicksortZbyRealPart
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortZbyRealPart                   (A,first,last)                                           !:.:.:
  complex(dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (REAL(a(i), kind=dp) > REAL(x   , kind=dp)   )
    i=i+1
   end do
   do while (REAL(x   , kind=dp) > REAL(a(j), kind=dp)   )
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortZbyRealPart(a, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortZbyRealPart(a, j+1  , last)

 end subroutine     array1_reversequicksortZbyRealPart
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortZbyImagPart                          (A,first,last)                                           !:.:.:
  complex(dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (AIMAG(a(i)) < AIMAG(x   ) )
    i=i+1
   end do
   do while (AIMAG(x   ) < AIMAG(a(j)) )
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortZbyImagPart(a, first, i-1 )
  if (j+1   < last)  call array1_quicksortZbyImagPart(a, j+1  , last)

 end subroutine     array1_quicksortZbyImagPart
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortZbyImagPart                   (A,first,last)                                           !:.:.:
  complex(dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (AIMAG(a(i)) > AIMAG(x   )   )
    i=i+1
   end do
   do while (AIMAG(x   ) > AIMAG(a(j))   )
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortZbyImagPart(a, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortZbyImagPart(a, j+1  , last)

 end subroutine     array1_reversequicksortZbyImagPart
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortDbyModulus_d                         (A,first,last)                                           !:.:.:
  real   (dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) < ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) < ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortDbyModulus_d(a, first, i-1 )
  if (j+1   < last)  call array1_quicksortDbyModulus_d(a, j+1  , last)

 end subroutine     array1_quicksortDbyModulus_d
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortDbyModulus_d                  (A,first,last)                                           !:.:.:
  real   (dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) > ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) > ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortDbyModulus_d(a, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortDbyModulus_d(a, j+1  , last)

 end subroutine     array1_reversequicksortDbyModulus_d
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortDbyValue_d                           (A,first,last)                                           !:.:.:
  real   (dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while ((a(i)) < (x)   )
    i=i+1
   end do
   do while ((x   ) < (a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortDbyValue_d(a, first, i-1 )
  if (j+1   < last)  call array1_quicksortDbyValue_d(a, j+1  , last)

 end subroutine     array1_quicksortDbyValue_d
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortDbyValue_d                    (A,first,last)                                           !:.:.:
  real   (dp), dimension(:)                       :: A
  integer                                         :: first, last
  integer                                         :: i,j
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while ((a(i)) > (x)   )
    i=i+1
   end do
   do while ((x   ) > (a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t = a(i);  a(i) = a(j);  a(j) = t
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortDbyValue_d(a, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortDbyValue_d(a, j+1  , last)

 end subroutine     array1_reversequicksortDbyValue_d
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortZbyModulus_pos                       (A,pos,first,last)                                       !:.:.:
  complex(dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) < ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) < ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortZbyModulus_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_quicksortZbyModulus_pos(a, pos, j+1  , last)

 end subroutine     array1_quicksortZbyModulus_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortZbyModulus_pos                (A,pos,first,last)                                       !:.:.:
  complex(dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) > ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) > ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortZbyModulus_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortZbyModulus_pos(a, pos, j+1  , last)

 end subroutine     array1_reversequicksortZbyModulus_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortZbyRealPart_pos                      (A,pos,first,last)                                       !:.:.:
  complex(dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (REAL(a(i), kind=dp) < REAL(x   , kind=dp)   )
    i=i+1
   end do
   do while (REAL(x   , kind=dp) < REAL(a(j), kind=dp)   )
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortZbyRealPart_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_quicksortZbyRealPart_pos(a, pos, j+1  , last)

 end subroutine     array1_quicksortZbyRealPart_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortZbyRealPart_pos               (A,pos,first,last)                                       !:.:.:
  complex(dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (REAL(a(i), kind=dp) > REAL(x   , kind=dp)   )
    i=i+1
   end do
   do while (REAL(x   , kind=dp) > REAL(a(j), kind=dp)   )
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortZbyRealPart_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortZbyRealPart_pos(a, pos, j+1  , last)

 end subroutine     array1_reversequicksortZbyRealPart_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortZbyImagPart_pos                      (A,pos,first,last)                                       !:.:.:
  complex(dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (AIMAG(a(i)) < AIMAG(x   ) )
    i=i+1
   end do
   do while (AIMAG(x   ) < AIMAG(a(j)) )
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortZbyImagPart_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_quicksortZbyImagPart_pos(a, pos, j+1  , last)

 end subroutine     array1_quicksortZbyImagPart_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! Re(A) sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortZbyImagPart_pos               (A,pos,first,last)                                       !:.:.:
  complex(dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  complex(dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (AIMAG(a(i)) > AIMAG(x   )   )
    i=i+1
   end do
   do while (AIMAG(x   ) > AIMAG(a(j))   )
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortZbyImagPart_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortZbyImagPart_pos(a, pos, j+1  , last)

 end subroutine     array1_reversequicksortZbyImagPart_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortDbyModulus_d_pos                     (A,pos,first,last)                                       !:.:.:
  real   (dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) < ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) < ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortDbyModulus_d_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_quicksortDbyModulus_d_pos(a, pos, j+1  , last)

 end subroutine     array1_quicksortDbyModulus_d_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortDbyModulus_d_pos              (A,pos,first,last)                                       !:.:.:
  real   (dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while (ABS(a(i)) > ABS(x)   )
    i=i+1
   end do
   do while (ABS(x   ) > ABS(a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortDbyModulus_d_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortDbyModulus_d_pos(a, pos, j+1  , last)

 end subroutine     array1_reversequicksortDbyModulus_d_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from smallest to largest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_quicksortDbyValue_d_pos                       (A,pos,first,last)                                       !:.:.:
  real   (dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while ((a(i)) < (x)   )
    i=i+1
   end do
   do while ((x   ) < (a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_quicksortDbyValue_d_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_quicksortDbyValue_d_pos(a, pos, j+1  , last)

 end subroutine     array1_quicksortDbyValue_d_pos
!-----------------------------------------------------------------------------------------------------------------------------------
! From quicksort.f -*-f90-*-, Author: t-nissie, License: GPLv3, Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! |A| sorted from largest to smallest
!-----------------------------------------------------------------------------------------------------------------------------------
 recursive subroutine array1_reversequicksortDbyValue_d_pos                (A,pos,first,last)                                       !:.:.:
  real   (dp), dimension(:)                       :: A
  integer    , dimension(:)                       :: pos
  integer                                         :: first, last
  integer                                         :: i,j,k
  real   (dp)                                     :: x,t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
   do while ((a(i)) > (x)   )
    i=i+1
   end do
   do while ((x   ) > (a(j)))
    j=j-1
   end do
   if (i >= j) exit
   t =   a(i);   a(i) =   a(j);   a(j) = t
   k = pos(i); pos(i) = pos(j); pos(j) = k
   i=i+1
   j=j-1
  end do
  if (first < i-1 )  call array1_reversequicksortDbyValue_d_pos(a, pos, first, i-1 )
  if (j+1   < last)  call array1_reversequicksortDbyValue_d_pos(a, pos, j+1  , last)

 end subroutine     array1_reversequicksortDbyValue_d_pos

!-----------------------------------------------------------------------------------------------------------------------------------
!
! Random sorting of arrays using the Fisher-Yates Shuffle algorithm
!
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_sort_i                                          (array)                                                  !:.:.:
  integer, dimension(:), intent(inout) :: array
  integer  :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  n = size(array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end subroutine     random_sort_i
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_sort_d                                          (array)                                                  !:.:.:
  real(dp), dimension(:), intent(inout) :: array
  real(dp) :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  n = size(array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end subroutine     random_sort_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_sort_c                                          (array)                                                  !:.:.:
  complex(dp), dimension(:), intent(inout) :: array
  complex(dp) :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  n = size(array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end subroutine     random_sort_c
!-----------------------------------------------------------------------------------------------------------------------------------
 function           random_list_n                                          (n1,n2)                  result(array)                   !:.:.:
  integer, allocatable, dimension(:) :: array
  integer                            :: temp
  integer, intent(in)                :: n1
  integer, intent(in), optional      :: n2
  !------------------
  integer  :: nmin, nmax
  integer  :: i, j, n
  real(dp) :: r
  character(1000) :: err
  !------------------

  if(present(n2))then
   nmin = n1
   nmax = n2
  else
   nmin = 1
   nmax = n1
  end if

  n = nmax - nmin + 1! Get the length of the array
  if( n < 1 )then
   write(err,*) 'random_sort: random_list_n: n<1 nmin= ',nmin,' nmax= ',nmax
   call matrix_error(err)
  else
   allocate(array(n))
   array = [(i, i=nmin,nmax)]
  end if

  ! Iterate through the array, starting from the first element
  do i = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end function       random_list_n
!-----------------------------------------------------------------------------------------------------------------------------------
 function           random_sort_fun_i                                      (array_in)               result(array)                   !:.:.:
  integer, dimension(:), intent(in)  :: array_in
  integer, dimension(:), allocatable :: array
  integer  :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  array = array_in
  n = size(array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end function       random_sort_fun_i
!-----------------------------------------------------------------------------------------------------------------------------------
 function           random_sort_fun_d                                      (array_in)               result(array)                   !:.:.:
  real(dp), dimension(:), intent(in)  :: array_in
  real(dp), dimension(:), allocatable :: array
  real(dp) :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  array = array_in
  n     = size(array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i  = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end function       random_sort_fun_d
!-----------------------------------------------------------------------------------------------------------------------------------
 function           random_sort_fun_c                                      (array_in)               result(array)                   !:.:.:
  complex(dp), dimension(:), intent(in)  :: array_in
  complex(dp), dimension(:), allocatable :: array
  complex(dp) :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  array = array_in
  n     = size(array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i  = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp     = array(i)
   array(i) = array(j)
   array(j) = temp
  end do
  
 end function       random_sort_fun_c
!-----------------------------------------------------------------------------------------------------------------------------------
 function           random_sort_fun_char                                   (array_in)               result(array)                   !:.:.:
  character(*), intent(in)  :: array_in
  character(:), allocatable :: array
  character(1) :: temp
  !------------------
  integer  :: i, j, n
  real(dp) :: r
  !------------------
  array = trim(array_in)
  n     = len (array)! Get the length of the array
  
  ! Iterate through the array, starting from the first element
  do i  = 1, n - 1
   
   call random_number(r)  
   j = i + int(r * (n - i + 1))  ! Generate a random index 'j' between i and n (inclusive)
   
   ! Swap elements at index i and j
   temp       = array(i:i)
   array(i:i) = array(j:j)
   array(j:j) = temp
  end do
  
 end function       random_sort_fun_char



!...................................................................................................................................
!.......................... random_number interface: ...............................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_complex_scalar                           (z)                                                      !:.:.:
  complex(dp), intent(in out)                     :: z
  real(dp)   , dimension(2)                       :: x

  call random_number(x)
  
  z = CMPLX( x(1),x(2),kind=dp)

 end subroutine     random_number_complex_scalar
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array3                                   (C)                                                      !:.:.:
  complex(dp), dimension(:,:,:)                         :: C
  real   (dp), dimension(size(C,1),size(C,2),size(C,3)) :: x,y
  
  call random_number(x); call random_number(y)
  C = CMPLX(x,y,kind=dp)

 end subroutine     random_number_array3
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array2                                   (C)                                                      !:.:.:
  complex(dp), dimension(:,:)                     :: C
  real   (dp), dimension(size(C,1),size(C,2))     :: x,y
  
  call random_number(x); call random_number(y)
  C = CMPLX(x,y,kind=dp)

 end subroutine     random_number_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array1                                   (C)                                                      !:.:.:
  complex(dp), dimension(:)                       :: C
  real   (dp), dimension(size(C,1))               :: x,y
  
  call random_number(x); call random_number(y)
  C = CMPLX(x,y,kind=dp)

 end subroutine     random_number_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_complex_scalar_gaussian                  (z,sigma)                                                !:.:.:
  complex(dp),                intent(out)         :: z
  real   (dp),                intent(in)          :: sigma
  real   (dp)                                     :: x,y

  call random_number(x); call random_number(y)
  z = sigma * sqrt(-2.0_dp * log(x)) * CMPLX(cos(TWOPI*y),sin(TWOPI*y),kind=dp)
  
 end subroutine     random_number_complex_scalar_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_real_scalar_gaussian                     (r,sigma)                                                !:.:.:
  real   (dp),                intent(out)         :: r
  real   (dp),                intent(in)          :: sigma
  real   (dp)                                     :: x,y

  call random_number(x);call random_number(y);
  r = sigma * sqrt(-2.0_dp * log(x)) * cos(TWOPI*y)
  
 end subroutine     random_number_real_scalar_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array3_gaussian                          (C,sigma)                                                !:.:.:
  complex(dp), dimension(:,:,:),intent(out)              :: C
  real   (dp),                  intent(in)               :: sigma
  real   (dp), dimension(size(C,1),size(C,2),size(C,3))  :: x,y
  
  call random_number(x); call random_number(y)
  C = sigma * sqrt(-2.0_dp * log(x)) * CMPLX(cos(TWOPI*y),sin(TWOPI*y),kind=dp)
  
 end subroutine     random_number_array3_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array3_gaussian_d                        (C,sigma)                                                !:.:.:
  real   (dp), dimension(:,:,:),intent(out)             :: C
  real   (dp),                  intent(in)              :: sigma
  real   (dp), dimension(size(C,1),size(C,2),size(C,3)) :: x,y

  call random_number(x); call random_number(y)
  C = sigma * sqrt(-2.0_dp * log(x)) * cos(TWOPI*y)

 end subroutine random_number_array3_gaussian_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array2_gaussian                          (C,sigma)                                                !:.:.:
  complex(dp), dimension(:,:),intent(out)         :: C
  real   (dp),                intent(in)          :: sigma
  real   (dp), dimension(size(C,1),size(C,2))     :: x,y
  
  call random_number(x); call random_number(y)
  C = sigma * sqrt(-2.0_dp * log(x)) * CMPLX(cos(TWOPI*y),sin(TWOPI*y),kind=dp)
  
 end subroutine     random_number_array2_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array2_gaussian_d                        (C,sigma)                                                !:.:.:
  real   (dp), dimension(:,:),intent(out)         :: C
  real   (dp),                intent(in)          :: sigma
  real   (dp), dimension(size(C,1),size(C,2))     :: x,y

  call random_number(x); call random_number(y)
  C = sigma * sqrt(-2.0_dp * log(x)) * cos(TWOPI*y)

 end subroutine random_number_array2_gaussian_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array1_gaussian                          (C,sigma)                                                !:.:.:
  complex(dp), dimension(:)  ,intent(out)         :: C
  real   (dp),                intent(in)          :: sigma
  real   (dp), dimension(size(C,1))               :: x,y
  
  call random_number(x); call random_number(y)
  C = sigma * sqrt(-2.0_dp * log(x)) * CMPLX(cos(TWOPI*y),sin(TWOPI*y),kind=dp)
  
 end subroutine     random_number_array1_gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         random_number_array1_gaussian_d                        (C,sigma)                                                !:.:.:
  real   (dp), dimension(:)  ,intent(out)         :: C
  real   (dp),                intent(in)          :: sigma
  real   (dp), dimension(size(C,1))               :: x,y

  call random_number(x); call random_number(y)
  C = sigma * sqrt(-2.0_dp * log(x)) * cos(TWOPI*y)

 end subroutine     random_number_array1_gaussian_d

!...................................................................................................................................
!.......................... array2 matmul operator overload ........................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_array2                                   (C1,C2)                  result(C3)                      !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C1
  complex(dp), dimension(:,:), intent(in)         :: C2
  complex(dp), dimension(size(C1,1),size(C2,2))   :: C3

  C3 = lmatmul(C1,C2)

 end  function      array2_matmul_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_array2_d                                 (C1,C2)                  result(C3)                      !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C1
  real   (dp), dimension(:,:), intent(in)         :: C2
  complex(dp), dimension(size(C1,1),size(C2,2))   :: C3
  complex(dp), dimension(size(C2,1),size(C2,2))   :: C4

  C4 = C2  ! real -> complex
  C3 = lmatmul(C1,C4)

 end  function      array2_matmul_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_d_matmul_array2                                 (C1,C2)                  result(C3)                      !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C1
  complex(dp), dimension(:,:), intent(in)         :: C2
  complex(dp), dimension(size(C1,1),size(C2,2))   :: C3
  complex(dp), dimension(size(C1,1),size(C1,2))   :: C4

  C4 = C1  ! real -> complex
  C3 = lmatmul(C4,C2)

 end  function      array2_d_matmul_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_d_matmul_array2_d                               (C1,C2)                  result(C3)                      !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C1
  real   (dp), dimension(:,:), intent(in)         :: C2
  real   (dp), dimension(size(C1,1),size(C2,2))   :: C3

  C3 = lmatmul(C1,C2)

 end  function      array2_d_matmul_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_array1                                   (C1,v2)                  result(v3)                      !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C1
  complex(dp), dimension(:  ), intent(in)         :: v2
  complex(dp), dimension(size(C1,1))              :: v3

  v3 = lmatmul( A = C1 , v = v2 )

 end  function      array2_matmul_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_matmul_array2                                   (v1,C2)                  result(v3)                      !:.:.:
  complex(dp), dimension(:  ), intent(in)         :: v1
  complex(dp), dimension(:,:), intent(in)         :: C2
  complex(dp), dimension(size(C2,2))              :: v3

  v3 = lmatmul( v = v1 , A = C2  , type='T')

 end  function      array1_matmul_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_d_matmul_array1_d                               (C1,v2)                  result(v3)                      !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C1
  real   (dp), dimension(:  ), intent(in)         :: v2
  real   (dp), dimension(size(C1,1))              :: v3

  v3 = lmatmul( A = C1 , v = v2 )

 end  function      array2_d_matmul_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_d_matmul_array2_d                               (v1,C2)                  result(v3)                      !:.:.:
  real   (dp), dimension(:  ), intent(in)         :: v1
  real   (dp), dimension(:,:), intent(in)         :: C2
  real   (dp), dimension(size(C2,2))              :: v3

  v3 = lmatmul( v = v1 , A = C2  , type='T')

 end  function      array1_d_matmul_array2_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_matmul_array1_d                                 (C1,v2)                  result(v3)                      !:.:.:
  complex(dp), dimension(:,:), intent(in)         :: C1
  real   (dp), dimension(:  ), intent(in)         :: v2
  complex(dp), dimension(size(C1,1))              :: v3
  complex(dp), dimension(size(v2,1))              :: v4

  v4 = v2

  v3 = lmatmul( A = C1 , v = v4 )

 end  function      array2_matmul_array1_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_d_matmul_array2                                 (v1,C2)                  result(v3)                      !:.:.:
  real   (dp), dimension(:  ), intent(in)         :: v1
  complex(dp), dimension(:,:), intent(in)         :: C2
  complex(dp), dimension(size(C2,2))              :: v3
  complex(dp), dimension(size(v1,1))              :: v4

  v4 = v1

  v3 = lmatmul( v = v4 , A = C2  , type='T')

 end  function      array1_d_matmul_array2
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_d_matmul_array1                                 (C1,v2)                  result(v3)                      !:.:.:
  real   (dp), dimension(:,:), intent(in)         :: C1
  complex(dp), dimension(:  ), intent(in)         :: v2
  complex(dp), dimension(size(C1,1))              :: v3
  complex(dp), dimension(size(C1,1),size(C1,2))   :: C4

  C4 = C1

  v3 = lmatmul( A = C4 , v = v2 )

 end  function      array2_d_matmul_array1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_matmul_array2_d                                 (v1,C2)                  result(v3)                      !:.:.:
  complex(dp), dimension(:  ), intent(in)         :: v1
  real   (dp), dimension(:,:), intent(in)         :: C2
  complex(dp), dimension(size(C2,2))              :: v3
  complex(dp), dimension(size(C2,1),size(C2,2))   :: C4

  C4 = C2

  v3 = lmatmul( v = v1 , A = C4  , type='T')

 end  function      array1_matmul_array2_d


!...................................................................................................................................
!.......................... utilities: .............................................................................................!:.:.:
!...................................................................................................................................

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call array2_print(C(:,:), unit=f_mout,fmt='G28.17',form='Mathematica',ips=12,ipe=24,jps=55,jpe=65)
!
!C(:,:) must be allocatable (so that we can pass the array bounds). If you wish to print a non allocatable array Cstatic,  allocate(C,source=Cstatic)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_print                                           (C,unit,fmt,form,name,ips,ipe,jps,jpe)                   !:.:.:
  complex (dp),allocatable, intent(in)            :: C(:,:)                                                                         ! allocatable: Fortran 2003 feature: passes lbound,ubound of C
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe,jps,jpe
  character(*),  optional, intent(in)             :: name
  character(name_len)                             :: namem
  character(25)                                   :: fmts,fmtn,myform,eol
  character(math_len)                             :: nleft,nright
  integer                                         :: n,m,i,j,is,ie,js,je,un
  character(2000)                                 :: s  ! Gamma matrix format
  character(3)                                    :: e
  real    (dp), parameter                         :: cut = 1.0e-10_dp

  un    = f_mout;   fmtn    = 'G28.17'
  if(present(unit)) un      = unit
  if(present(fmt )) fmtn    = fmt

  m     = size  (C,1); n    = size  (C,2)
  is    = lbound(C,1); ie   = ubound(C,1)
  js    = lbound(C,2); je   = ubound(C,2)

  namem = 'AMAT'
  if(present(name)) namem = name

  if(present(ips) .and. ips >= is) is = ips  ! if asked, print only the ips:ipe, jps:jpe section of the arrays
  if(present(ipe) .and. ipe <= ie) ie = ipe
  if(present(jps) .and. jps >= js) js = jps
  if(present(jpe) .and. jpe <= je) je = jpe

  myform = 'default'
  if(present(form)) myform = form

  select case (myform(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('d','D')
   write(un,'(A,I6,A3,I6,A,I6,A3,I6,A2,I6,A3,I6,A2)') &
        '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
        '#pmstart: ' // TRIM(namem) // ' =[                    (is a ',m,' x ',n,'  matrix) (printing bounds: ',    &
        is,' : ',ie,'  ',js,' : ',je,' )'

   write(fmts,'(A,I12,2A)') '(',2*n, TRIM(fmtn),')'
   do i = is, ie
    write(un,fmts) (C(i,j), j = js, je)
   end do;

   write(un,'(A           )') &
       '#pmend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('m','M') ! Mathematica Format
   nleft=Mathematica_nleft; nright=Mathematica_nright ! from matrix_mod_common: a mathematica number in exponential notation must be between nleft and nright
   write(un,'(A)') trim(namem) // '={'

   !--------------------------------------------------------------------------------------------------------------------------------
   do  i = is, ie
    write(un,'(A)') '{'
    do j = js, je - 1
     write(un,'(A,'//trim(fmtn)//',A,'//trim(fmtn)//',A)')  &
          trim(nleft),REAL(C(i,j), kind=dp),trim(nright)//' + I ('//trim(nleft),IMAG(C(i,j)),trim(nright)//'),'
    end do
    eol             = ') },'
    if(i == ie) eol = ') }'
    write (un,'(A,'//trim(fmtn)//',A,'//trim(fmtn)//',A)')  &
         trim(nleft), REAL(C(i,j), kind=dp),trim(nright)//' + I ('//trim(nleft),IMAG(C(i,j)),trim(nright)//trim(eol)
   end do !i = is, ie
   !--------------------------------------------------------------------------------------------------------------------------------
   write(un,'(A)') '};'
   !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I') ! Info form
   write(un,'(A,I6,A,I6,A,I6,A,I6,A,I6,A,I6,A)') '#pminfo: ' // TRIM(namem) // ' is a ',m,' x ',n,'  matrix, printing bounds: ',   &
        is,' : ',ie,'  ',js,' : ',je,' '
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('g','G')
   write(un,'(A,I6,A3,I6,A,I6,A3,I6,A2,I6,A3,I6,A2)') &
        '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
        '#pmstart: ' // TRIM(namem) // ' =[                    (is a ',m,' x ',n,'  matrix) (printing bounds: ',    &
        is,' : ',ie,'  ',js,' : ',je,' )'

   do  i = is , ie
    s    = ''
    do j = js, je
     if     ( abs(C(i,j)-ZERO) < cut)then
      e = '  0'
     else if( abs(C(i,j)- ONE) < cut)then
      e = '  1'
     else if( abs(C(i,j)+ ONE) < cut)then
      e = ' -1'
     else if( abs(C(i,j)- IMU) < cut)then
      e = '  i'
     else if( abs(C(i,j)+ IMU) < cut)then
      e = ' -i'
     else
      e = ' NN'
     end if
     s = trim(s) // e
    end do
    print '(A)', trim(s)
   end do

   write(un,'(A           )') &
       '#pmend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
   !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call matrix_error('matrix_print: unknown form for printing MAT')
  !---------------------------------------------------------------------------------------------------------------------------------
  end select

 end subroutine array2_print
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call array1_print(C(:)), unit=f_mout,fmt='G28.17',form='Mathematica',ips=12,ipe=24)
!
!C(:)   must be allocatable (so that we can pass the array bounds). If you wish to print a non allocatable array Cstatic,  allocate(C,source=Cstatic)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_print                                           (C,unit,fmt,form,name,ips,ipe)                           !:.:.:
  complex(dp),allocatable, intent(in)             :: C(:)                                                                           ! allocatable: Fortran 2003 feature: passes lbound,ubound of C
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe
  character(*),  optional, intent(in)             :: name
  character(name_len)                             :: namem
  character(25)                                   :: fmts,fmtn,myform,eol
  character(math_len)                             :: nleft,nright
  integer                                         :: n,m,i,is,ie,un

  un    = f_mout;   fmtn    = 'G28.17'
  if(present(unit)) un      = unit
  if(present(fmt )) fmtn    = fmt

  m     = size  (C,1);
  is    = lbound(C,1); ie   = ubound(C,1)

  namem = 'AMAT'
  if(present(name)) namem = name

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

   write(fmts,'(A,I12,2A)') '(',2, TRIM(fmtn),')'
   do i = is, ie
    write(un,fmts) C(i) 
   end do;

   write(un,'(A           )') &
       '#pvend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('m','M') ! Mathematica Format
   nleft=Mathematica_nleft; nright=Mathematica_nright ! from module_mod_common, a mathematica number in exponential notation must be between nleft and nright
   write(un,'(A)') trim(namem) // '={'

   !--------------------------------------------------------------------------------------------------------------------------------
   do  i = is, ie

    eol = '),'
    if(i == ie) eol = ')'
    write(un,'(A,'//trim(fmtn)//',A,'//trim(fmtn)//',A)') &
         trim(nleft),REAL(C(i), kind=dp),trim(nright)//' + I ('//trim(nleft),AIMAG(C(i)),trim(nright)//trim(eol)

   end do !i = is, ie
   !--------------------------------------------------------------------------------------------------------------------------------
   write(un,'(A)') '};'
   !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I') ! info mode
   write(un,'(A,I6,A,I6,A,I6,A)') '#pminfo: ' // TRIM(namem) // ' is a ',m,' vector, printing bounds: ',   &
        is,' : ',ie,'  '
   !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call matrix_error('matrix_print: unknown form for printing MAT')
  !---------------------------------------------------------------------------------------------------------------------------------
  end select

 end subroutine array1_print
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call arra2_print(C(:,:), unit=f_mout,fmt='G28.17',form='Mathematica',ips=12,ipe=24,jps=55,jpe=65)
!
!C(:,:) must be allocatable (so that we can pass the array bounds). If you wish to print a non allocatable array Cstatic,  allocate(C,source=Cstatic)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_print_d                                         (C,unit,fmt,form,name,ips,ipe,jps,jpe)                   !:.:.:   
  real   (dp),allocatable, intent(in)             :: C(:,:)                                                                         ! allocatable: Fortran 2003 feature: passes lbound,ubound of C
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe,jps,jpe
  character(*),  optional, intent(in)             :: name
  character(name_len)                             :: namem
  character(25)                                   :: fmts,fmtn,myform,eol
  character(math_len)                             :: nleft,nright
  integer                                         :: n,m,i,j,is,ie,js,je,un

  un    = f_mout;   fmtn    = 'G28.17'
  if(present(unit)) un      = unit
  if(present(fmt )) fmtn    = fmt

  m     = size  (C,1); n    = size  (C,2)
  is    = lbound(C,1); ie   = ubound(C,1)
  js    = lbound(C,2); je   = ubound(C,2)

  namem = 'AMAT'
  if(present(name)) namem = name

  if(present(ips) .and. ips >= is) is = ips  ! if asked, print only the ips:ipe, jps:jpe section of the arrays
  if(present(ipe) .and. ipe <= ie) ie = ipe
  if(present(jps) .and. jps >= js) js = jps
  if(present(jpe) .and. jpe <= je) je = jpe

  myform = 'default'
  if(present(form)) myform = form

  select case (myform(1:1))
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('d','D')
   write(un,'(A,I6,A3,I6,A,I6,A3,I6,A2,I6,A3,I6,A2)') &
        '#------------------------------------------------------------------------------------------'//NEW_LINE('A')//&
        '#pmstart: ' // TRIM(namem) // ' =[                    (is a ',m,' x ',n,'  matrix) (printing bounds: ',    &
        is,' : ',ie,'  ',js,' : ',je,' )'

   write(fmts,'(A,I12,2A)') '(',  n, TRIM(fmtn),')'
   do i = is, ie
    write(un,fmts) (C(i,j), j = js, je)
   end do;

   write(un,'(A           )') &
       '#pmend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('m','M') ! Mathematica Format
   nleft=Mathematica_nleft; nright=Mathematica_nright ! from matrix_mod_common, a mathematica number in exponential notation must be between nleft and nright
   write(un,'(A)') trim(namem) // '={'

   !--------------------------------------------------------------------------------------------------------------------------------
   do  i = is, ie
    write(un,'(A)') '{'
    do j = js, je - 1
     write(un,'(A,'//trim(fmtn)//                   ',A)')  &
          trim(nleft),C(i,j),trim(nright)//','
    end do
    eol             = '},'
    if(i == ie) eol = '}'
    write (un,'(A,'//trim(fmtn)//                   ',A)')  &
         trim(nleft), C(i,j),trim(nright)//trim(eol)
   end do !i = is, ie
   !--------------------------------------------------------------------------------------------------------------------------------
   write(un,'(A)') '};'
   !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I') ! Info form
   write(un,'(A,I6,A,I6,A,I6,A,I6,A,I6,A,I6,A)') '#pminfo: ' // TRIM(namem) // ' is a ',m,' x ',n,'  matrix, printing bounds: ',   &
        is,' : ',ie,'  ',js,' : ',je,'  '
   !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call matrix_error('matrix_print_d: unknown form for printing MAT')
  !---------------------------------------------------------------------------------------------------------------------------------
  end select

 end subroutine array2_print_d

!-----------------------------------------------------------------------------------------------------------------------------------
!Usage: call array1_print_d(C(:)), unit=f_mout,fmt='G28.17',form='Mathematica',ips=12,ipe=24)
!
!C(:)   must be allocatable (so that we can pass the array bounds). If you wish to print a non allocatable array Cstatic,  allocate(C,source=Cstatic)
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_print_d                                         (C,unit,fmt,form,name,ips,ipe)                           !:.:.:
  real   (dp),allocatable, intent(in)             :: C(:)                                                                           ! allocatable: Fortran 2003 feature: passes lbound,ubound of C
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe
  character(*),  optional, intent(in)             :: name
  character(name_len)                             :: namem
  character(25)                                   :: fmts,fmtn,myform,eol
  character(math_len)                             :: nleft,nright
  integer                                         :: n,m,i,is,ie,un

  un    = f_mout;   fmtn    = 'G28.17'
  if(present(unit)) un      = unit
  if(present(fmt )) fmtn    = fmt

  m     = size  (C,1);
  is    = lbound(C,1); ie   = ubound(C,1)

  namem = 'AMAT'
  if(present(name)) namem = name

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

   write(fmts,'(A)') '('//TRIM(fmtn)//')'
   do i = is, ie
    write(un,fmts) C(i) 
   end do;

   write(un,'(A           )') &
       '#pvend: ]'                                                                                   //NEW_LINE('A')//&
       '#------------------------------------------------------------------------------------------'
  !---------------------------------------------------------------------------------------------------------------------------------
  case ('m','M') ! Mathematica Format
   nleft=Mathematica_nleft; nright=Mathematica_nright ! from matrix_mod_common, a mathematica number in exponential notation must be between nleft and nright
   write(un,'(A)') trim(namem) // '={'

   !--------------------------------------------------------------------------------------------------------------------------------
   do  i = is, ie

    eol = ','
    if(i == ie) eol = ''
    write(un,'(A,'//trim(fmtn)//',A)') &
          trim(nleft),C(i),trim(nright)//trim(eol)

   end do !i = is, ie
   !--------------------------------------------------------------------------------------------------------------------------------
   write(un,'(A)') '};'
   !---------------------------------------------------------------------------------------------------------------------------------
  case ('i','I') ! info mode
   write(un,'(A,I6,A,I6,A,I6,A)') '#pminfo: ' // TRIM(namem) // ' is a ',m,' vector, printing bounds: ',   &
        is,' : ',ie,' '
   !---------------------------------------------------------------------------------------------------------------------------------
  case default
   call matrix_error('matrix_print: unknown form for printing MAT')
  !---------------------------------------------------------------------------------------------------------------------------------
  end select

 end subroutine array1_print_d
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_print_nonallocatable                            (C,unit,fmt,form,name,ips,ipe,jps,jpe)                   !:.:.:
  complex(dp)            , intent(in)             :: C(:,:)                                                                         
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe,jps,jpe
  character(*),  optional, intent(in)             :: name
  complex(dp),   allocatable                      :: Ca(:,:)             

  ALLOCATE(Ca,source=C)

  call array2_print(Ca,unit,fmt,form,name,ips,ipe,jps,jpe) 

  DEALLOCATE(Ca)

 end subroutine     array2_print_nonallocatable
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_print_nonallocatable_d                          (C,unit,fmt,form,name,ips,ipe,jps,jpe)                   !:.:.:
  real   (dp)            , intent(in)             :: C(:,:)                                                                         
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe,jps,jpe
  character(*),  optional, intent(in)             :: name
  real   (dp),   allocatable                      :: Ca(:,:)             

  ALLOCATE(Ca,source=C)

  call array2_print_d(Ca,unit,fmt,form,name,ips,ipe,jps,jpe) 

  DEALLOCATE(Ca)

 end subroutine     array2_print_nonallocatable_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_print_nonallocatable                            (C,unit,fmt,form,name,ips,ipe,jps,jpe)                   !:.:.:
  complex(dp)            , intent(in)             :: C(:)                                                                         
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe,jps,jpe
  character(*),  optional, intent(in)             :: name
  complex(dp),   allocatable                      :: Ca(:)             

  ALLOCATE(Ca,source=C)

  call array1_print(Ca,unit,fmt,form,name,ips,ipe) 

  DEALLOCATE(Ca)

 end subroutine     array1_print_nonallocatable
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_print_nonallocatable_d                          (C,unit,fmt,form,name,ips,ipe,jps,jpe)                   !:.:.:
  real   (dp)            , intent(in)             :: C(:)                                                                         
  integer     ,  optional, intent(in)             :: unit
  character(*),  optional, intent(in)             :: fmt
  character(*),  optional, intent(in)             :: form
  integer     ,  optional, intent(in)             :: ips,ipe,jps,jpe
  character(*),  optional, intent(in)             :: name
  real   (dp),   allocatable                      :: Ca(:)             

  ALLOCATE(Ca,source=C)

  call array1_print_d(Ca,unit,fmt,form,name,ips,ipe) 

  DEALLOCATE(Ca)

 end subroutine     array1_print_nonallocatable_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array3_save_matrix                                     (C,unit,fmt)                                             !:.:.:
  complex(dp) , dimension(:,:,:), intent(in)      :: C
  integer     , optional        , intent(in)      :: unit
  character(*), optional        , intent(in)      :: fmt
  character(20)                                   :: fmtout
  integer                                         :: n,m,k,i,j,l,un

  un     = f_mout
  fmtout = 'G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt

  m = size(C,1)  ; n = size(C,2) ; k = size(C,3)

  do   l = 1, k  !the outermost index is the last one: save as A^l_ij
   do  i = 1, m
    do j = 1, n
     write(un,'(2'//fmtout//')') C(i,j,l)
    end do
   end do
  end do

 end subroutine     array3_save_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_save_matrix                                     (C,unit,fmt)                                             !:.:.:
  complex(dp) , dimension(:,:), intent(in)        :: C
  integer     , optional      , intent(in)        :: unit
  character(*), optional      , intent(in)        :: fmt
  character(20)                                   :: fmtout
  integer                                         :: n,m,i,j,un

  un     = f_mout
  fmtout = 'G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt

  m = size(C,1)  ; n = size(C,2)

  do  i = 1, m
   do j = 1, n
    write(un,'(2'//fmtout//')') C(i,j)
   end do
  end do

 end subroutine     array2_save_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_save_matrix                                     (C,unit,fmt)                                             !:.:.:
  complex(dp) , dimension(:)  , intent(in)        :: C
  integer     , optional      , intent(in)        :: unit
  character(*), optional      , intent(in)        :: fmt
  character(20)                                   :: fmtout
  integer                                         :: n,m,i,j,un

  un     = f_mout
  fmtout = 'G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt

  m = size(C,1)

  do  i = 1, m
   write(un,'(2'//fmtout//')') C(i)
  end do

 end subroutine     array1_save_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array3_save_matrix_d                                   (C,unit,fmt)                                             !:.:.:
  real   (dp) , dimension(:,:,:), intent(in)      :: C
  integer     , optional        , intent(in)      :: unit
  character(*), optional        , intent(in)      :: fmt
  character(20)                                   :: fmtout
  integer                                         :: n,m,k,i,j,l,un

  un     = f_mout
  fmtout = 'G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt

  m = size(C,1)  ; n = size(C,2) ; k = size(C,3)

  do   l = 1, k  !the outermost index is the last one: save as A^l_ij
   do  i = 1, m
    do j = 1, n
     write(un,'('//fmtout//')') C(i,j,l)
    end do
   end do
  end do

 end subroutine     array3_save_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_save_matrix_d                                   (C,unit,fmt)                                             !:.:.:
  real   (dp) , dimension(:,:), intent(in)        :: C
  integer     , optional      , intent(in)        :: unit
  character(*), optional      , intent(in)        :: fmt
  character(20)                                   :: fmtout
  integer                                         :: n,m,i,j,un

  un     = f_mout
  fmtout = 'G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt

  m = size(C,1)  ; n = size(C,2)

  do  i = 1, m
   do j = 1, n
    write(un,'('//fmtout//')') C(i,j)
   end do
  end do

 end subroutine     array2_save_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_save_matrix_d                                   (C,unit,fmt)                                             !:.:.:
  real   (dp) , dimension(:)  , intent(in)        :: C
  integer     , optional      , intent(in)        :: unit
  character(*), optional      , intent(in)        :: fmt
  character(20)                                   :: fmtout
  integer                                         :: n,m,i,j,un

  un     = f_mout
  fmtout = 'G28.17'
  if(present(unit)) un     = unit
  if(present(fmt )) fmtout = fmt

  m = size(C,1)

  do  i = 1, m
   write(un,'('//fmtout//')') C(i)
  end do

 end subroutine     array1_save_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array3_read_matrix                                     (C,unit)                                                 !:.:.:
  complex(dp) , dimension(:,:,:), intent(in out)  :: C
  integer     , optional      , intent(in)        :: unit
  integer                                         :: n,m,k,i,j,l,un
  real   (dp)                                     :: x, y

  un     = f_minput
  if(present(unit)) un     = unit

  m = size(C,1)  ; n = size(C,2) ; k = size(C,3)

  do   l = 1, k !the outermost index is the last one: read as A^l_ij
   do  i = 1, m
    do j = 1, n
     read (un,*)    x,y
     C(i,j,l) = CMPLX(x,y,kind=dp)
    end do
   end do
  end do

 end subroutine     array3_read_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_read_matrix                                     (C,unit)                                                 !:.:.:
  complex(dp) , dimension(:,:), intent(in out)    :: C
  integer     , optional      , intent(in)        :: unit
  integer                                         :: n,m,i,j,un
  real   (dp)                                     :: x, y

  un     = f_minput
  if(present(unit)) un     = unit

  m = size(C,1)  ; n = size(C,2)

  do  i = 1, m
   do j = 1, n
    read (un,*)    x,y
    C(i,j) = CMPLX(x,y,kind=dp)
   end do
  end do

 end subroutine     array2_read_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_read_matrix                                     (C,unit)                                                 !:.:.:
  complex(dp) , dimension(:)  , intent(in out)    :: C
  integer     , optional      , intent(in)        :: unit
  integer                                         :: n,m,i,j,un
  real   (dp)                                     :: x, y

  un     = f_minput
  if(present(unit)) un     = unit

  m = size(C,1)

  do  i = 1, m
   read (un,*)    x,y
   C(i) = CMPLX(x,y,kind=dp)
  end do

 end subroutine     array1_read_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array3_read_matrix_d                                   (C,unit)                                                 !:.:.:
  real   (dp) , dimension(:,:,:), intent(in out)  :: C
  integer     , optional      , intent(in)        :: unit
  integer                                         :: n,m,k,i,j,l,un

  un     = f_minput
  if(present(unit)) un     = unit

  m = size(C,1)  ; n = size(C,2) ; k = size(C,3)

  do   l = 1, k
   do  i = 1, m
    do j = 1, n
     read (un,*) C(i,j,l) 
    end do
   end do
  end do

 end subroutine     array3_read_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array2_read_matrix_d                                   (C,unit)                                                 !:.:.:
  real   (dp) , dimension(:,:), intent(in out)    :: C
  integer     , optional      , intent(in)        :: unit
  integer                                         :: n,m,i,j,un

  un     = f_minput
  if(present(unit)) un     = unit

  m = size(C,1)  ; n = size(C,2)

  do  i = 1, m
   do j = 1, n
    read (un,*) C(i,j) 
   end do
  end do

 end subroutine     array2_read_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine         array1_read_matrix_d                                   (C,unit)                                                 !:.:.:
  real   (dp) , dimension(:)  , intent(in out)    :: C
  integer     , optional      , intent(in)        :: unit
  integer                                         :: n,m,i,j,un

  un     = f_minput
  if(present(unit)) un     = unit

  m = size(C,1) 

  do  i = 1, m
   read (un,*) C(i) 
  end do

 end subroutine     array1_read_matrix_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array3_is_nan                                          (C)                      result(itis)                    !:.:.:
  complex(dp) , dimension(:,:,:), intent(in)      :: C
  logical                                         :: itis
  integer                                         :: i,j,k

  itis = .false.

  do concurrent(i=1:size(C,1),j=1:size(C,2),k=1:size(C,3))
   if( matrix_return_is_nan_complex8_z(C(i,j,k)) ) itis = .true.
  end do

 end  function      array3_is_nan
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_nan                                          (C)                      result(itis)                    !:.:.:
  complex(dp) , dimension(:,:), intent(in)        :: C
  logical                                         :: itis
  integer                                         :: i,j

  itis = .false.

  do concurrent(i=1:size(C,1),j=1:size(C,2))
   if( matrix_return_is_nan_complex8_z(C(i,j)) )   itis = .true.
  end do

 end  function       array2_is_nan
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function       array1_is_nan                                         (C)                      result(itis)                    !:.:.:
  complex(dp) , dimension(:)  , intent(in)        :: C
  logical                                         :: itis
  integer                                         :: i

  itis = .false.

  do concurrent(i=1:size(C,1))
   if( matrix_return_is_nan_complex8_z(C(i)) )     itis = .true.
  end do

 end  function      array1_is_nan
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array3_is_nan_d                                        (C)                      result(itis)                    !:.:.:
  real   (dp) , dimension(:,:,:), intent(in)      :: C
  logical                                         :: itis
  integer                                         :: i,j,k

  itis = .false.

  do concurrent(i=1:size(C,1),j=1:size(C,2),k=1:size(C,3))
   if( matrix_return_is_nan_real8_d(C(i,j,k)) ) itis = .true.
  end do

 end function       array3_is_nan_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array2_is_nan_d                                        (C)                      result(itis)                    !:.:.:
  real   (dp) , dimension(:,:), intent(in)        :: C
  logical                                         :: itis
  integer                                         :: i,j

  itis = .false.

  do concurrent(i=1:size(C,1),j=1:size(C,2))
   if( matrix_return_is_nan_real8_d(C(i,j  )) ) itis = .true.
  end do

 end  function      array2_is_nan_d
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      array1_is_nan_d                                        (C)                      result(itis)                    !:.:.:
  real   (dp) , dimension(:)  , intent(in)        :: C
  logical                                         :: itis
  integer                                         :: i

  itis = .false.

  do concurrent(i=1:size(C,1))
   if( matrix_return_is_nan_real8_d(C(i    )) ) itis = .true.
  end do

 end  function      array1_is_nan_d
 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module          matrix_mod_array
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================



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
