!...................................................................................................................................
! function name:      column 20,    ::                 : column 40   , comments: column 132.                  Function information: ! : . : . :   
! function arguments: column 75     function result(..): column 100
! complex  (dp), optional, intent(in)  ::
!          ^ 11  ^ 17      ^ 28        ^ 40
!...................................................................................................................................!:.:.:
!.......................... File: matrix_mod_matrix.f90 ............................................................................!:.:.:
!...................................................................................................................................!:.:.:
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
module              matrix_mod_matrix                                                                                               !:.:.:
 use                matrix_mod_common                                                                                               !:.:.:
 use                matrix_mod_array                                                                                                !:.:.:
 implicit none
 private; save
!-----------------------------------------------------------------------------------------------------------------------------------!:.:.:
 public                                 :: Matrix, DMatrix, Vector, DVector                                                         !:.:.:
 public                                 :: mmmult, mvmult , vmmult                                                                  !:.:.:
 public                                 :: random_number, mcmplx, real, aimag, conjg, transpose, hermitian, norm, symmetric         !:.:.:
 public                                 :: isNaN, sort, dot_product, maxval, minval, trace, trace2, trace2c, traceless              !:.:.:
 public                                 :: diagonal, diagonalMatrix, metadata_copy, inverse, determinant, eigenvalues, eigenvectors !:.:.:
 public                                 :: lndet, Pfaffian, lnPfaffian, isHermitian, isSymmetric, isAntiSymmetric                   !:.:.:
 public                                 :: traceless_set, hermitian_set, symmetric_set, antisymmetric_set                           !:.:.:
 public                                 :: abs, sin, cos, exp, log, sqrt                                                            !:.:.:
 public                                 :: matrix_random_init, NaN, f_mout                                                          !:.:.:
 public                                 :: assignment(=), operator(+), operator(-), operator(*), operator(/) , operator(**)         !:.:.:
!-----------------------------------------------------------------------------------------------------------------------------------!:.:.:
!                                                                                                                                         ! Basic type with common data for all Matrix Classes:
 type               MatrixClass                                                                                                     !:.:.:
  integer                               :: m=0, n =0                                                                                !:.:.:! an m x n matrix M_ij;  m=0 indicates that it is not initiated
  integer                               :: is=0,ie=0,js=0,je=0                                                                      !:.:.:! i= is, is+1, ..., ie   j= js, js+1, ..., je, ie=is+m-1, je=js+n-1
  character(mtype_len)                  :: mtype='GG'                                                                               !:.:.:! Lapack convention: ZG, ZH, DG, DS, ... if GG: Not initialized
  character( name_len)                  :: name=''                                                                                  !:.:.:! name, further info
 contains
  procedure, private, pass(MAT)         :: matrix_save
  procedure, private, pass(MAT)         :: matrix_read 
  procedure, private, pass(MAT)         :: matrix_print
  procedure, private, pass(MAT)         :: matrix_random_set
  procedure, private, pass(MAT)         :: matrix_gaussian_set                            
  generic                               :: save           => matrix_save                                                            !:.:.:! call MAT % save (unit) makes a raw save      in (optional) unit (default f_mout)  
  generic                               :: read           => matrix_read                                                            !:.:.:! call MAT % read (unit) reads matrix        from (optional) unit (default f_minput)    
  generic                               :: print          => matrix_print                                                           !:.:.:! call MAT % print(unit) makes pretty printing in (optional) unit (default f_mout)
  generic                               :: random         => matrix_random_set                                                      !:.:.:! call MAT % random      makes matrix random.  If mtype='ZH'/'DS', it is hermitian/symmetric
  generic                               :: gaussian       => matrix_gaussian_set                                                    !:.:.:
 end type           MatrixClass
!-----------------------------------------------------------------------------------------------------------------------------------
 type,      extends(MatrixClass)        :: Matrix                                                                                   !:.:.:! Complex(dp) matrices 
  complex(dp), allocatable              :: v(:,:)                                                                                   !:.:.:! MAT % v(is:ie, js:je): values of matrices
 contains
  procedure, private, pass(MAT)         :: matrix_hermitian_set
  procedure, private, pass(MAT)         :: matrix_traceless_set
  procedure, private, pass(MATA)        :: matrix_return_transpose
  procedure, private, pass(MATA)        :: matrix_return_hermitian
  procedure, private, pass(MATA)        :: matrix_return_conjg
  procedure, private, pass(MATA)        :: matrix_return_real_dmatrix
  procedure, private, pass(MATA)        :: matrix_return_imag_dmatrix
  generic                               :: hermitian_set  => matrix_hermitian_set                                                   !:.:.:! call MAT % hermitian_set makes matrix hermitian and mtype 'ZH'
  generic                               :: traceless_set  => matrix_traceless_set
  generic                               :: conjg          => matrix_return_conjg                                                    !:.:.:! m2 = m1  % conjg    () returns the complex conjugate of matrix m1
  generic                               :: transpose      => matrix_return_transpose                                                !:.:.:! m2 = m1  % transpose() returns the transpose         of matrix m1
  generic                               :: hermitian      => matrix_return_hermitian                                                !:.:.:! m2 = m1  % hermitian() returns the hermitian         of matrix m1
  generic                               :: dagger         => matrix_return_hermitian                                                !:.:.:! m2 = m1  % hermitian() returns the hermitian         of matrix m1
  generic                               :: re             => matrix_return_real_dmatrix                                             !:.:.:
  generic                               :: im             => matrix_return_imag_dmatrix                                             !:.:.:
 end type           Matrix
 interface          Matrix                                                                                                          ! MAT=Matrix(m,n,is,js,mtype='ZG',name='MMAT')
  module procedure                      :: matrix_construct_zero  , matrix_construct_real, matrix_construct_complex                 ! MAT=Matrix(n,mtype='ZH'), MAT=Matrix(PI,n), MAT=Matrix(IMU,n)
  module procedure                      :: matrix_construct_random, matrix_construct_array2                                         ! MAT=Matrix('uniform',n,mtype='ZH'), MAT=Matrix('gaussian',n,sigma=PI), MAT=Matrix(C(:,:))
 end interface      Matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 type,      extends(MatrixClass)        :: DMatrix                                                                                  !:.:.:! Real(dp)    matrices
  real   (dp), allocatable              :: v(:,:)                                                                                   !:.:.:! MAT % v(is:ie, js:je): values of matrices
 contains
  procedure, private, pass(MAT)         :: matrix_symmetric_set_d 
  procedure, private, pass(MAT)         :: dmatrix_traceless_set 
  procedure, private, pass(MATA)        :: matrix_return_transpose_d
  generic                               :: symmetric_set  => matrix_symmetric_set_d                                                 !:.:.:! call MAT % symmetric makes matrix symmetric and mtype 'DS'
  generic                               :: traceless_set  => dmatrix_traceless_set
  generic                               :: transpose      => matrix_return_transpose_d                                              !:.:.:! m2 = m1  % transpose() returns the transpose         of matrix m1
  generic                               :: symmetric      => matrix_return_transpose_d                                              !:.:.:! m2 = m1  % symmetric() returns the symmetric         of matrix m1
 end type           DMatrix                                                                                                         
 interface          DMatrix                                                                                                         ! DMT=DMatrix(m,n,is,js,mtype='DG',name='DMAT')
  module procedure                      :: matrix_construct_zero_d  , matrix_construct_real_d  ,matrix_construct_complex_d          ! DMT=DMatrix(    n,mtype='DS'), DMT=Dmatrix(PI,n)
  module procedure                      :: matrix_construct_random_d, matrix_construct_array2_d                                     ! DMT=DMatrix('u',n,mtype='DS'), DMT=Dmatrix('g',n,sigma=PI),MAT=DMatrix(D(:,:),mtype='DS',name='D')
 end interface      DMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 type               Vector                                                                                                          !:.:.:! Complex(dp) vectors
  integer                               :: n =0                                                                                     !:.:.:! An n vector V_i
  integer                               :: is=0,ie=0                                                                                !:.:.:! i= is, is+1, ..., ie    ie=is+n-1
  character( name_len)                  :: name=''                                                                                  !:.:.:! Free use by the user
  complex  (dp), allocatable            :: v(:) 
 contains
  procedure, private, pass(vec)         :: vector_print
  procedure, private, pass(vec)         :: vector_save
  procedure, private, pass(vec)         :: vector_read
  procedure, private, pass(vec)         :: vector_random_set
  procedure, private, pass(vec)         :: vector_gaussian_set
  procedure, private, pass(vec)         :: vector_norm 
  procedure, private, pass(vec)         :: vector_return_real_dvector
  procedure, private, pass(vec)         :: vector_return_imag_dvector
  procedure, private, pass(vec)         :: vector_return_return_conjg
  generic                               :: conjg          => vector_return_return_conjg
  generic                               :: re             => vector_return_real_dvector
  generic                               :: im             => vector_return_imag_dvector
  generic                               :: print          => vector_print
  generic                               :: read           => vector_read
  generic                               :: save           => vector_save
  generic                               :: random         => vector_random_set
  generic                               :: gaussian       => vector_gaussian_set
 end type           Vector
 interface          Vector
  module procedure                      :: vector_construct_zero  ,vector_construct_real  ,vector_construct_complex
  module procedure                      :: vector_construct_random,vector_construct_array1
 end interface      Vector
!-----------------------------------------------------------------------------------------------------------------------------------
 type               DVector                                                                                                         !:.:.:! Real(dp) vectors
  integer                               :: n =0                                                                                     !:.:.:! An n vector V_i
  integer                               :: is=0,ie=0                                                                                !:.:.:! i= is, is+1, ..., ie    ie=is+n-1
  character( name_len)                  :: name=''                                                                                  !:.:.:! Free use by the user
  real     (dp), allocatable            :: v(:) 
 contains
  procedure, private, pass(vec)         :: dvector_print
  procedure, private, pass(vec)         :: dvector_save
  procedure, private, pass(vec)         :: dvector_read
  procedure, private, pass(vec)         :: dvector_random_set
  procedure, private, pass(vec)         :: dvector_gaussian_set
  generic                               :: print          => dvector_print
  generic                               :: read           => dvector_read
  generic                               :: save           => dvector_save
  generic                               :: random         => dvector_random_set
  generic                               :: gaussian       => dvector_gaussian_set
 end type           DVector
 interface          DVector
  module procedure                      :: dvector_construct_zero  ,dvector_construct_real  ,dvector_construct_complex
  module procedure                      :: dvector_construct_random,dvector_construct_array1
 end interface      DVector
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          assignment(=)                                                                                                   !:.:.:
  module procedure                      :: matrix_assignFrom_dmatrix                                                                                          
  module procedure                      :: matrix_assignFrom_real   ,matrix_assignFrom_complex  ,matrix_assignFrom_array2            
  module procedure                      :: matrix_assignFrom_real_d ,matrix_assignFrom_complex_d,matrix_assignFrom_array2_d                                                 
  module procedure                      :: vector_assignFrom_real   ,vector_assignFrom_complex  ,vector_assignFrom_array1                                                 
  module procedure                      :: vector_assignFrom_dvector                                                                          
  module procedure                      :: dvector_assignFrom_real  ,dvector_assignFrom_complex ,dvector_assignFrom_array1                                                 
  module procedure                      :: dvector_assignFrom_vector,dvector_assignFrom_array1_d
 end interface assignment(=)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          operator(+)                                                                                                     !:.:.: 
  module procedure                      :: real_plus_matrix     , matrix_plus_real    , complex_plus_matrix , matrix_plus_complex     
  module procedure                      :: matrix_plus_array2   , array2_plus_matrix  , matrix_plus_matrix                           
  module procedure                      :: matrix_plus_dmatrix  , dmatrix_plus_matrix                                                
  module procedure                      :: real_plus_matrix_d   , matrix_plus_real_d                                                 
  module procedure                      :: matrix_plus_array2_d , array2_plus_matrix_d, matrix_plus_matrix_d                                      
  module procedure                      :: real_plus_vector     , vector_plus_real    , complex_plus_vector , vector_plus_complex
  module procedure                      :: vector_plus_array1   , array1_plus_vector  , vector_plus_array1_d, array1_d_plus_vector
  module procedure                      :: vector_plus_vector   , vector_plus_dvector , dvector_plus_vector  
  module procedure                      :: real_plus_dvector    , dvector_plus_real    
  module procedure                      :: dvector_plus_array1_d, array1_d_plus_dvector
  module procedure                      :: dvector_plus_dvector 
 end interface      operator(+)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          operator(-)                                                                                                     !:.:.:                        
  module procedure                      :: real_subtract_matrix      , matrix_subtract_real
  module procedure                      :: complex_subtract_matrix   , matrix_subtract_complex
  module procedure                      :: matrix_subtract_array2    , array2_subtract_matrix       , matrix_subtract_matrix
  module procedure                      :: matrix_subtract_dmatrix   , dmatrix_subtract_matrix
  module procedure                      :: real_subtract_matrix_d    , matrix_subtract_real_d
  module procedure                      :: matrix_subtract_array2_d  , array2_subtract_matrix_d     , matrix_subtract_matrix_d
  module procedure                      :: matrix_return_minus_matrix, matrix_return_minus_matrix_d
  module procedure                      :: real_subtract_vector      , vector_subtract_real
  module procedure                      :: complex_subtract_vector   , vector_subtract_complex
  module procedure                      :: vector_subtract_array1    , array1_subtract_vector
  module procedure                      :: vector_subtract_array1_d  , array1_d_subtract_vector
  module procedure                      :: dvector_subtract_vector   , vector_subtract_dvector      , vector_subtract_vector
  module procedure                      :: vector_return_minus_vector
  module procedure                      :: real_subtract_dvector     , dvector_subtract_real  
  module procedure                      :: dvector_subtract_array1_d , array1_d_subtract_dvector
  module procedure                      :: dvector_subtract_dvector  , dvector_return_minus_dvector  
 end interface      operator(-)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          operator(*)                                                                                                     !:.:.: 
  module procedure                      :: matrix_mult_real    , real_mult_matrix    , matrix_mult_complex  , complex_mult_matrix
  module procedure                      :: matrix_mult_array2  , array2_mult_matrix  , matrix_mult_matrix
  module procedure                      :: matrix_mult_dmatrix , dmatrix_mult_matrix
  module procedure                      :: matrix_mult_real_d  , real_mult_matrix_d  , matrix_mult_complex_d, complex_mult_matrix_d
  module procedure                      :: matrix_mult_array2_d, array2_mult_matrix_d, matrix_mult_matrix_d
  module procedure                      :: real_mult_vector    , vector_mult_real    , complex_mult_vector  , vector_mult_complex
  module procedure                      :: real_mult_dvector   , dvector_mult_real    
  module procedure                      :: matrix_mult_vector  , vector_mult_matrix
  module procedure                      :: dmatrix_mult_dvector, dvector_mult_dmatrix
 end interface      operator(*)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          operator(/)                                                                                                     !:.:.:                                                                                          
  module procedure                      :: matrix_divide_real   , matrix_divide_complex
  module procedure                      :: matrix_divide_real_d 
  module procedure                      :: vector_divide_real   , vector_divide_complex
  module procedure                      :: dvector_divide_real
 end interface      operator(/)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          operator(**)                                                                                                    !:.:.: 
  module procedure                      ::  matrix_power_integer,  matrix_power_real, matrix_power_complex
  module procedure                      :: dmatrix_power_integer, dmatrix_power_real
  module procedure                      ::  vector_power_integer,  vector_power_real, vector_power_complex
  module procedure                      :: dvector_power_integer, dvector_power_real
 end interface      operator(**)
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          random_number                                                                                                   !:.:.:
  module procedure                      :: random_number_matrix ,random_number_matrix_gaussian
  module procedure                      :: random_number_dmatrix,random_number_dmatrix_gaussian
  module procedure                      :: vector_random_set    ,random_number_vector_gaussian_set
  module procedure                      :: dvector_random_set   ,random_number_dvector_gaussian_set
 end interface      random_number
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          mmmult                                                                                                          !:.:.:
  module procedure                      :: matrix_mult_matrix_sub, matrix_mult_matrix_sub_d
 end interface      mmmult
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          mvmult                                                                                                          !:.:.:
  module procedure                      :: matrix_mult_vector_sub, matrix_mult_vector_sub_d  
 end interface      mvmult
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          vmmult                                                                                                          !:.:.:
  module procedure                      :: vector_mult_matrix_sub, vector_mult_matrix_sub_d 
 end interface      vmmult
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          conjg                                                                                                           !:.:.:
  module procedure                      :: matrix_return_conjg
 end interface      conjg
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          transpose                                                                                                       !:.:.:
  module procedure                      :: matrix_return_transpose
  module procedure                      :: matrix_return_transpose_d
 end interface      transpose
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          hermitian                                                                                                       !:.:.:
  module procedure                      :: matrix_return_hermitian
 end interface      hermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          symmetric                                                                                                       !:.:.:
  module procedure                      :: matrix_return_transpose
  module procedure                      :: matrix_return_transpose_d       
 end interface      symmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          trace                                                                                                           !:.:.:
  module procedure                      :: matrix_trace  , dmatrix_trace
 end interface      trace
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          trace2                                                                                                          !:.:.:
  module procedure                      :: matrix_trace2 , dmatrix_trace2
 end interface      trace2
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          trace2c                                                                                                         !:.:.:
  module procedure                      :: matrix_trace2c, dmatrix_trace2c
 end interface      trace2c
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          inverse                                                                                                         !:.:.:
  module procedure                      :: matrix_inverse      , dmatrix_inverse
 end interface      inverse
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          eigenvalues                                                                                                     !:.:.:
  module procedure                      :: matrix_eigenvalues  , dmatrix_eigenvalues
 end interface      eigenvalues
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          eigenvectors                                                                                                    !:.:.:
  module procedure                      :: matrix_eigenvectors , dmatrix_eigenvectors
 end interface      eigenvectors
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          determinant                                                                                                     !:.:.:
  module procedure                      :: matrix_determinant  , dmatrix_determinant
 end interface      determinant
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          lndet                                                                                                           !:.:.:
  module procedure                      :: matrix_lndet        , dmatrix_lndet
 end interface      lndet
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          Pfaffian                                                                                                        !:.:.:
  module procedure                      :: matrix_Pfaffian
 end interface      Pfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          lnPfaffian                                                                                                      !:.:.:
  module procedure                      :: matrix_lnPfaffian
 end interface      lnPfaffian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          diagonal                                                                                                        !:.:.:
  module procedure                      ::  matrix_diagonal_get ,  matrix_diagonal_set_from_vector 
  module procedure                      :: dmatrix_diagonal_get , dmatrix_diagonal_set_from_dvector
 end interface      diagonal
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          diagonalMatrix                                                                                                  !:.:.:
  module procedure                      ::  matrix_diagonal_set_from_matrix  , dmatrix_diagonal_set_from_dmatrix
 end interface      diagonalMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          sort                                                                                                            !:.:.:
  module procedure                      :: vector_sort, dvector_sort
 end interface      sort
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          norm                                                                                                            !:.:.:
  module procedure                      :: matrix_norm,matrix_norm_d, vector_norm, dvector_norm
 end interface      norm
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          real                                                                                                            !:.:.:
  module procedure                      :: matrix_return_real_dmatrix
 end interface      real
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          aimag                                                                                                           !:.:.:
  module procedure                      :: matrix_return_imag_dmatrix
 end interface      aimag
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          mcmplx                                                                                                          !:.:.:
  module procedure                      :: dmatrix_dmatrix_complex_return_matrix, dvector_dvector_complex_return_vector
 end interface      mcmplx
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          dot_product                                                                                                     !:.:.:
  module procedure                      :: vector_dot_product, dvector_dot_product
 end interface      dot_product
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          maxval                                                                                                          !:.:.:
  module procedure                      :: dvector_maxval
 end interface      maxval
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          minval                                                                                                          !:.:.:
  module procedure                      :: dvector_minval
 end interface      minval
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          traceless                                                                                                       !:.:.:
  module procedure                      :: matrix_traceless_get, dmatrix_traceless_get
 end interface      traceless
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          traceless_set                                                                                                   !:.:.:
  module procedure                      :: matrix_traceless_set, dmatrix_traceless_set
 end interface      traceless_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          hermitian_set                                                                                                   !:.:.:
  module procedure                      :: matrix_hermitian_set
 end interface      hermitian_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          symmetric_set                                                                                                   !:.:.:
  module procedure                      :: matrix_symmetric_set, matrix_symmetric_set_d
 end interface      symmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          antisymmetric_set                                                                                               !:.:.:
  module procedure                      :: matrix_antisymmetric_set, dmatrix_antisymmetric_set
 end interface      antisymmetric_set
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isHermitian                                                                                                     !:.:.:
  module procedure                      :: matrix_is_hermitian
 end interface      isHermitian
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isSymmetric                                                                                                     !:.:.:
  module procedure                      :: matrix_is_symmetric    , dmatrix_is_symmetric
 end interface      isSymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isAntiSymmetric                                                                                                 !:.:.:
  module procedure                      :: matrix_is_antisymmetric, dmatrix_is_antisymmetric
 end interface      isAntiSymmetric
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          isNaN                                                                                                           !:.:.:
  module procedure                      :: vector_is_nan, dvector_is_nan, matrix_is_nan, dmatrix_is_nan
 end interface      isNaN
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          metadata_copy                                                                                                   !:.:.:
  module procedure                      :: vector_metadata_copy_dvector, dvector_metadata_copy_vector, dvector_metadata_copy_dvector
  module procedure                      :: matrix_metadata_copy        , vector_metadata_copy_vector
 end interface      metadata_copy
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          abs                                                                                                             !:.:.:
  module procedure                      :: matrix_abs  , dmatrix_abs   , vector_abs  , dvector_abs
 end interface      abs
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          sin                                                                                                             !:.:.:
  module procedure                      :: matrix_sin  , dmatrix_sin   , vector_sin  , dvector_sin
 end interface      sin
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          cos                                                                                                             !:.:.:
  module procedure                      :: matrix_cos  , dmatrix_cos   , vector_cos  , dvector_cos
 end interface      cos
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          exp                                                                                                             !:.:.:
  module procedure                      :: matrix_exp  , dmatrix_exp   , vector_exp  , dvector_exp
 end interface      exp
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          log                                                                                                             !:.:.:
  module procedure                      :: matrix_log  , dmatrix_log   , vector_log  , dvector_log
 end interface      log
!-----------------------------------------------------------------------------------------------------------------------------------
 interface          sqrt                                                                                                            !:.:.:
  module procedure                      :: matrix_sqrt , dmatrix_sqrt  , vector_sqrt , dvector_sqrt
 end interface      sqrt
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
 include 'matrix_mod_matrix_matrixClass.f90'
 include 'matrix_mod_matrix_matrix.f90'
 include 'matrix_mod_matrix_dmatrix.f90'
 include 'matrix_mod_matrix_vector.f90'
 include 'matrix_mod_matrix_dvector.f90'
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module          matrix_mod_matrix
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
