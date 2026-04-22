!...................................................................................................................................!:.:.:
!.......................... File: tensorprod_mod.f90   .............................................................................!:.:.:
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
module              tensorprod_mod                                                                                                  !:.:.:
 use                matrix_mod_common
 use                matrix_mod_array
 implicit none
 private ; save

!-----------------------------------------------------------------------------------------------------------------------------------!:.:.:
 public                                 :: tensorprod, vec2row, vec2col, vectorize                                                  !:.:.:
!-----------------------------------------------------------------------------------------------------------------------------------!:.:.:


!-----------------------------------------------------------------------------------------------------------------------------------
 interface          tensorprod                                                                                                      !:.:.:
  module procedure                      :: tensorprod_complex_complex
  module procedure                      :: tensorprod_complex_3,tensorprod_complex_4,tensorprod_complex_5,tensorprod_complex_6
  module procedure                      :: tensorprod_real_real
  module procedure                      :: tensorprod_real_3   ,tensorprod_real_4   ,tensorprod_real_5   ,tensorprod_real_6
  module procedure                      :: tensorprod_complex_vec_vec, tensorprod_real_vec_vec
 end interface      tensorprod

 interface          vec2col
  module procedure                      :: vec2col_complex, vec2col_real
 end interface      vec2col

 interface          vec2row
  module procedure                      :: vec2row_complex, vec2row_real 
 end interface      vec2row

 interface          vectorize
  module procedure                      :: vectorize_complex ,  vectorize_real
 end interface      vectorize

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

! Tensor product: Given X(Nm,Nn) and Y(Np,Nq), then Z = X  x Y, where Z(Ni,Nj), Ni=Nm*Np, Nj=Nn*Nq, such that:
!
! X(m,n) m=1,...,Nm  n=1,...,Nn
! Y(p,q) p=1,...,Np  q=1,...,Nq
! Z(i,j) i=1,...,Ni  j=1,...,Nj    Ni=Nm*Np,     Nj=Nn*Nq
!
! Z(i,j) = X( (i/Np) + 1, (j/Nq) + 1 )   *   Y( (i-1) % Np + 1 , (j-1) % Nq  + 1 )
!
!        = X(m,n) * Y(p,q)  with i=Np*(m-1)+p   j=Nq*(n-1)+q
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of two complex matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_complex_complex                             (X,Y)                    result(Z)                       !:.:.:
!     function      tensorprod_complex_complex                             (X,Y)                    result(Z)                       
  complex(dp), intent(in)               ::         X(:,:), Y(:,:)
  complex(dp), allocatable              :: Z(:,:)
  integer                               :: Ni, Nj, Nm, Nn, Np, Nq
  integer                               ::  i,  j,  m,  n,  p,  q 

  Nm = size(X,1) ; Nn = size(X,2)
  Np = size(Y,1) ; Nq = size(Y,2)
  Ni = Nm * Np   ; Nj = Nn * Nq

  allocate(Z(Ni,Nj))

  j           = 0

  do    n     = 1, Nn
   do   q     = 1, Nq


    j         = j + 1                                   !; j         = Nq * (n-1) + q    !;     print *,'XX: j=',j, Nq * (n-1) + q - j

    i         = 0

    do  m     = 1, Nm
     do p     = 1, Np


      i       = i + 1                                   !;  i        = Np * (m-1) + p    !;     print *,'YY: i=',i, Np * (m-1) + p - i

      Z(i,j)  = X(m,n) * Y(p,q)

     end do
    end do
   end do
  end do

 end  function      tensorprod_complex_complex
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 3 complex matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_complex_3                                   (X1,X2,X3)               result(Z)                       !:.:.:
  complex(dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:)
  complex(dp), allocatable              :: XX(:,:)
  complex(dp), allocatable              :: Z (:,:)

  XX = tensorprod_complex_complex(X2,X3)
  Z  = tensorprod_complex_complex(X1,XX)

 end function       tensorprod_complex_3
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 4 complex matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_complex_4                                   (X1,X2,X3,X4)            result(Z)                       !:.:.:
  complex(dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:), X4(:,:)
  complex(dp), allocatable              :: XX(:,:)
  complex(dp), allocatable              :: Z (:,:)

  XX = tensorprod_complex_complex(X3,X4)
  XX = tensorprod_complex_complex(X2,XX)
  Z  = tensorprod_complex_complex(X1,XX)

 end function       tensorprod_complex_4
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 5 complex matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_complex_5                                   (X1,X2,X3,X4,X5)         result(Z)                       !:.:.:
  complex(dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:), X4(:,:), X5(:,:)
  complex(dp), allocatable              :: XX(:,:)
  complex(dp), allocatable              :: Z (:,:)

  XX = tensorprod_complex_complex(X4,X5)
  XX = tensorprod_complex_complex(X3,XX)
  XX = tensorprod_complex_complex(X2,XX)
  Z  = tensorprod_complex_complex(X1,XX)

 end function       tensorprod_complex_5
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 6 complex matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_complex_6                                   (X1,X2,X3,X4,X5,X6)      result(Z)                       !:.:.:
  complex(dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:), X4(:,:), X5(:,:), X6(:,:)
  complex(dp), allocatable              :: XX(:,:)
  complex(dp), allocatable              :: Z (:,:)

  XX = tensorprod_complex_complex(X5,X6)
  XX = tensorprod_complex_complex(X4,XX)
  XX = tensorprod_complex_complex(X3,XX)
  XX = tensorprod_complex_complex(X2,XX)
  Z  = tensorprod_complex_complex(X1,XX)

 end function       tensorprod_complex_6
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of two real matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_real_real                                   (X,Y)                    result(Z)                       !:.:.:
  real   (dp), intent(in)               ::         X(:,:), Y(:,:)
  real   (dp), allocatable              :: Z(:,:)
  integer                               :: Ni, Nj, Nm, Nn, Np, Nq
  integer                               ::  i,  j,  m,  n,  p,  q 

  Nm = size(X,1) ; Nn = size(X,2)
  Np = size(Y,1) ; Nq = size(Y,2)
  Ni = Nm * Np   ; Nj = Nn * Nq

  allocate(Z(Ni,Nj))

  j           = 0

  do    n     = 1, Nn
   do   q     = 1, Nq


    j         = j + 1                                   !;  j       = Nq * (n-1) + q

    i         = 0

    do  m     = 1, Nm
     do p     = 1, Np


      i       = i + 1                                   !;  i       = Np * (m-1) + p 

      Z(i,j)  = X(m,n) * Y(p,q)

     end do
    end do
   end do
  end do

 end  function      tensorprod_real_real
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 3 real matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_real_3                                      (X1,X2,X3)               result(Z)                       !:.:.:
  real   (dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:)
  real   (dp), allocatable              :: XX(:,:)
  real   (dp), allocatable              :: Z (:,:)

  XX = tensorprod_real_real(X2,X3)
  Z  = tensorprod_real_real(X1,XX)

 end function       tensorprod_real_3
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 4 real matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_real_4                                      (X1,X2,X3,X4)            result(Z)                       !:.:.:
  real   (dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:), X4(:,:)
  real   (dp), allocatable              :: XX(:,:)
  real   (dp), allocatable              :: Z (:,:)

  XX = tensorprod_real_real(X3,X4)
  XX = tensorprod_real_real(X2,XX)
  Z  = tensorprod_real_real(X1,XX)

 end function       tensorprod_real_4
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 5 real matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_real_5                                      (X1,X2,X3,X4,X5)         result(Z)                       !:.:.:
  real   (dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:), X4(:,:), X5(:,:)
  real   (dp), allocatable              :: XX(:,:)
  real   (dp), allocatable              :: Z (:,:)

  XX = tensorprod_real_real(X4,X5)
  XX = tensorprod_real_real(X3,XX)
  XX = tensorprod_real_real(X2,XX)
  Z  = tensorprod_real_real(X1,XX)

 end function       tensorprod_real_5
!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of 6 real matrices:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_real_6                                      (X1,X2,X3,X4,X5,X6)      result(Z)                       !:.:.:
  real   (dp), intent(in)               :: X1(:,:), X2(:,:), X3(:,:), X4(:,:), X5(:,:), X6(:,:)
  real   (dp), allocatable              :: XX(:,:)
  real   (dp), allocatable              :: Z (:,:)

  XX = tensorprod_real_real(X5,X6)
  XX = tensorprod_real_real(X4,XX)
  XX = tensorprod_real_real(X3,XX)
  XX = tensorprod_real_real(X2,XX)
  Z  = tensorprod_real_real(X1,XX)

 end function       tensorprod_real_6
!-----------------------------------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of two complex vectors: u is considered 1 x Nn matrix (=uT) and v a Np x 1 matrix. 
! returns: Z = uT \otimes v,  a Np x Nn matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_complex_vec_vec                             (u,v)                    result(Z)                       !:.:.:
  complex(dp), intent(in)               ::         u(:), v(:)
  complex(dp), allocatable              :: Z(:,:)
  integer                               ::         Nn  , Np
  integer                               ::          n  ,  p

  Nn = size(u) ; Np = size(v)
  allocate(Z(Np,Nn))

  do  n = 1, Nn
   do p = 1, Np

    Z(p,n) = u(n) * v(p)

   end do
  end do

 end function       tensorprod_complex_vec_vec

!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of two real vectors: u is considered 1 x Nn matrix (=uT) and v a Np x 1 matrix. 
! returns: Z = uT \otimes v,  a Np x Nn matrix
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      tensorprod_real_vec_vec                                (u,v)                    result(Z)                       !:.:.:
  real   (dp), intent(in)               ::         u(:), v(:)
  real   (dp), allocatable              :: Z(:,:)
  integer                               ::         Nn  , Np
  integer                               ::          n  ,  p

  Nn = size(u) ; Np = size(v)
  allocate(Z(Np,Nn))

  do  n = 1, Nn
   do p = 1, Np

    Z(p,n) = u(n) * v(p)

   end do
  end do

 end function       tensorprod_real_vec_vec


!-----------------------------------------------------------------------------------------------------------------------------------
! Tensor product of vector and matrices can be ambiguous in other cases. The solution is to convert a vector to an array with two indices
! and then use the tensor product at will.
!-----------------------------------------------------------------------------------------------------------------------------------
! a vector of size Np is returned as a Np x 1 array:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vec2col_complex                                        (u)                      result(Z)                       !:.:.:
  complex(dp), intent(in)               :: u(:)
  complex(dp), allocatable              :: Z(:,:)
  integer                               :: Np, p

  Np     = size(u)
  allocate(Z(Np,1))

  Z(:,1) = u
  
 end function       vec2col_complex
!-----------------------------------------------------------------------------------------------------------------------------------
! a vector of size Nn is returned as a 1 x Nn array:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vec2row_complex                                        (u)                      result(Z)                       !:.:.:
  complex(dp), intent(in)               :: u(:)
  complex(dp), allocatable              :: Z(:,:)
  integer                               :: Nn, n

  Nn     = size(u)
  allocate(Z(1,Nn))

  Z(1,:) = u
  
 end function       vec2row_complex
!-----------------------------------------------------------------------------------------------------------------------------------
! a vector of size Np is returned as a Np x 1 array:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vec2col_real                                           (u)                      result(Z)                       !:.:.:
  real   (dp), intent(in)               :: u(:)
  real   (dp), allocatable              :: Z(:,:)
  integer                               :: Np, p

  Np     = size(u)
  allocate(Z(Np,1))

  Z(:,1) = u
  
 end function       vec2col_real
!-----------------------------------------------------------------------------------------------------------------------------------
! a vector of size Nn is returned as a 1 x Nn array:
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vec2row_real                                           (u)                      result(Z)                       !:.:.:
  real   (dp), intent(in)               :: u(:)
  real   (dp), allocatable              :: Z(:,:)
  integer                               :: Nn, n

  Nn     = size(u)
  allocate(Z(1,Nn))

  Z(1,:) = u
  
 end function       vec2row_real

!-----------------------------------------------------------------------------------------------------------------------------------
! vec(Z) is a column vector that is constructed by putting the columns of Z the one below the other:
!
! Z(Nm,Nn) -> u(Nm * Nn)
!
! u(i,j) = Z(m,n)   with  i = (n-1) * Nm + m ,       m = (i-1) % Nm + 1       n = (i/Nm) + 1
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vectorize_complex                                      (Z)                      result(u)                       !:.:.:
  complex(dp), intent(in)               :: Z(:,:)
  complex(dp), allocatable              :: u(:)
  integer                               :: Nm, Nn, Ni
  integer                               ::  m,  n,  i

  Nm = size(Z,1) ; Nn = size(Z,2); Ni = Nm * Nn
  allocate(u(Ni))

  i       = 0

  do   n  = 1, Nn
   do  m  = 1, Nm

    i     = i + 1
    u(i)  = Z(m,n)

   end do
  end do

 end function       vectorize_complex
!-----------------------------------------------------------------------------------------------------------------------------------
 pure function      vectorize_real                                         (Z)                      result(u)                       !:.:.:
  real   (dp), intent(in)               :: Z(:,:)
  real   (dp), allocatable              :: u(:  )
  integer                               :: Nm, Nn, Ni
  integer                               ::  m,  n,  i

  Nm = size(Z,1) ; Nn = size(Z,2); Ni = Nm * Nn
  allocate(u(Ni))

  i       = 0

  do   n  = 1, Nn
   do  m  = 1, Nm

    i     = i + 1
    u(i)  = Z(m,n)

   end do
  end do

 end function       vectorize_real
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module          tensorprod_mod
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================


!-----------------------------------------------------------------------------------------------------------------------------------
!  Copyright by Konstantinos N. Anagnostopoulos, Physics Department, National Technical University of Athens, 2022
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
