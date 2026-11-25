program testme
 integer, parameter   :: dp=8
 real(8), allocatable :: A(:,:), B(:,:)
 real(8)              :: C(4,3)
 
 allocate(A(5:9,4:7))
 A = 3.0_dp
 C = 2.0_dp

 allocate(B,source=A) ! same effect as B = A
!B =  A
 print *, size(B),lbound(B),ubound(B)
 B =  A+C

 do i = lbound(B,1), ubound(B,1)
  print '(2000F8.3)',( B(i,j), j=lbound(B,2),ubound(B,2))
 end do
 print *, size(B),lbound(B),ubound(B)

 print '(A,2000F8.3)','B=', B

end program testme

!-------------------------------------------------------------------
!!!! BUG?? !!!!
!gfortran hi.f90 -o hi;./hi
!          20           5           4           9           7
!   5.000   5.000   5.000
!   5.000   5.000   5.000
!   5.000   5.000   5.000
!   5.000   5.000   5.000
!          12           1           1           4           3
! B=   5.000   5.000   5.000   5.000   5.000   5.000   5.000   5.000   5.000   5.000   5.000   5.000
!-------------------------------------------------------------------
! ifort hi.f90 -o hi;./hi
!          20           5           4           9           7
!   5.000   5.000   5.000   3.000
!   5.000   5.000   5.000   3.000
!   5.000   5.000   5.000   3.000
!   5.000   5.000   5.000   3.000
!   3.000   3.000   3.000   3.000
!          20           5           4           9           7
! B=   5.000   5.000   5.000   5.000   3.000   5.000   5.000   5.000   5.000   3.000   5.000   5.000   5.000   5.000   3.000   3.000   3.000   3.000   3.000   3.000




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
