!...................................................................................................................................!:.:.:
!.......................... File: array_mod.f90        .............................................................................!:.:.:
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
module              array_mod                                                                                                       !:.:.:
 use                matrix_mod_common
 use                matrix_mod_array
 use                tensorprod_mod
 implicit none
 private;save
!-----------------------------------------------------------------------------------------------------------------------------------!:.:.:
 public                                 :: mmmult , mvmult, vmmult                                                                  !:.:.:
 public                                 :: lmatmul, eigenvalues, eigenvectors, determinant, lndet, pfaffian, lnPfaffian, inverse    !:.:.:
 public                                 :: hermitian, hermitian_set, symmetric, symmetric_set, antisymmetric, antisymmetric_set     !:.:.:
 public                                 :: diagonal, diagonalMatrix, trace, trace2, trace2c, traceless, traceless_set               !:.:.:
 public                                 :: identitymatrix, cidentitymatrix, didentitymatrix, paulimatrix                            !:.:.:
 public                                 :: norm, isHermitian, isSymmetric, isAntisymmetric, sort                                    !:.:.:
 public                                 :: random_number, matrix_random_init                                                        !:.:.:
 public                                 :: print, printna, save, read, isNaN, NaN                                                   !:.:.:
 public                                 :: operator(.mm.)                                                                           !:.:.:
 public                                 :: tensorprod, vec2row, vec2col, vectorize                                                  !:.:.:
!-----------------------------------------------------------------------------------------------------------------------------------!:.:.:
 public                                 :: f_mout, f_minput                                                                         !:.:.:


!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module          array_mod
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
