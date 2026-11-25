!
! Verification test for BUG-07 (unallocated allocatable array check in print -> expected matrix_error abort)
! gfortran -fopenmp -I../lib/gnu/ test_ag06.f90 -L../lib/gnu -lMatrix -llapack -lblas -o t; ./t
!
program test_bug07_unalloc
  use array_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  complex(dp), allocatable :: c(:,:)
  call print(c)
end program test_bug07_unalloc
