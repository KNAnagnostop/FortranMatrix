program             testme
 use                matrix_mod_common
 use                matrix_mod_basic
 use                matrix_mod_array
 implicit none
 integer,parameter :: f_out = 13
 type( Matrix)     :: mm
 type(DMatrix)     :: dd
 complex(dp)       :: aa(2,2)
 real   (dp)       :: bb(2,2)

 open(unit=f_out,file='test.dat')
 f_mout   = f_out
 call matrix_random_init

 write(f_out,*) '------------------------1-----------------------------'
 mm = Matrix((2.3_dp,-1.1_dp),3,3,0,-1,mtype='ZG',name='MMAT')
 write(f_out,*) mm % m, mm % n, mm % is, mm % ie, mm % js, mm % je
 write(f_out,*) mm % mtype,trim(mm % name)
 call mm % print(f_out,fmt='F6.2')
 call mm % hermitian
 call mm % print(f_out,fmt='F6.2')
 write(f_out,*) '------------------------2-----------------------------'
 dd = DMatrix(7.3_dp,3,3,0,-1,mtype='DG',name='DMAT')
 write(f_out,*) dd % m, dd % n, dd % is, dd % ie, dd % js, dd % je
 write(f_out,*) dd % mtype,trim(dd % name)
 call dd % print(f_out,fmt='F6.2')
 write(f_out,*) '------------------------3-----------------------------'
 mm % mtype = 'ZH'
 call mm % gaussian
 call mm % print(f_out,fmt='F8.3')
 write(f_out,*) '------------------------4-----------------------------'
 dd % mtype = 'DS'
 call dd % gaussian(sigma=PI)
 call dd % print(f_out,fmt='F8.3')
 write(f_out,*) '------------------------5-----------------------------'
 dd = DMatrix('g',4,mtype='DG',sigma=PI)
 call dd % print(f_out,fmt='F8.3')
 write(f_out,*) '------------------------6-----------------------------'
 mm =  Matrix('g',4,3,js=3,mtype='ZG',name='MMAT')
 call mm % print(f_out,fmt='F8.3')
 write(f_out,*) '------------------------7-----------------------------'
 call array2_random_set(aa)
 print *, aa
end program         testme
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
