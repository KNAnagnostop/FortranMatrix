!...................................................................................................................................!:.:.:
!.......................... File: timer_mod.f90        .............................................................................!:.:.:
!
! see https://fortranwiki.org/fortran/show/Timing for more details and a more involved timing module
!-----------------------------------------------------------------------------------------------------------------------------------
! Usage:
! use  timer_mod
! call timer_mod_init    ; call timer_mod_init (un)                                                     ! prints output message in unit un
! call timer_mod_print   ; call timer_mod_print(un,messg);call timer_mod_print(messg='Timing message'); 
!
! Control output:
! call timer_mod_now
! print *, timer_count,timer_cputot,timer_cpudel,timer_wtimetot,timer_wtimedel,trim(timer_date)
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
module      timer_mod
 use, intrinsic :: iso_fortran_env
 implicit none
 private
!-----------------------------------------------------------------------------------------------------------------------------------
 public              :: timer_mod_init, timer_mod_now,timer_mod_print
 public              :: timer_cputot, timer_cpudel, timer_count, timer_date, timer_wtimetot, timer_wtimedel
!-----------------------------------------------------------------------------------------------------------------------------------
 integer , parameter :: dp = real64, ip= int64
!-----------------------------------------------------------------------------------------------------------------------------------
 real     (dp)       :: timer_cputot, timer_cpudel, timer_wtimetot, timer_wtimedel
 integer             :: timer_count
!-----------------------------------------------------------------------------------------------------------------------------------
 real     (dp)       :: cpu_start  , cpu_mark  , cpu_now              ! _mark variables store time since last call to timer_mod_now
 integer  (ip)       :: wtime_start, wtime_mark, wtime_now, wtime_rate
 character(100)      :: timer_date
!-----------------------------------------------------------------------------------------------------------------------------------
 real(dp), parameter :: ZERO=0.0_dp


!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine     timer_mod_print(un,messg)                                                                                           !:.:.:
  integer     , intent(in), optional :: un
  character(*), intent(in), optional :: messg
  integer                            :: un_
  character(:), allocatable          :: messg_

  if(present(un))then
   un_ = un
  else
   un_ = output_unit
  end if

  if(present(messg))then
   messg_ = messg
  else
   messg_ = ' '
  end if

  call timer_mod_now

  write(un_,'(A,I6,A,4G28.16,A)') '#TIMERMOD ', &
       timer_count,' '//trim(messg_)//' ',timer_cputot,timer_cpudel,timer_wtimetot,timer_wtimedel,trim(timer_date)

 end subroutine timer_mod_print
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine     timer_mod_now                                                                                                       !:.:.:
 
  timer_count   = timer_count + 1

  cpu_mark      = cpu_now
  call cpu_time(cpu_now)

  timer_cputot  = cpu_now - cpu_start
  timer_cpudel  = cpu_now - cpu_mark

  wtime_mark    = wtime_now   
  call system_clock(wtime_now,wtime_rate)
  timer_wtimetot = real(wtime_now-wtime_start,kind=dp)/real(wtime_rate,kind=dp)
  timer_wtimedel = real(wtime_now-wtime_mark ,kind=dp)/real(wtime_rate,kind=dp)

  call timer_mod_date

 end subroutine timer_mod_now
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine     timer_mod_init(un)                                                                                                  !:.:.:
  integer, intent(in), optional :: un
  integer                       :: un_

  if(present(un))then
   un_ = un
  else
   un_ = output_unit
  end if

  call cpu_time(cpu_start)

  timer_count  = 0
  cpu_mark     = cpu_start
  cpu_now      = cpu_start

  call system_clock(wtime_start,wtime_rate)
  wtime_mark  = wtime_start
  wtime_now   = wtime_start
  
  call timer_mod_date

  write(un_,'(A)') &
       '#TIMERINIT:: ' // trim(timer_date)
  write(un_,'(A)') &
       '#TIMERINIT::   '// &
       'i         cpu_tot                     cpu_del                     wtime_tot                   wtime_del            Date'
  
 end subroutine timer_mod_init
!-----------------------------------------------------------------------------------------------------------------------------------
 subroutine     timer_mod_date                                                                                                      !:.:.:
  integer d(8)

  call date_and_time(values=d)

  write(timer_date,'(A,I0.4,A1,I0.2,A1,I0.2,A,I0.2,A,I0.2,A,I0.2,A,I4)') &
       'Date: ',d(1),'-',d(2),'-',d(3),'    ',d(5),':',d(6),':',d(7),'   ms: ',d(8)

 end subroutine timer_mod_date
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
end module timer_mod
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
!  Copyright by Konstantinos N. Anagnostopoulos, Physics Department, National Technical University of Athens, 2025
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
!program     testme
! use timer_mod
! implicit none
! integer :: N, i
! real(8), allocatable :: A(:,:), B(:,:)
!
! N=4000
! allocate(A(N,N)); A = 0.999_8
! call timer_mod_init; 
! print *,'init on: ',trim(timer_date)
! B = matmul(A,A)
! call timer_mod_now
! print *,'timer1: ',timer_count,timer_cputot,timer_cpudel,timer_wtimetot,timer_wtimedel,trim(timer_date)
! B = matmul(A,A)
! call timer_mod_print(messg='End of program:')
!end program testme

