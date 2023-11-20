!=================================================================================!
!                                     RANDOM                                      !
!=================================================================================!
!                                                                                 !
! This module is part of SHEAP --- for generating random points                   !
! Copyright (C) 2019-2022 Ben Shires, Chris Pickard                               !
!                                                                                 !
! This program is free software; you can redistribute it and/or                   !
! modify it under the terms of the GNU General Public License                     !
! as published by the Free Software Foundation; either version 2                  !
! of the License, or (at your option) any later version.                          !
!                                                                                 !
! This program is distributed in the hope that it will be useful,                 !
! but WITHOUT ANY WARRANTY; without even the implied warranty of                  !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   !
! GNU General Public License for more details.                                    !
!                                                                                 !
! You should have received a copy of the GNU General Public License               !
! along with this program; if not, write to the Free Software                     !
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. !
!                                                                                 !
!---------------------------------------------------------------------------------!

module random

  ! * Prerequisite modules
  use constants
  use rng

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public subroutines
  public :: random_add_noise
  public :: random_point_in_hypersphere
  
contains

  subroutine random_add_noise(vec,sig)

    integer                                    :: n

    real(kind=dp)                              :: U,V
    real(kind=dp), intent(in)                  :: sig
    real(kind=dp), dimension(:), intent(inout) :: vec

    do n=1,size(vec)

       U=random_single()
       V=random_single()

       vec(n)=vec(n)+sqrt(-2*log(U))*cos(tpi*V)*sig

    end do

  end subroutine random_add_noise

  subroutine random_point_in_hypersphere(vec,r) 

    integer                                  :: n

    real(kind=dp)                            :: d
    real(kind=dp),               intent(in)  :: r
    real(kind=dp), dimension(:), intent(out) :: vec

    n=size(vec)

    vec=0.0_dp

    d=r*random_single()**(1/real(n,dp))

    call random_add_noise(vec,1.0_dp)

    vec=d*vec/norm2(vec)

  end subroutine random_point_in_hypersphere
  
end module random
