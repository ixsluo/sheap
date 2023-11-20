!=================================================================================!
!                                    CONSTANTS                                    !
!=================================================================================!
!                                                                                 !
! This module is part of SHEAP - contains various constants                       !
! Copyright (C) 2019-2021  Ben Shires, Chris Pickard                              !
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

module constants

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public variables
  integer, public, parameter       :: dp=8 ! Single = 4 Double = 8 Quad = 16
  integer, public, parameter       :: stdin=5
  integer, public, parameter       :: stdout=6
  integer, public, parameter       :: stderr=0

  real(kind=dp), public, parameter :: pi=3.141592653589793238462643383279502884197_dp
  real(kind=dp), public, parameter :: tpi=2.0_dp*pi
  real(kind=dp), public, parameter :: gr=(sqrt(5.0_dp)+1.0_dp)/2.0_dp
  real(kind=dp), public, parameter :: dgrd = pi/180.0_dp
  real(kind=dp), public, parameter :: evbyang3=160.2176487_dp
  real(kind=dp), public, parameter :: bohr2ang = 0.529177210903_dp
  real(kind=dp), public, parameter :: delta = 1e-13_dp

  real(kind=dp), public, parameter, dimension(3,3) :: &
       ident=reshape((/1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp/),(/3,3/))

end module constants
