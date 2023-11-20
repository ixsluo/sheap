!=================================================================================!
!                                      LOSS                                       !
!=================================================================================!
!                                                                                 !
! This module is part of SHEAP --- routines for computing loss function gradient  !
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

module loss

  ! * Prerequisite modules
  use constants
  use control
  use model
  use rng

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public subroutines
  public :: loss_gradient

  ! * Public variables
  logical, public :: core =.false.

  ! * Private variables
  integer                                            :: n,ni,nj,neg

  real(kind=dp)                                      :: qij,Z,coeff,temp,dist
  real(kind=dp), save, allocatable, dimension(:)     :: vec
  real(kind=dp), save, allocatable, dimension(:,:,:) :: vecs
  real(kind=dp), save, allocatable, dimension(:,:)   :: rij2

contains

  subroutine loss_gradient(grad,cost,exag)

    real(kind=dp), dimension(ldim,npoints), intent(out) :: grad
    real(kind=dp),                          intent(out) :: cost
    real(kind=dp),                          intent(in)  :: exag

    if(.not.allocated(vecs)) allocate(vec(ldim),vecs(ldim,npoints,npoints),rij2(npoints,npoints))

    ! * Cost function contribution to gradient
    if(use_perplexity) then
       if(use_kl) then
          if(.not.use_sgd) then
             call perp_kl_grad(grad,cost,exag)
          else
             write(stderr,*) 'ERROR: SGD not compatible with KL cost function'; stop
          end if
       else if(use_ce) then
          if(.not.use_sgd) then
             call perp_ce_grad(grad,cost,exag)
          else
             call perp_ce_sgd(grad,exag)
          end if
       end if
    else if(use_knn) then
       if(use_ce) then
          if(.not.use_sgd) then
             call knn_ce_grad(grad,cost,exag)
          else
             call knn_ce_sgd(grad,exag)
          end if
       end if
    end if

    ! * Hard sphere contribution to gradient
    if(core.and..not.use_sgd) then
       !$omp parallel do private(dist,vec) reduction(+:grad,cost) schedule(dynamic)
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          do nj=ni+1,npoints
             if(point_count(nj).eq.0) cycle
             rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
             dist=sqrt(rij2(nj,ni))
             if(dist.lt.point_radius(ni)+point_radius(nj)) then
                vec(:)=-vecs(:,nj,ni)/dist
                dist=(point_radius(ni)+point_radius(nj)-dist)/2.0_dp
                grad(:,ni)=grad(:,ni)+vec(:)*dist*core_strength/2.0_dp
                grad(:,nj)=grad(:,nj)-vec(:)*dist*core_strength/2.0_dp
                cost=cost+dist**2/2.0_dp*core_strength
             end if
          end do
       end do
       !$omp end parallel do
    end if

  end subroutine loss_gradient

  subroutine perp_kl_grad(grad,cost,exag)

    real(kind=dp), dimension(ldim,npoints), intent(out) :: grad
    real(kind=dp),                          intent(out) :: cost
    real(kind=dp),                          intent(in)  :: exag

    Z=0.0_dp
    cost=czero
    grad=0.0_dp
    !$omp parallel

    !$omp do reduction(+:Z) schedule(dynamic)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=ni+1,npoints
          if(point_count(nj).eq.0) cycle
          vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
          rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
          Z=Z+2.0_dp/(1.0_dp+rij2(nj,ni))
       end do
    end do
    !$omp end do

    !$omp do private(qij,vec) reduction(+:grad,cost) schedule(dynamic)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=ni+1,npoints
          if(point_count(nj).eq.0) cycle
          qij=1.0_dp/(1.0_dp+rij2(nj,ni))/Z
          vec(:)=4.0_dp*(exag*pij(nj,ni)-qij)*qij*vecs(:,nj,ni)*Z
          grad(:,ni)=grad(:,ni)+vec(:)
          grad(:,nj)=grad(:,nj)-vec(:)
          cost=cost-pij(nj,ni)*log(qij)*2.0_dp
       end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine perp_kl_grad

  subroutine perp_ce_grad(grad,cost,exag)

    real(kind=dp), dimension(ldim,npoints), intent(out) :: grad
    real(kind=dp),                          intent(out) :: cost
    real(kind=dp),                          intent(in)  :: exag

    Z=1/(real(reduced_point_count,dp)*(real(reduced_point_count,dp)-1.0_dp))
    cost=czero
    grad=0.0_dp
    !$omp parallel do private(qij,vec) reduction(+:grad,cost) schedule(dynamic)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=ni+1,npoints
          if(point_count(nj).eq.0) cycle
          vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
          rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
          qij=Z/(1.0_dp+rij2(nj,ni))
          vec(:)=4.0_dp*(exag*pij(nj,ni)-(1-pij(nj,ni))/(1-qij)*qij)*qij*vecs(:,nj,ni)/Z
          grad(:,ni)=grad(:,ni)+vec(:)
          grad(:,nj)=grad(:,nj)-vec(:)
          cost=cost-pij(nj,ni)*log(qij)*2.0_dp-(1-pij(ni,nj))*log(1-qij)*2.0_dp
       end do
    end do
    !$omp end parallel do

  end subroutine perp_ce_grad

  subroutine perp_ce_sgd(grad,exag)

    real(kind=dp), dimension(ldim,npoints), intent(out) :: grad
    real(kind=dp),                          intent(in)  :: exag

    Z=1.0_dp/(real(reduced_point_count,dp)*(real(reduced_point_count,dp)-1.0_dp))
    grad=0.0_dp
    !$omp parallel do private(ni,nj,qij,vec,temp,coeff) reduction(+:grad) schedule(dynamic)
    do n=1,n_simplices

       ! * Attraction
       ni=head(n)
       nj=tail(n)
       vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
       rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
       qij=1.0_dp/(1.0_dp+rij2(nj,ni))
    
       vec(:)=4.0_dp*min(1.0_dp,max(-1.0_dp,simplices(n)*qij*vecs(:,nj,ni))) ! Cap magnitude of each grad component
       grad(:,ni)=grad(:,ni)+vec(:)
       grad(:,nj)=grad(:,nj)-vec(:)
    
       ! * Repulsion
       do neg=1,nneg
    
         do
           nj=random_integer(npoints)
           if(nj.ne.ni) exit
         end do
         vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
         rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
         temp=1.0_dp+rij2(nj,ni)
         coeff=exag*(1.0_dp-pij(nj,ni))/((temp/Z-1.0_dp)*temp) ! (1-pij)*qij*qij/(1-qij)
    
         vec(:)=-4.0_dp*min(1.0_dp,max(-1.0_dp,coeff*vecs(:,nj,ni))) ! Cap magnitude of each grad component
         grad(:,ni)=grad(:,ni)+vec(:)
         grad(:,nj)=grad(:,nj)-vec(:)
    
       end do
    
    end do
    !$omp end parallel do

  end subroutine perp_ce_sgd

  subroutine knn_ce_grad(grad,cost,exag)

    real(kind=dp), dimension(ldim,npoints), intent(out) :: grad
    real(kind=dp),                          intent(out) :: cost
    real(kind=dp),                          intent(in)  :: exag

    cost=czero
    grad=0.0_dp
    !$omp parallel do private(qij,vec) reduction(+:grad,cost) schedule(dynamic)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=ni+1,npoints
          if(point_count(nj).eq.0) cycle
          vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
          rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
          qij=1.0_dp/(1.0_dp+rij2(nj,ni))
          vec(:)=4.0_dp*(exag*pij(nj,ni)-(1-pij(nj,ni))/(1-qij)*qij)*qij*vecs(:,nj,ni)
          grad(:,ni)=grad(:,ni)+vec(:)
          grad(:,nj)=grad(:,nj)-vec(:)
          if(pij(nj,ni).gt.epsilon(1.0_dp).and.pij(nj,ni).lt.(1.0_dp-epsilon(1.0_dp))) then
            cost=cost-pij(nj,ni)*log(qij)*2.0_dp-(1-pij(ni,nj))*log(1-qij)*2.0_dp
          else if(pij(nj,ni).gt.epsilon(1.0_dp).and.pij(nj,ni).gt.(1.0_dp-epsilon(1.0_dp))) then
            cost=cost-pij(nj,ni)*log(qij)*2.0_dp
          else if(pij(nj,ni).lt.epsilon(1.0_dp).and.pij(nj,ni).lt.(1.0_dp-epsilon(1.0_dp))) then
            cost=cost-(1-pij(ni,nj))*log(1-qij)*2.0_dp
          end if
       end do
    end do
    !$omp end parallel do

  end subroutine knn_ce_grad

  subroutine knn_ce_sgd(grad,exag)

    real(kind=dp), dimension(ldim,npoints), intent(out) :: grad
    real(kind=dp),                          intent(in)  :: exag

    grad=0.0_dp
    !$omp parallel do private(ni,nj,qij,vec,coeff) reduction(+:grad) schedule(dynamic)
    do n=1,n_simplices

       ! * Attraction
       ni=head(n)
       nj=tail(n)
       vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
       rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
       qij=1.0_dp/(1.0_dp+rij2(nj,ni))

       vec(:)=4.0_dp*exag*min(1.0_dp,max(-1.0_dp,simplices(n)*qij*vecs(:,nj,ni))) ! Cap magnitude of each grad component
       grad(:,ni)=grad(:,ni)+vec(:)
       grad(:,nj)=grad(:,nj)-vec(:)

       ! * Repulsion
       do neg=1,nneg

         do
           nj=random_integer(npoints)
           if(nj.ne.ni) exit
         end do
         vecs(:,nj,ni)=point_pos(:,ni)-point_pos(:,nj)
         rij2(nj,ni)=dot_product(vecs(:,nj,ni),vecs(:,nj,ni))
         coeff=(1.0_dp-pij(nj,ni))/(rij2(nj,ni)*(1.0_dp+rij2(nj,ni))) ! (1-pij)*qij*qij/(1-qij)
         vec(:)=-4.0_dp*min(1.0_dp,max(-1.0_dp,coeff*vecs(:,nj,ni))) ! Cap magnitude of each grad component
         grad(:,ni)=grad(:,ni)+vec(:)
         grad(:,nj)=grad(:,nj)-vec(:)

       end do

    end do
    !$omp end parallel do

  end subroutine knn_ce_sgd

end module loss
