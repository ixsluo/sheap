!=================================================================================!
!                                      OPTIM                                      !
!=================================================================================!
!                                                                                 !
! This module is part of SHEAP --- routines for optimising embedding              !
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

module optim

  ! * Prerequisite modules
  use random
  use constants
  use control
  use model
  use loss

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public subroutines
  public :: optim_ise
  public :: optim_test

  ! * Public variables
  integer, public       :: steps,total_steps=0
  integer, public       :: start_growth

  real(kind=dp), public :: final_cost

  ! * Private variables
  integer, parameter :: unit_track=21
  
contains

  subroutine optim_ise(thresh,maxsteps)

    integer,       intent(in) :: maxsteps
    integer                   :: ni,point_counter,e_num,info,iwork(5*ldim),ifail(ldim)


    real(kind=dp), intent(in) :: thresh
    real(kind=dp)             :: com(ldim),point_red(ldim,reduced_point_count),cov(ldim,ldim)
    real(kind=dp)             :: e_val(ldim),e_vec(ldim,ldim),vl,vu,lwork
    real(kind=dp), allocatable, dimension(:) :: work

    if(track) open(unit=unit_track,file='sheap.track',form='formatted',status='unknown')

    ! * Call optimisation scheme of choice
    if(use_tpsd) then
       call tpsd(thresh,maxsteps)
    else if(use_sgd) then
       call sgd(thresh,maxsteps)
    end if

    if(verbose) write(stderr,*) "doing PCA on result"

    ! * Compute centre-of-mass of map
    com=0.0_dp
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       com=com+point_pos(:,ni)
    end do
    com=com/real(reduced_point_count,dp)

    ! * Move CoM to the origin
    sphere_radius=0.0_dp
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       point_pos(:,ni)=point_pos(:,ni)-com(:)
    end do

    ! * Do PCA on result
    if(reduced_point_count.ne.npoints) then
       point_counter=0
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          point_counter=point_counter+1
          point_red(:,point_counter)=point_pos(:,ni)
       end do
       cov=matmul(point_red(1:ldim,1:reduced_point_count),transpose(point_red(1:ldim,1:reduced_point_count)))
    else
       cov=matmul(point_pos(1:ldim,1:reduced_point_count),transpose(point_pos(1:ldim,1:reduced_point_count)))
    end if
    
    ifail=0
    call dsyevx('V','I','U',ldim,cov,ldim,vl,vu,1,ldim,epsilon(1.0),&
            e_num,e_val,e_vec,ldim,lwork,-1,iwork,ifail,info)
    allocate(work(int(lwork)))
    
    ifail=0
    call dsyevx('V','I','U',ldim,cov,ldim,vl,vu,1,ldim,epsilon(1.0),&
            e_num,e_val,e_vec,ldim,lwork,int(lwork),iwork,ifail,info)
    
    point_pos=matmul(transpose(e_vec(1:ldim,ldim:1:-1)),point_pos(1:ldim,1:npoints))

    if(track) close(unit_track)

  end subroutine optim_ise

  subroutine sgd(thresh,maxsteps)

    integer,       intent(in) :: maxsteps

    real(kind=dp), intent(in) :: thresh
    real(kind=dp)             :: cost,initial_alpha,exag
    real(kind=dp)             :: gmat(ldim,npoints),gvec(ldim*npoints),pos_vec(ldim*npoints)

    if(.not.quiet) write (stderr,'(a)',advance='NO') ' '

    ! * Initialise variables
    steps=0
    exag=1.0_dp

    if(alpha.gt.0) then
      initial_alpha=alpha
    else if(use_perplexity) then
      initial_alpha=real(reduced_point_count,dp)
    else
      initial_alpha=1.0_dp
    end if
    alpha=initial_alpha

    pos_vec=reshape(point_pos,(/ldim*npoints/))

    ! * If tracing optimisation, write initial map layout output
    if(trace) then
       call loss_gradient(gmat,cost,exag)
       gvec=reshape(gmat,(/ldim*npoints/))
       call model_write(cost)
    end if

    ! * Optimisation loop
    do while(steps.lt.(maxsteps))

       ! * Increment optimisation step 
       steps=steps+1 ; total_steps=total_steps+1

       ! * If using perplexity scheme with sgd, switch to full gradient optimisation with tpsd following final step of sgd
       if(use_perplexity.and.steps.ge.maxsteps) then
         use_sgd=.false.
         nexag=0
         sexag=1.0_dp
         if(verbose) then
            write(stderr,'(a)') ''
            write(stderr,'(a)') ' doing tpsd'
         end if
         call tpsd(thresh,maxsteps) 
         exit
       end if

       ! * Compute gradient
       call loss_gradient(gmat,cost,exag)
       gvec=reshape(gmat,(/ldim*npoints/))

       ! * Update map points
       pos_vec=pos_vec-gvec*alpha

       point_pos=reshape(pos_vec,(/ldim,npoints/))

       ! * If tracing optimisation, write updated map points to output
       if(trace.and.(mod(steps-1,trace_step).eq.0)) then
          call model_write(cost)
       end if

       ! * Update step size
       alpha=initial_alpha*(1.0_dp-real(steps)/real(maxsteps))

       ! * Write simple illustration of progress of optimisation to stderr
       if(.not.quiet.and.mod(steps,100).eq.0) write (stderr,'(a)',advance='NO') '.'
       if(.not.quiet.and.steps.ge.(maxsteps)) write (stderr,'(a)',advance='YES') ''

    end do

    ! * Compute and update final cost of map layout 
    use_sgd=.false.
    call loss_gradient(gmat,cost,exag)
    final_cost=cost
    use_sgd=.true.

    if(.not.quiet) write (stderr,*)

  end subroutine sgd

  subroutine tpsd(thresh,maxsteps)

    logical                   :: still_growing=.false.

    integer,       intent(in) :: maxsteps
    integer                   :: ni,glogint

    real(kind=dp), intent(in) :: thresh
    real(kind=dp)             :: cost,step,sgg,grun,exag
    real(kind=dp)             :: gmat(ldim,npoints),gvec(ldim*npoints),gptb(ldim*npoints)
    real(kind=dp)             :: pos_vec(ldim*npoints),ptb_vec(ldim*npoints),com(ldim)

    if(.not.quiet) write (stderr,'(a)',advance='NO') ' '

    ! * Initialise variables
    steps=0
    sgg=huge(1.0_dp)
    grun=0.0_dp
    glogint=huge(1)
    exag=sexag

    pos_vec=reshape(point_pos,(/ldim*npoints/))

    ! * If tracing optimisation, write initial map layout to output
    if(trace) then
       call loss_gradient(gmat,cost,exag)
       gvec=reshape(gmat,(/ldim*npoints/))
       call model_write(cost)
    end if

    ! * Prepare for hard sphere growth, if enabled
    if(grow) still_growing=.true.

    ! * Optimisation loop
    do while((((grun.gt.log10(thresh)).or.(steps.lt.(nexag+100)).or.(steps.lt.100))&
               .and.((steps.lt.maxsteps).or.(maxsteps.lt.0))).or.still_growing)

       ! * Increment optimisation step
       steps=steps+1 ; total_steps=total_steps+1

       ! * Switch off exaggeration following nexag steps
       if(steps.ge.nexag) exag=1.0_dp

       ! * Compute gradient
       call loss_gradient(gmat,cost,exag)
       gvec=reshape(gmat,(/ldim*npoints/))

       ! * Apply small perturbation to gradient vector
       call random_point_in_hypersphere(ptb_vec,1e-2_dp)        

       gptb=gvec*(1.0_dp+ptb_vec)

       ! * Compute step size
       step=abs(tpsd_step(pos_vec,gptb,init=((steps.eq.1).or.(steps.eq.nexag))))

       ! * Update map points
       pos_vec=pos_vec-step*gptb

       point_pos=reshape(pos_vec,(/ldim,npoints/))

       ! * Compute grun - for monitoring convergence
       sgg=abs(dot_product(step*gptb,gvec))

       grun=grun+(log10(sgg)-grun)/min(steps,100)

       ! * Write simple illustration of progress of optimisation to stderr
       if((int(grun).lt.glogint).and.(steps.gt.10)) then
          glogint=int(grun)
          if(.not.quiet) write (stderr,'(a)',advance='NO') ':'
       end if

       if((int(grun).gt.glogint).and.(steps.gt.10)) then
          glogint=int(grun)
          if(.not.quiet) write (stderr,'(a)',advance='NO') '^'
       end if

       ! * If tracing optimisation, write updated map points to output
       if(trace.and.(mod(steps-1,trace_step).eq.0)) then
          call model_write(cost)
       end if

       ! * If tracking optimisation, write to sheap.track
       if(track) then
          if(steps.eq.1) write (unit_track,'(a10,a25,a25,a25,a25)') 'steps','cost','gmag','log10(sgg)','grun'
          write (unit_track,'(i10,f25.20,f25.20,f25.20,f25.20)') steps,cost,sqrt(dot_product(gvec,gvec)),log10(sgg),grun
          flush(unit_track)
       end if

       ! * Set HS radii for start of growth period
       if(grow.and.(.not.core).and.(grun.lt.log10(thresh*growth_tol_coeff)).and.steps.ge.(nexag+100)) then
!       if(grow.and.(.not.core).and.(grun.lt.log10(thresh*growth_tol_coeff))) then
          core=.true.
          start_growth=steps
          if(.not.quiet) write (stderr,'(a)') '|'

          ! * If input minimum sphere radius was negative, compute value from the layout of map points at the start of growth
          if(recalc_spheres) then

             if(.not.quiet) then
                write (stderr,'(a)') ' computing minimum hard sphere radius'
                write (stderr,*)
             end if

             ! * Compute centre-of-mass of map
             com=0.0_dp
             do ni=1,npoints
                if(point_count(ni).eq.0) cycle
                com=com+point_pos(:,ni)
             end do
             com=com/real(reduced_point_count,dp)

             ! * Move CoM to the origin and compute sphere radius as average distance to CoM
             sphere_radius=0.0_dp
             do ni=1,npoints
                if(point_count(ni).eq.0) cycle
                point_pos(:,ni)=point_pos(:,ni)-com(:)
                sphere_radius=sphere_radius+dot_product(point_pos(:,ni),point_pos(:,ni))
             end do
             sphere_radius=sqrt(sphere_radius)/reduced_point_count*(real(reduced_point_count)&
                            /real(total_point_count))**(1.0_dp/real(ldim))/10.0_dp*2.0_dp**(ldim)

             ! * Write computed minimum sphere radius to stderr
             if(.not.quiet) then
                write (stderr,'(a25,a3,f14.8)') 'sphere_radius',' : ', sphere_radius
                write (stderr,*)
             end if

          end if

          ! * Compute sphere radii at start of HS growth
          do ni=1,npoints
             if(point_count(ni).eq.0) cycle
             point_radius(ni)=sphere_radius*real(point_count(ni),dp)**(1.0_dp/real(ldim,dp))
          end do
          point_radius=point_radius/100.0_dp

          ! * Report beginning of HS growth regime
          if(.not.quiet) then
             write (stderr,'(a)') ' continuing tpsd with hard sphere growth'
             write (stderr,'(a)',advance='NO') ' '
          end if

       end if

       ! * Increase size of spheres at each step during growth period
       if(core.and.(steps.gt.(start_growth)).and.(steps.lt.(start_growth+growth_steps))) then
          if(steps.lt.(start_growth+2)) then
             point_radius=point_radius*((real(steps)-real(start_growth))*100.0_dp/real(growth_steps))
          else
             point_radius=point_radius*((real(steps)-real(start_growth))*100.0_dp/real(growth_steps))&
                                  /((real(steps)-1.0_dp-real(start_growth))*100.0_dp/real(growth_steps))
          end if
       end if

       if(core.and.steps.ge.(start_growth+growth_steps+100)) still_growing=.false.

    end do

    ! * Store final map cost
    final_cost=cost

    if(.not.quiet) write (stderr,*)

  end subroutine tpsd

  function tpsd_step(v,g,init)

    logical, optional, intent(in)                  :: init
    logical                                        :: initialise

    real(kind=dp), dimension(:), intent(in)        :: v,g
    real(kind=dp), save, allocatable, dimension(:) :: v0,g0
    real(kind=dp), dimension(size(v))              :: dv,dg
    real(kind=dp)                                  :: tpsd_step
    real(kind=dp), parameter                       :: default_step=1e-1_dp

    tpsd_step=-1.0_dp

    if(present(init)) then
       initialise=init.or.(.not.allocated(v0))
    else
       initialise=.not.allocated(v0)
    end if

    if(initialise) then

       if(allocated(v0)) deallocate(v0,g0)
       allocate(v0(size(v)),g0(size(v)))
       tpsd_step=default_step 

    else

       dv=v-v0
       dg=g-g0

       tpsd_step=tpsd_alpha(dv,dg)

    end if

    v0 = v
    g0 = g

  end function tpsd_step

  function tpsd_alpha(dvec,dgrd)

    real(kind=dp)                           :: dgdg,dvdg,tpsd_alpha
    real(kind=dp), dimension(:), intent(in) :: dvec,dgrd

    dgdg=dot_product(dgrd,dgrd)
    dvdg=dot_product(dvec,dgrd)

    tpsd_alpha=dvdg/dgdg

  end function tpsd_alpha

  subroutine optim_test()

    integer                                    :: ni,nj

    real(kind=dp)                              :: cost,eplus,eminus,step=1e-8_dp,exag
    real(kind=dp), dimension(:,:), allocatable :: grad,gnum,pos_tmp

    allocate(grad(ldim,npoints),gnum(ldim,npoints),pos_tmp(ldim,npoints))

    pos_tmp=point_pos

    exag=1.0_dp

    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=1,ldim

          point_pos=pos_tmp
          point_pos(nj,ni)=pos_tmp(nj,ni)+step
          call loss_gradient(grad,eplus,exag)

          point_pos=pos_tmp
          point_pos(nj,ni)=pos_tmp(nj,ni)-step
          call loss_gradient(grad,eminus,exag)

          gnum(nj,ni)=(eplus-eminus)/2.0_dp/step
          
       end do

    end do

    point_pos=pos_tmp

    call loss_gradient(grad,cost,exag)

    write (*,*) 'grad' 
    write (*,'(3f20.10)') grad
    write (*,*)
    write (*,*) 'gnum'
    write (*,'(3f20.10)') gnum
    write (*,*)
    write (*,*) 'grad/gnum'
    write (*,'(3f20.10)') grad/gnum
    write (stderr,*)
    write (stderr,'(a25,a3,e25.15)') '|gnum-grad|',' : ',norm2(gnum-grad)
    write (stderr,'(a25,a3,e25.15)') '|gnum-grad|/|grad|',' : ',norm2(gnum-grad)/norm2(grad)
    write (stderr,'(a25,a3,f25.15)') 'cost',' : ',cost
    
  end subroutine optim_test

end module optim
