!=================================================================================!
!                                     CONTROL                                     !
!=================================================================================!
!                                                                                 !
! This module is part of SHEAP --- controls the command line arguments            !
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

module control

  ! * Prerequisite modules
  use constants

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public subroutines
  public :: control_arguments
  public :: control_report

  ! * Public variables

  ! * Read control
  character(len=3),  public :: file_type='xyz'
  character(len=30), public :: meta_data=''
  character(len=30), public :: resume_file=''

  logical, public :: resume=.false.

  integer, public           :: npoints_max=120000

  ! * Output control
  logical, public           :: verbose   =.false.
  logical, public           :: quiet     =.false.
  logical, public           :: grad_test =.false.

  ! * Dimensionality of map
  integer, public           :: ldim=2

  ! * Normalize source data
  logical, public           :: normalize=.false.

  ! * Thresholds for combining structures
  real(kind=dp), public     :: sim_thresh=0.0_dp 
  real(kind=dp), public     :: energy_thresh=1e8_dp

  ! * Manifold learning scheme
  logical, public           :: use_perplexity =.true.
  logical, public           :: use_knn        =.false.

  integer, public           :: knn=15

  real(kind=dp), public     :: perplexity=15.0_dp

  ! * Embedding initialisation
  logical, public           :: init_only         =.false.
  logical, public           :: uniform_packing   =.false.
  logical, public           :: random_projection =.true.
  logical, public           :: pca               =.false.

  real(kind=dp), public     :: packing_fraction=0.5_dp
  real(kind=dp), public     :: projection_compression=10.0_dp

  ! * Cost function
  logical, public           :: use_kl     =.false.
  logical, public           :: use_ce     =.true.
  logical, public           :: use_stress =.false.

  ! * Optimisation scheme
  logical, public           :: use_tpsd   =.true.
  logical, public           :: use_sgd    =.false.
  logical, public           :: trace      =.false.
  logical, public           :: track      =.false.

  integer, public           :: trace_step=1
  integer, public           :: niter=10000
  integer, public           :: nexag=200
  integer, public           :: nneg=1

  real(kind=dp), public     :: sexag=5.0_dp
  real(kind=dp), public     :: alpha=-1.0_dp
  real(kind=dp), public     :: tolerance=1e-8_dp

  ! * Hard sphere growth    
  logical, public           :: grow =.false.

  integer, public           :: growth_steps=500

  real(kind=dp), public     :: sphere_radius=0.01_dp
  real(kind=dp), public     :: core_strength=1.0_dp
  real(kind=dp), public     :: growth_tol_coeff=2.0_dp

contains

  subroutine control_arguments()

    integer                                     :: i,iargc,num_args

    character(len=80), allocatable,dimension(:) :: cargs

    num_args = iargc()

    allocate(cargs(num_args))

    do i=1,num_args
       call getarg(i,cargs(i))
    end do

    ! * Read control
    if(any(cargs.eq."-read")) then               ! * Format of input file to be read
       do i=1,num_args-1
          if(cargs(i).eq."-read") exit
       end do
       read(cargs(i+1),*) file_type
    end if
    if(any(cargs.eq."-m")) then                  ! * Location of seperate meta-data file
       do i=1,num_args-1
          if(cargs(i).eq."-m") exit
       end do
       read(cargs(i+1),*) meta_data
    end if
    if(any(cargs.eq."-resume")) then             ! * Location of resume file
       do i=1,num_args-1
          if(cargs(i).eq."-resume") exit
       end do
       read(cargs(i+1),*) resume_file
       resume=.true.
    end if
    if(any(cargs.eq."-np")) then                 ! * Max number of points to read from .vec file
       do i=1,num_args-1
          if(cargs(i).eq."-np") exit
       end do
       read(cargs(i+1),*) npoints_max
    end if

    ! * Output control
    if(any(cargs.eq."-v")) verbose=.true.        ! * Verbose
    if(any(cargs.eq."-q")) quiet=.true.          ! * Quiet
    if(any(cargs.eq."-grad")) grad_test=.true.   ! * Test gradient against numerical scheme - does not carry out optimisation

    ! * Dimensionality of map                    ! * Dimensionality of map
    if(any(cargs.eq."-dim")) then
       do i=1,num_args-1
          if(cargs(i).eq."-dim") exit
       end do
       read(cargs(i+1),*) ldim
    end if

    ! * Normalize source data
    if(any(cargs.eq."-scale")) normalize=.true.  ! * Scale source data such that each variable has stddev of 1

    ! * Thresholds for combining structures
    if(any(cargs.eq."-st")) then                 ! * Structure similarity threshold
       do i=1,num_args-1
          if(cargs(i).eq."-st") exit
       end do
       read(cargs(i+1),*) sim_thresh
    end if
    if(any(cargs.eq."-et")) then                 ! * Cost threshold
       do i=1,num_args-1
          if(cargs(i).eq."-et") exit
       end do
       read(cargs(i+1),*) energy_thresh
    end if

    ! * Standard embedding algorithms
    if(any(cargs.eq."-tsne")) then               ! * t-SNE
       use_perplexity=.true.
       use_knn=.false.
       use_kl=.true.
       use_ce=.false.
       use_stress=.false.
       use_tpsd=.true.
       use_sgd=.false.
       nexag=250
       sexag=4
    end if
    if(any(cargs.eq."-umap")) then               ! * UMAP
       use_perplexity=.false.
       use_knn=.true.
       use_kl=.false.
       use_ce=.true.
       use_stress=.false.
       use_tpsd=.false.
       use_sgd=.true.
    end if

    ! * Manifold learning scheme
    if(any(cargs.eq."-p")) then                  ! * Use perplexity algorithm for computing weights
       use_perplexity=.true.
       use_knn=.false.
       do i=1,num_args-1
          if(cargs(i).eq."-p") exit
       end do
       read(cargs(i+1),*) perplexity
    end if
    if(any(cargs.eq."-k")) then                  ! * Use K-NN algorithm for computing weights
       use_perplexity=.false.
       use_knn=.true.
       do i=1,num_args-1
          if(cargs(i).eq."-k") exit
       end do
       read(cargs(i+1),*) knn
    end if

    ! * Embedding initialisation
    if(any(cargs.eq."-init")) then               ! * Stop after initialising low-dimensional embedding - for doing just pca or random projection
       init_only=.true.
    end if
    if(any(cargs.eq."-up")) then                 ! * Use random uniform packing to initialise low-dimensional embedding
       uniform_packing=.true.
       random_projection=.false.
       pca=.false.
    end if
    if(any(cargs.eq."-f")) then                  ! * Packing fraction when using random uniform packing to initialise low-dimensional embedding
       uniform_packing=.true.
       random_projection=.false.
       pca=.false.
       do i=1,num_args-1
          if(cargs(i).eq."-f") exit
       end do
       read(cargs(i+1),*) packing_fraction
    end if
    if(any(cargs.eq."-pca")) then                ! * Use pca to initialise low-dimensional embedding
       uniform_packing=.false.
       random_projection=.false.
       pca=.true.
    end if
    if(any(cargs.eq."-pc")) then                 ! * Compression factor to apply to initialisation of low-dimensional embedding
       do i=1,num_args-1
          if(cargs(i).eq."-pc") exit
       end do
       read(cargs(i+1),*) projection_compression
    end if

    ! * Cost function
    if(any(cargs.eq."-kl")) then                 ! * Use KL-divergence loss function
       use_kl=.true.
       use_ce=.false.
       use_stress=.false.
    end if
    if(any(cargs.eq."-ce")) then                 ! * Use cross-entropy loss function
       use_kl=.false.
       use_ce=.true.
       use_stress=.false.
    end if

    ! * Optimisation scheme
    if(any(cargs.eq."-tpsd")) then               ! * Use two-point steepest descent algorithm
       use_tpsd=.true.
       use_sgd=.false.
    end if
    if(any(cargs.eq."-sgd")) then                ! * Use stochastic gradient descent
       use_tpsd=.false.
       use_sgd=.true.
    end if
    if(any(cargs.eq."-t")) trace=.true.          ! * Trace the optimisation trajectory
    if(any(cargs.eq."-ts")) then                 ! * Output every trace_step steps
       do i=1,num_args-1
          if(cargs(i).eq."-ts") exit
       end do
       read(cargs(i+1),*) trace_step
    end if
    if(any(cargs.eq."-track")) track=.true.      ! * Track the optimisation - produces sheap.track file containing progression of cost, gradient, etc. during optimisation
    if(any(cargs.eq."-n")) then                  ! * Maximum number of iterations
       do i=1,num_args-1
          if(cargs(i).eq."-n") exit
       end do
       read(cargs(i+1),*) niter
    end if
    if(any(cargs.eq."-ne")) then                 ! * Number of steps to apply early exaggeration
       do i=1,num_args-1
          if(cargs(i).eq."-ne") exit
       end do
       read(cargs(i+1),*) nexag
       if(nexag.gt.niter) then
          write(stderr,*) "Error: number of exaggerated steps cannot be greater than total number of steps."
          stop
       end if
    end if
    if(any(cargs.eq."-se")) then                 ! * Size of exaggeration
       do i=1,num_args-1
          if(cargs(i).eq."-se") exit
       end do
       read(cargs(i+1),*) sexag
    end if
    if(any(cargs.eq."-neg")) then                ! * Number of negative samples per positive sample when using SGD
       do i=1,num_args-1
          if(cargs(i).eq."-neg") exit
       end do
       read(cargs(i+1),*) nneg
    end if
    if(any(cargs.eq."-al")) then                 ! * Initial alpha for SGD
       do i=1,num_args-1
          if(cargs(i).eq."-al") exit
       end do
       read(cargs(i+1),*) alpha
    end if
    if(any(cargs.eq."-tol")) then                ! * Gradient tolerance for convergence of optimisation
       do i=1,num_args-1
          if(cargs(i).eq."-tol") exit
       end do
       read(cargs(i+1),*) tolerance
    end if

    ! * Hard sphere growth
    if(any(cargs.eq."-hs")) then                 ! * Grow hard spheres
       grow=.true.
    end if
    if(any(cargs.eq."-cs")) then                 ! * Core strength
       do i=1,num_args-1
          if(cargs(i).eq."-cs") exit
       end do
       read(cargs(i+1),*) core_strength
    end if
    if(any(cargs.eq."-gs")) then                 ! * Number of steps over which to grow hard spheres
       do i=1,num_args-1
          if(cargs(i).eq."-gs") exit
       end do
       read(cargs(i+1),*) growth_steps
       if(growth_steps.le.0) then
          grow=.false.
       end if
    end if
    if(any(cargs.eq."-grtol")) then              ! * Coefficient to scale tolerance by to produce cost error below which to start growing hard spheres
       do i=1,num_args-1
          if(cargs(i).eq."-grtol") exit
       end do
       read(cargs(i+1),*) growth_tol_coeff
    end if
    if(any(cargs.eq."-rs")) then                 ! * Sphere radius
       do i=1,num_args-1
          if(cargs(i).eq."-rs") exit
       end do
       read(cargs(i+1),*) sphere_radius
    else
       sphere_radius=sphere_radius*2.0_dp**(real(ldim,dp)-2.0_dp)
    end if

    ! * Help argument
    if(any(cargs.eq."-h")) goto 100

    deallocate(cargs)

    return

    ! * Help message

    ! * Usage
100 write (*,*) 'Usage    : sheap [OPTIONS] < in.xyz 2> [MESSAGES] > [out.xyz]'

    ! * Header
    write (*,*)
    write (*,*) 'Options list'
    write (*,*) '============'

    ! * Help message
    write (*,*)
    write (*,*) 'Help'
    write (*,*) '-h       : Output this message'

    ! * Input control

    write(*,*)
    write(*,*) 'Input control'
    write(*,*) '-read     : Format of input file. Options: vec, old, xyz [xyz]'
    write(*,*) '-m        : Location of seperate meta-data file, format:'
    write(*,*) '            ...'
    write(*,*) '            file_label num_atoms formula "symmetry" volume energy count'
    write(*,*) '            ...'

    ! * Output control
    write (*,*)
    write (*,*) 'Output control'
    write (*,*) '-v       : Verbose'
    write (*,*) '-q       : Quiet'
    write (*,*) '-grad    : Test gradient against numerical scheme - does not carry out optimisation'

    ! * Dimensionality of map
    write (*,*)
    write (*,*) 'Dimensionality of map'
    write (*,*) '-dim N   : Dimensionality of map [2]'

    ! * Normalize source data
    write(*,*)
    write(*,*) 'Normalize source data such that each variable has stddev of 1'
    write(*,*) '-scale    : Set normalize to true [false]'

    ! * Thresholds for combining structures
    write (*,*)
    write (*,*) 'Thresholds for combining structures'
    write (*,*) '-st R    : Similarity threshold [0.0]'
    write (*,*) '-et R    : Cost threshold [0.0]'

    ! * Standard algorithms
    write (*,*)
    write (*,*) 'Standard algorithms'
    write (*,*) '-tsne    : Use standard t-SNE algorithm'
    write (*,*) '-umap    : Use standard UMAP algorithm'

    ! * Manifold learning scheme
    write (*,*)
    write (*,*) 'Manifold learning scheme'
    write (*,*) '-p R     : Use perplexity algorithm to compute weights and set perplexity [15.0]'
    write (*,*) '-k N     : Use K-NN algorithm to compute weights and set K [15]'

    ! * Embedding initialisation
    write (*,*)
    write (*,*) 'Embedding initialisation'
    write (*,*) '-init    : Stop after initialising low-dimensional embedding',&
                            ' - for doing random projection or pca'
    write (*,*) '-up      : Use random uniform packing to initialise low-dimensional embedding',&
                            ' [default is to use random projection]'
    write (*,*) '-f R     : Packing fraction for random uniform packing initialisation of low-dimensional embedding [0.5]'
    write (*,*) '-pca     : Use pca to initialise low-dimensional embedding',&
                            ' [default is to use random projection]'
    write (*,*) '-pc R    : Compression factor to apply to initialisation of low_dimensional embedding [10.0]'

    ! * Cost function
    write (*,*)
    write (*,*) 'Cost function'
    write (*,*) '-kl      : Use KL-divergence loss function - only compatable with perplexity initialisation of weights'
    write (*,*) '-ce      : Use cross-entropy loss function [default]'

    ! * Optimisation scheme
    write (*,*)
    write (*,*) 'Optimisation scheme'
    write (*,*) '-tpsd    : Use two point steepest descent optimisation algorithm [default]'
    write (*,*) '  -tol R : Tolerance for optimisation [1e-8]'
    write (*,*) '  -ne N  : Number of steps to apply early exaggeration [200]'
    write (*,*) '  -se R  : Size of exaggeration [5.0]'
    write (*,*) '  -track : Track the optimisation [false]'
    write (*,*) '            - generates sheap.track file'
    write (*,*) '-sgd     : Use stochastic gradient descent optimisation algorithm'
    write (*,*) '  -neg N : Number of negative samples per positive sample [1]'
    write (*,*) '  -al R  : Alpha [1.0]'
    write (*,*) '            - initial step size'
    write (*,*) '-n  N    : Maximum number of iterations [10000]'
    write (*,*) '            - sets total for SGD'
    write (*,*) '-t       : Trace the optimisation [false]'
    write (*,*) '-ts N    : Output layout every number of steps [1]'

    ! * Hard sphere growth
    write (*,*)
    write (*,*) 'Hard sphere growth (TPSD only)'
    write (*,*) '-hs      : Grow hard spheres'
    write (*,*) '-cs R    : Hard sphere core strength [1.0]'
    write (*,*) '-gs N    : Number of steps over which to grow hard spheres [500]'
    write (*,*) '-grtol R : Tolerence for starting hard sphere growth (multiplies -tol) [2.0]'
    write (*,*) '            - must be greater than 1.0'
    write (*,*) '-rs R    : Sphere radius [0.01 (2D), 0.02 (3D)]'
    write (*,*) '            - if negative, computes suitable value from map layout at onset of growth'

    stop

  end subroutine control_arguments

  subroutine control_report()

    ! * Header
    write (stderr,*) 'control parameters'
    write (stderr,*)

    ! * Dimensionality of map
    write (stderr,'(a25,a3,i5)') 'ldim',' : ',ldim

    ! * Thresholds for combining structures
    write (stderr,'(a25,a3,f10.7)') 'sim_thresh',' : ',sim_thresh
    write (stderr,'(a25,a3,e10.5)') 'energy_thresh',' : ',energy_thresh

    ! * Manifold learning scheme
    if(use_perplexity) write (stderr,'(a25,a3,f10.4)') 'perplexity',' : ',perplexity
    if(use_knn) write (stderr,'(a25,a3,i5)') 'knn',' : ',knn

    ! * Embedding initialisation
    if(uniform_packing) write (stderr,'(a25,a3,f10.4)') 'packing_fraction',' : ',packing_fraction
    if(random_projection) write (stderr,'(a25,a3,f10.4)') 'projection_compression',' : ',projection_compression

    ! * Optimisation scheme
    write (stderr,'(a25,a3,l5)') 'trace',' : ',trace
    if(trace) write (stderr,'(a25,a3,l5)') 'trace_step',' : ',trace_step
    write (stderr,'(a25,a3,l5)') 'track',' : ',track
    write (stderr,'(a25,a3,i5)') 'niter',' : ',niter
    write (stderr,'(a25,a3,i5)') 'nexag',' : ',nexag
    if(nexag.gt.0) write (stderr,'(a25,a3,f10.4)') 'sexag',' : ',sexag
    if(use_sgd) write (stderr,'(a25,a3,i5)') 'nneg',' : ',nneg
    if(alpha.gt.0)write (stderr,'(a25,a3,f10.4)') 'alpha',' : ',alpha
    write (stderr,'(a25,a3,e10.4)') 'tolerance',' : ',tolerance

    ! * Hard sphere growth
    write (stderr,'(a25,a3,l5)') 'grow',' : ',grow
    if(grow) write (stderr,'(a25,a3,i5)') 'growth_steps',' : ',growth_steps
    if(grow) write (stderr,'(a25,a3,f10.4)') 'core_strength',' : ',core_strength
    if(grow) write (stderr,'(a25,a3,f10.4)') 'growth_tol_coeff',' : ',growth_tol_coeff
    if(sphere_radius.gt.0) write (stderr,'(a25,a3,f10.4)') 'sphere_radius',' : ',sphere_radius

    write (stderr,*)

  end subroutine control_report

end module control
