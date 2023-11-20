!=================================================================================!
!                                      MODEL                                      !
!=================================================================================!
!                                                                                 !
! This module is part of SHEAP --- routines for reading input and intialising     !
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

module model

  ! * Prerequisite modules
  use random
  use constants
  use control
  use rng

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public subroutines
  public :: model_read_vec
  public :: model_read_old
  public :: model_read_xyz
  public :: model_write
  public :: model_init
  public :: model_uniform_packing
  public :: model_random_projection
  public :: model_pca
  public :: model_resume
  public :: model_report
   
  ! * Public variables
  logical,       public                              :: recalc_spheres=.false.

  integer,       public                              :: npoints,reduced_point_count,total_point_count,n_simplices
  integer,       public, allocatable, dimension(:)   :: point_count,head,tail
  integer,       public, allocatable, dimension(:,:) :: nn_vec_indx

  real(kind=dp), public                              :: sig_ave,czero
  real(kind=dp), public, allocatable, dimension(:)   :: simplices,sigma,point_radius
  real(kind=dp), public, allocatable, dimension(:,:) :: nn_vec_dist,pij,point_pos

  ! * Private variables
  integer                                    :: hdim,pairs

  real(kind=dp), allocatable, dimension(:)   :: hcom,stddev
  real(kind=dp), allocatable, dimension(:,:) :: data_vec,dij2

  ! * AIRSS meta-data
  integer,    public, allocatable, dimension(:) :: ion_num

  real(kind=dp),      allocatable, dimension(:) :: ion_energy,ion_volume

  character(len=240), allocatable, dimension(:) :: ion_label
  character(len=40),  allocatable, dimension(:) :: ion_symmetry,ion_formula

  !! * For paper: computes #contacts with distances descriptor
  real(kind=dp),      allocatable, dimension(:) :: ion_nn
  real(kind=dp)                                 :: ion_nn_dist
    
contains

  subroutine model_read_vec()
    
    integer :: ni,length,min_length
    integer                       :: indx,stat,met2,met7
    character(len=50)             :: met1,met3,met4
    real(kind=dp)                 :: met5,met6

    min_length=huge(1)
   
    npoints=0
    do ni=1,npoints_max

       ! * Read length of vector, update minimum length if necessary
       read(stdin,*,end=1001) length
       if(length.lt.min_length) min_length=length
       if(ni.eq.1) then
          allocate(data_vec(min_length,npoints_max),ion_energy(npoints_max),ion_volume(npoints_max),ion_label(npoints_max),&
               ion_symmetry(npoints_max),ion_formula(npoints_max),ion_num(npoints_max),point_count(npoints_max))
       end if

       ! * Read descriptor vector
       data_vec(:,ni)=0.0_dp
       read(stdin,*) data_vec(1:min_length,ni)
       npoints=npoints+1

       ! * Read meta-data from seperate file, if specified
       if(len_trim(meta_data).gt.0) then

          read(stdin,*) ion_label(ni)

          open(unit=44,file=meta_data,iostat=stat,status='old')
          if(stat.ne.0) then
             write(stderr,*) 'ERROR: meta-data file does not exist'
             stop
          end if
          indx=0
          do
             indx=indx+1
             read(44,*,end=1003) met1,met2,met3,met4,met5,met6,met7
             if(trim(adjustl(met1)).eq.trim(adjustl(ion_label(ni)))) then
                ion_num(ni)=met2
                ion_formula(ni)=met3
                ion_symmetry(ni)=met4
                ion_volume(ni)=met5
                ion_energy(ni)=met6
                point_count(ni)=met7
                exit
             end if
          end do
          1003 continue
          close(44)

       ! * Read meta-data from .vec file
       else

          read(stdin,*) ion_label(ni),ion_num(ni),ion_formula(ni),ion_symmetry(ni),ion_volume(ni),ion_energy(ni),point_count(ni)

       end if
       ion_volume(ni)=ion_volume(ni)/real(ion_num(ni),dp)
       ion_energy(ni)=ion_energy(ni)/real(ion_num(ni),dp)
    end do
    1001 continue

    if(npoints.eq.npoints_max) then
       write(stderr,*) 'ERROR: number of points in .vec file exceeds maximum number allocated. Increase this with -np.'
    end if

    hdim=min_length
    allocate(point_pos(ldim,npoints),point_radius(npoints),dij2(npoints,npoints))
   
    total_point_count=sum(point_count(1:npoints))

  end subroutine model_read_vec

  subroutine model_read_old()

    integer :: ni
    integer                       :: indx,stat,met2,met7
    character(len=50)             :: met1,met3,met4
    real(kind=dp)                 :: met5,met6

    ! * Read number of points and dimensionality of original space
    read(stdin,*) npoints,hdim

    allocate(data_vec(hdim,npoints),point_pos(ldim,npoints),ion_energy(npoints),ion_volume(npoints),ion_label(npoints),&
         ion_symmetry(npoints),ion_formula(npoints),ion_num(npoints),point_radius(npoints),point_count(npoints),&
         dij2(npoints,npoints))

    ! * Read descriptor vectors
    do ni=1,npoints
       read(stdin,*) data_vec(:,ni)
    end do

    ! * Read meta-data
    do ni=1,npoints

       ! * Read meta-data from seperate file, if specified
       if(len_trim(meta_data).gt.0) then

          read(stdin,*) ion_label(ni)

          open(unit=44,file=meta_data,iostat=stat,status='old')
          if(stat.ne.0) then
             write(stderr,*) 'ERROR: meta-data file does not exist'
             stop
          end if
          indx=0
          do
             indx=indx+1
             read(44,*,end=1003) met1,met2,met3,met4,met5,met6,met7
             if(trim(adjustl(met1)).eq.trim(adjustl(ion_label(ni)))) then
                ion_num(ni)=met2
                ion_formula(ni)=met3
                ion_symmetry(ni)=met4
                ion_volume(ni)=met5
                ion_energy(ni)=met6
                point_count(ni)=met7
                exit
             end if
          end do
          1003 continue
          close(44)

       ! * Read meta-data from .xyz file
       else

          read(stdin,*) ion_label(ni),ion_num(ni),ion_formula(ni),ion_symmetry(ni),ion_volume(ni),ion_energy(ni),point_count(ni)

       end if
       ion_volume(ni)=ion_volume(ni)/real(ion_num(ni),dp)
       ion_energy(ni)=ion_energy(ni)/real(ion_num(ni),dp)
    end do

    total_point_count=sum(point_count)

  end subroutine model_read_old

  subroutine model_read_xyz()

    integer                       :: i,ni,indx,quote1,quote2,counts(15),stat,met2,met7
    character(len=:), allocatable :: line
    character(len=4)              :: element,elements(15)
    character(len=50)             :: met1,met3,met4
    real(kind=dp)                 :: met5,met6,lattice(9)
    real(kind=dp), allocatable, dimension(:) :: soap1,soap2

    allocate(character(len=10000000) :: line)

    allocate(ion_energy(npoints_max),ion_volume(npoints_max),ion_label(npoints_max),&
               ion_symmetry(npoints_max),ion_formula(npoints_max),ion_num(npoints_max),point_count(npoints_max))

    npoints=0
    do ni=1,npoints_max
       read(stdin,*,end=1002) ion_num(ni)
       npoints=npoints+1

       read(stdin,'(A)') line

       ! * Read structure label
       indx=index(line, "name=")
       read(line(indx+5:),*) ion_label(ni)

       ! * Read meta-data from seperate file, if specified
       if(len_trim(meta_data).gt.0) then

          open(unit=44,file=meta_data,iostat=stat,status='old')
          if(stat.ne.0) then
             write(stderr,*) 'ERROR: meta-data file does not exist'
             stop
          end if
          indx=0
          do
             indx=indx+1
             read(44,*,end=1003) met1,met2,met3,met4,met5,met6,met7
             if(trim(adjustl(met1)).eq.trim(adjustl(ion_label(ni)))) then
                ion_label(ni)=met1
                if(met2.ne.ion_num(ni)) then
                   write(stderr,*) 'ERROR: Number of atoms in input does not match meta-data for: ',trim(adjustl(ion_label(ni)))
                   stop
                end if
                ion_formula(ni)=met3
                ion_symmetry(ni)=met4
                ion_volume(ni)=met5
                ion_energy(ni)=met6
                point_count(ni)=met7
                exit
             end if
          end do
          1003 continue
          close(44)

       ! * Read meta-data from .xyz file
       else

          indx=index(line, "Lattice=")
          quote1=indx+index(line(indx:), '"')-1
          quote2=quote1+index(line(quote1+1:), '"')
          read(line(quote1+1:quote2-1),*) lattice(:)
          ion_volume(ni)=lattice(1)*(lattice(5)*lattice(9)-lattice(8)*lattice(6))+&
                         lattice(4)*(lattice(8)*lattice(3)-lattice(2)*lattice(9))+&
                         lattice(7)*(lattice(2)*lattice(6)-lattice(5)*lattice(3))

          indx=index(line, "energy=")
          read(line(indx+7:),*) ion_energy(ni)

          indx=index(line, "spacegroup=")
          read(line(indx+11:),*) ion_symmetry(ni)

          indx=index(line, "times_found=")
          read(line(indx+12:),*) point_count(ni)

       end if
       ion_volume(ni)=ion_volume(ni)/real(ion_num(ni),dp)
       ion_energy(ni)=ion_energy(ni)/real(ion_num(ni),dp)

       ! * Read SOAP vector(s)
       indx=index(line, "SOAP-")
       quote1=indx+index(line(indx:), '"')-1
       quote2=quote1+index(line(quote1+1:), '"')

       hdim=COUNT([(line(i:i),i=quote1+1,quote2-1)].eq.' ')+1
       if(ni.eq.1) allocate(soap1(hdim),soap2(hdim))

       read(line(quote1+1:quote2-1),*) soap1

       indx=quote2+index(line(quote2+1:), "SOAP-")
       if(indx.ne.quote2) then
          quote1=indx+index(line(indx:), '"')-1
          quote2=quote1+index(line(quote1+1:), '"')

          hdim=hdim*2

          read(line(quote1+1:quote2-1),*) soap2
       end if

       if(ni.eq.1) allocate(data_vec(hdim,npoints_max))

       if(indx.eq.quote2) then
          data_vec(1:hdim,ni)=soap1
       else
          data_vec(1:hdim/2,ni)=soap1
          data_vec(hdim/2+1:hdim,ni)=soap2
       end if

       ! * Deduce atomic formula from structure
       elements="aaaa"
       counts=0
       indx=0
       do i=1,ion_num(ni)
          read(stdin,*) element
          if(any(elements.eq.element)) then
             counts(indx)=counts(indx)+1
          else
             indx=indx+1
             elements(indx)=element
             counts(indx)=counts(indx)+1
          end if
       end do
       ion_formula(ni)=''

       do i=1,indx
          write(element,'(i0)') counts(i)
          ion_formula(ni)=trim(adjustl(ion_formula(ni)))//trim(adjustl(elements(i)))//trim(adjustl(element))
       end do

    end do
    1002 continue

    allocate(point_pos(ldim,npoints),point_radius(npoints),dij2(npoints,npoints))
   
    total_point_count=sum(point_count(1:npoints))

  end subroutine model_read_xyz

  subroutine model_write(cost)

    integer                   :: ni,i
    character(len=20)         :: fmt
    real(kind=dp), intent(in) :: cost

    !! * For paper: computes #contacts with distances descriptor
!    call compute_ion_nn

    write(stdout,'(i0)') reduced_point_count
    write(stdout,'(a,i0,a,f16.10)') 'dim ',ldim,' : cost', cost
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       write(fmt,*) ldim
       write(stdout,'(a2,'//trim(adjustl(fmt))//'g25.13,a240,i6,a15,a10,g25.13,g25.13,i6,g25.13)') 'H',(point_pos(i,ni),i=1,ldim),&
            trim(ion_label(ni)),ion_num(ni),trim(ion_formula(ni)),trim(ion_symmetry(ni)),ion_volume(ni),ion_energy(ni),&
            point_count(ni),max(sphere_radius,point_radius(ni))
    end do

    flush(stdout)

  end subroutine model_write

  !! * For paper: compute #contacts for distances descriptor
  subroutine compute_ion_nn
  
    integer :: ni,nj
  
    if(.not.allocated(ion_nn)) allocate(ion_nn(npoints))
  
    ion_nn_dist=2.3
  
    ion_nn=0
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=1,hdim
          if(data_vec(nj,ni).lt.ion_nn_dist) then
             ion_nn(ni)=ion_nn(ni)+1
          else
             exit
          end if
       end do
    ion_nn(ni)=ion_nn(ni)/real(ion_num(ni),dp)
    end do
  
  end subroutine compute_ion_nn

  subroutine model_init()

    integer                                  :: ni,nj,pairs

    real(kind=dp)                            :: length2(npoints)

    if(verbose) write(stderr,*) 'initialising'

    ! * Needed for UMAP
    pairs=npoints/2*(npoints-1)

    allocate(pij(npoints,npoints),sigma(npoints),simplices(pairs),head(pairs),tail(pairs))
    allocate(hcom(hdim),stddev(hdim))

    ! * Centre and normalise source data

    if(verbose) write(stderr,*) 'centring source data'

    ! * Compute centre-of-mass
    hcom=0.0_dp
    !$omp parallel do reduction(+:hcom) schedule(dynamic) 
    do ni=1,npoints
       hcom=hcom+data_vec(1:hdim,ni)
    end do
    !$omp end parallel do
    hcom=hcom/real(npoints,dp)

    ! * Move CoM to the origin - parallelise??
    ! * Compute standard deviation for each variable if normalize is set
    stddev=0.0_dp
    do ni=1,npoints
       data_vec(1:hdim,ni)=data_vec(1:hdim,ni)-hcom(:)
       stddev=stddev+data_vec(1:hdim,ni)**2 
    end do
    stddev=sqrt(stddev/real(npoints-1,dp))
!    stddev=sqrt(stddev)/real(npoints-1,dp)

    ! * Normalise such that each variable has a standard deviation of 1 across source data
    if(normalize) then
       do ni=1,hdim
          if(stddev(ni).eq.0.0_dp) cycle
          data_vec(ni,:)=data_vec(ni,:)/stddev(ni)
       end do
    end if

    ! * Calculate distances

    if(verbose) write(stderr,*) 'computing length2'

    do ni=1,npoints
       length2(ni)=dot_product(data_vec(1:hdim,ni),data_vec(1:hdim,ni))
    end do

    if(verbose) write(stderr,*) 'computing dij2'

    !$omp parallel do
    do ni=1,npoints
       do nj=1,npoints
          dij2(nj,ni)=length2(ni)+length2(nj)+1e-10_dp
       end do
       dij2(ni,ni)=0
    end do
    !$omp end parallel do

    if(verbose) write(stderr,*) 'calling dgemm'

    call dgemm('T','N',npoints,npoints,hdim,-2.0_dp,data_vec(1:hdim,1:npoints),&
               hdim,data_vec(1:hdim,1:npoints),hdim,1.0_dp,dij2,npoints)

    !$omp parallel do
    do ni=1,npoints
       dij2(ni,ni)=0.0_dp
    end do
    !$omp end parallel do

    ! * Combine similar structures

    if(sim_thresh.gt.0.0) then

       write (stderr,'(a)',advance='NO') ' checking for duplicate structures:'

       if(normalize) then
               sim_thresh=sim_thresh*real(hdim,dp)/100_dp
       else
               sim_thresh=sim_thresh*sum(stddev)/100_dp
       end if
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          do nj=ni+1,npoints
             if(point_count(nj).eq.0) cycle
             if((dij2(nj,ni).lt.(sim_thresh)**2).and.(abs(ion_energy(ni)-ion_energy(nj)).lt.energy_thresh)) then
                if(ion_energy(ni).le.ion_energy(nj)) then
                  point_count(ni)=point_count(ni)+point_count(nj)
                  point_count(nj)=0
                else
                  point_count(nj)=point_count(ni)+point_count(nj)
                  point_count(ni)=0
                  exit
                end if
             end if
          end do
       end do

       reduced_point_count=count(point_count(1:npoints).gt.0)
   
       if(reduced_point_count.eq.npoints) then
          write (stderr,'(a20)') ' no duplicates found'
       else
          write(stderr,'(a41,i15)') ' duplicates found, reduced_point_count = ',reduced_point_count
       end if

    else

       reduced_point_count=count(point_count(1:npoints).gt.0)

    end if

    ! * Compute pij using desired scheme for manifold learning

    if(use_perplexity) call perplexity_weights()
    if(use_knn)        call knn_weights()

    ! * Compute average sigma, needed for scaling random projection

    if(verbose) write (stderr,*) 'computing average sigma'

    sig_ave=sum(sigma)/real(reduced_point_count,dp)

    ! * Set czero

    if(verbose) write (stderr,*) 'setting czero'

    czero=0.0_dp
    
    if(.not.use_stress) then
      !$omp parallel do reduction(+:czero) schedule(dynamic)
      do ni=1,npoints
         if(point_count(ni).eq.0) cycle
         do nj=ni+1,npoints
            if(point_count(nj).eq.0) cycle
            if(use_kl) then
              if(pij(nj,ni).gt.epsilon(1.0_dp)) czero=czero+pij(nj,ni)*log(pij(nj,ni))*2.0_dp
            else if(use_ce) then
              if(pij(nj,ni).gt.epsilon(1.0_dp).and.pij(nj,ni).lt.(1.0_dp-epsilon(1.0_dp))) then
                czero=czero+pij(nj,ni)*log(pij(nj,ni))*2.0_dp+(1-pij(nj,ni))*log(1-pij(nj,ni))*2.0_dp
              else if(pij(nj,ni).gt.epsilon(1.0_dp).and.pij(nj,ni).gt.(1.0_dp-epsilon(1.0_dp))) then
                czero=czero+pij(nj,ni)*log(pij(nj,ni))*2.0_dp
              else if(pij(nj,ni).lt.epsilon(1.0_dp).and.pij(nj,ni).lt.(1.0_dp-epsilon(1.0_dp))) then
                czero=czero+(1-pij(nj,ni))*log(1-pij(nj,ni))*2.0_dp
              end if
            end if
         end do
      end do
      !$omp end parallel do
    end if

    deallocate(dij2)

    if(verbose) write (stderr,*) 'pij generated'

  end subroutine model_init

  subroutine perplexity_weights

    integer :: ni,nj

    ! * Set sigma according to perplexity

    if(verbose.and.use_perplexity) write(stderr,*) 'using perplexity'

    if(verbose) write (stderr,*) 'setting sigma'
    if((abs(perplexity).lt.1.0_dp).or.(abs(perplexity).gt.real(reduced_point_count-1,dp))) &
         stop 'perplexity must be greater than 1, and less than N-1'

    sigma=0.0_dp

    !$omp parallel do schedule(static)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       sigma(ni)=find_sigma(ni)
    end do
    !$omp end parallel do
    
    if (perplexity.lt.0.0_dp) then
       sigma=sum(sigma)/real(reduced_point_count,dp)
       if(verbose) write (stderr,'(a,f10.5)') ' using a mean sigma: ',sum(sigma)/real(reduced_point_count,dp)
    end if
    
    ! * Set pij with the sigma determined above

    if(verbose) write(stderr,'(a)') ' computing pij using these sigma'

    pij=0.0_dp
    !$omp parallel do schedule(dynamic)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=1,npoints
          if((ni.eq.nj).or.(point_count(nj).eq.0)) cycle
          pij(nj,ni)=exp(-dij2(nj,ni)/sigma(ni)**2/2.0_dp)
       end do
       pij(:,ni)=pij(:,ni)/sum(pij(:,ni))
    end do
    !$omp end parallel do

    ! * Symmetrise

    if(use_sgd) then

       n_simplices=0
       do ni=1,npoints
         if(point_count(ni).eq.0) cycle
         do nj=ni+1,npoints
           if(point_count(nj).eq.0) cycle
           pij(ni,nj)=(pij(ni,nj)+pij(nj,ni))/real(2*reduced_point_count,dp)
           if(pij(ni,nj).lt.1e-6) pij(ni,nj)=0.0_dp
           pij(nj,ni)=pij(ni,nj)
           if(pij(ni,nj).gt.0) then
             n_simplices=n_simplices+1
             simplices(n_simplices)=pij(ni,nj)
             head(n_simplices)=ni
             tail(n_simplices)=nj
           end if
         end do
       end do

       if(verbose) then
          write(stderr,'(a,i10)') '   number of pairs  :',pairs
          write(stderr,'(a,i10)') '   non-zero weights :',n_simplices
       end if

       if(sexag.le.0.0_dp) sexag=(npoints*(npoints-1)/2)/n_simplices/nneg

    else

       !$omp parallel do schedule(dynamic)
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          do nj=ni+1,npoints
             if(point_count(nj).eq.0) cycle
             pij(ni,nj)=(pij(ni,nj)+pij(nj,ni))/real(2*reduced_point_count,dp)
             pij(nj,ni)=pij(ni,nj)
          end do
       end do
       !$omp end parallel do

    end if

  end subroutine perplexity_weights

  subroutine knn_weights()

    integer                                  :: ni,nj,k,first_nn_indx

    real(kind=dp)                            :: first_nn_dist
    real(kind=dp), allocatable, dimension(:) :: temp_vec1,temp_vec2

    ! * Find k-nn to each ion

    if(verbose.and.use_knn) write(stderr,*) 'using K-NN'

    if(verbose) write (stderr,*) 'finding k-NNs to each ion'

    if(knn.gt.(npoints-1).or.knn.lt.0) &
       stop 'number of NN must be an integer greater than 0, and less than N-1'

    allocate(nn_vec_indx(npoints,knn),nn_vec_dist(npoints,knn),temp_vec1(npoints),temp_vec2(knn))

    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       k=1
       temp_vec1(:)=dij2(ni,:)
       do nj=1,npoints
         if(point_count(nj).eq.0) temp_vec1(nj)=0 
       end do
       nn_vec_indx(ni,k)=ni
       nn_vec_dist(ni,k)=0.0_dp
       do
         k=k+1
         nn_vec_indx(ni,k)=minloc(temp_vec1,dim=1,mask=temp_vec1.gt.0)
         nn_vec_dist(ni,k)=sqrt(minval(temp_vec1,mask=temp_vec1.gt.0))
         temp_vec1(nn_vec_indx(ni,k))=0
         if(k.eq.knn) exit
       end do
    end do

    ! * Set local distance metric, sigma, for each ion, according to its k-NNs

    if(verbose) write (stderr,*) 'setting sigmas'

    !$omp parallel do schedule(static)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       sigma(ni)=find_sigma(ni)
    end do
    !$omp end parallel do

    ! * Set pij with the sigmas determined above

    write(stderr,'(a)') ' computing pij using these sigma'

    pij(:,:)=0.0_dp

    do ni=1,npoints
       temp_vec2=nn_vec_dist(ni,:)
       first_nn_indx=minloc(temp_vec2,dim=1,mask=temp_vec2.gt.0)
       first_nn_dist=minval(temp_vec2,mask=temp_vec2.gt.0)
       if(point_count(ni).eq.0) cycle
       pij(ni,ni)=0.0_dp
       pij(ni,nn_vec_indx(ni,first_nn_indx))=1.0_dp
       do nj=2,knn
         if(nj.ne.first_nn_indx) then
            pij(ni,nn_vec_indx(ni,nj))=exp(-(nn_vec_dist(ni,nj)-first_nn_dist)/sigma(ni))
         end if
       end do
    end do

    ! * Take union of probabilities for each pair of ions

    n_simplices=0
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       do nj=ni+1,npoints
          if(point_count(nj).eq.0) cycle
          pij(ni,nj)=pij(ni,nj)+pij(nj,ni)-pij(ni,nj)*pij(nj,ni)
          pij(nj,ni)=pij(ni,nj)
          if(pij(ni,nj).gt.0) then
             n_simplices=n_simplices+1
             simplices(n_simplices)=pij(ni,nj)
             head(n_simplices)=ni
             tail(n_simplices)=nj
          end if
       end do
    end do

    if(verbose) then
       write(stderr,'(a,i10)') '   number of pairs  :',pairs
       write(stderr,'(a,i10)') '   non-zero weights :',n_simplices
    end if

    if(sexag.le.0.0_dp) sexag=pairs/n_simplices/nneg

  end subroutine knn_weights

  function find_sigma(i) result(sss)

    integer, intent(in) :: i
    integer             :: n

    real(kind=dp)       :: s,p,ss,pp,sss,ppp,factor

    ! * Bracket

    s=1.0_dp

    factor=gr

    if(use_perplexity) then
       p=deltaperp(i,s)
    else if(use_knn) then
       p=deltalog2k(i,s)
    end if    

    if(p.gt.0.0_dp) factor=1.0_dp/factor

    ss=s*factor

    if(use_perplexity) then
       pp=deltaperp(i,ss)
    else if(use_knn) then
       pp=deltalog2k(i,ss)
    end if

    n=0
    do while (p*pp.gt.0.0_dp)
       n=n+1
       s=ss
       p=pp
       ss=ss*factor
       if(use_perplexity) then
          pp=deltaperp(i,ss)
       else if(use_knn) then
          pp=deltalog2k(i,ss)
       end if
       if(n.gt.100) then
          write (stderr,*) s,p,ss,pp
          stop 'bracket failed'
       end if
    end do

    ! * Bisect

    ppp=huge(1.0_dp)
    do while(abs(s-ss).gt.1e-13_dp.and.abs(ppp).gt.1e-13)

       sss=(s+ss)/2.0_dp; if((sss.eq.s).or.(sss.eq.ss)) exit
       if(use_perplexity) then
          ppp=deltaperp(i,sss)
       else if(use_knn) then
          ppp=deltalog2k(i,sss)
       end if

       if(p*ppp.lt.0.0_dp) then
          ss=sss
          pp=ppp
       else
          s=sss
          p=ppp
       end if

    end do

  end function find_sigma
  
  function deltaperp(i,sig_update)

    integer,       intent(in) :: i
    integer                   :: j

    real(kind=dp), intent(in) :: sig_update    
    real(kind=dp)             :: deltaperp,pi(npoints)

    pi=0.0_dp
    do j=1,npoints
       if((i.eq.j).or.(point_count(j).eq.0)) cycle
       pi(j)=exp(-dij2(j,i)/sig_update**2/2.0_dp)
    end do

    pi(:)=pi(:)/sum(pi)

    deltaperp=0.0_dp
    do j=1,npoints
       if((i.eq.j).or.(point_count(j).eq.0)) cycle
       if(pi(j).gt.tiny(1.0_dp)) deltaperp=deltaperp-pi(j)*log(pi(j))/log(2.0_dp)
    end do
    deltaperp=2.0_dp**deltaperp-abs(perplexity)

  end function deltaperp

  function deltalog2k(i,sig_update)

    integer,       intent(in) :: i
    integer                   :: k

    real(kind=dp), intent(in) :: sig_update
    real(kind=dp)             :: first_nn,dist,temp_vec(knn),exp_sum,log2_knn,deltalog2k

    log2_knn=log(real(knn))/log(2.0_dp)

    exp_sum=0.0_dp
    temp_vec=nn_vec_dist(i,:)
    first_nn=minval(temp_vec,mask=temp_vec.gt.0)
    do k=2,knn
       dist=nn_vec_dist(i,k)-first_nn
       if(dist.gt.0.0_dp) then
          exp_sum=exp_sum+exp(-dist/sig_update)
       else
          exp_sum=exp_sum+1.0_dp
       end if
    end do
    deltalog2k=exp_sum-log2_knn

  end function deltalog2k

  subroutine model_uniform_packing()

    integer       :: ni,nj

    real(kind=dp) :: vec(ldim),dist,overlap,grad(ldim,npoints),vt(ldim,npoints),gg,com(ldim)
    real(kind=dp) :: radius_separate,radius_pack

    if(verbose) write (stderr,*) 'generating starting points'

!    radius_pack=(real(npoints,dp)/packing_fraction)**(1.0_dp/real(ldim,dp))*sphere_radius
    radius_pack=(real(npoints,dp)/packing_fraction)**(1.0_dp/real(ldim,dp))


!    radius_separate=sphere_radius
    radius_separate=0.0_dp
    point_radius=0.0_dp
    
    overlap=huge(1.0_dp)
    do while(overlap.gt.1e-6_dp)

!write(stderr,*) 1

       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          call random_point_in_hypersphere(point_pos(:,ni),radius_pack-point_radius(ni))
!write(stderr,*) point_pos(:,ni)
       end do
      
       gg=huge(1.0_dp)
       vt=0.0_dp
       do while(gg.gt.1e-13_dp)
!write(stderr,*) gg,'too big'
          overlap=0.0_dp
          grad=0.0_dp
          !$omp parallel do private(vec,dist) reduction(+:overlap,grad) schedule(dynamic)
          do ni=1,npoints
             if(point_count(ni).eq.0) cycle

             ! * Other atoms
             do nj=ni+1,npoints
                if(point_count(nj).eq.0) cycle

                vec=point_pos(:,ni)-point_pos(:,nj)
                dist=norm2(vec)
                if(dist.lt.2.0_dp*radius_separate) then
                   vec=vec/dist
                   dist=radius_separate-dist/2.0_dp
                   grad(:,ni)=grad(:,ni)+vec*dist
                   grad(:,nj)=grad(:,nj)-vec*dist
                   overlap=overlap+dist**2/2
                end if
             end do
!write(stderr,*) 'oa',grad

!             ! * Boundary
!             vec=point_pos(:,ni)
!             dist=norm2(vec)
!             vec=vec/dist
!             if(dist.gt.(radius_pack-point_radius(ni))) then
!                dist=radius_pack-radius_separate-dist
!                grad(:,ni)=grad(:,ni)+vec*dist
!                overlap=overlap+dist**2/2
!             end if
!write(stderr,*) 'boundary',grad

          end do
          !$omp end parallel do

          gg=0.0_dp
          do ni=1,npoints
             if(point_count(ni).eq.0) cycle
             gg=gg+dot_product(grad(:,ni),grad(:,ni))
          end do

          vt=0.9_dp*vt+0.5_dp*grad

          point_pos=point_pos+vt

       end do

    end do

    if(verbose) write (stderr,'(a,f8.3,a,f6.3)') ' radius_pack :',radius_pack,' radius_separate :',radius_separate

    ! * Compute centre-of-mass of initial embedding
    com=0.0_dp
    !$omp parallel do reduction(+:com) schedule(dynamic)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       com=com+point_pos(:,ni)
    end do
    !$omp end parallel do
    com=com/real(reduced_point_count,dp)

    ! * Move CoM to the origin
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       point_pos(:,ni)=point_pos(:,ni)-com(:)
    end do
    
  end subroutine model_uniform_packing

  subroutine model_random_projection

    integer       :: i,j,k,ni

    real(kind=dp) :: random_vecs(ldim,hdim),coeff,com(ldim)!,hcom(hdim)

    ! * Compute (ldim x hdim) matrix of random orthonormal vectors in high-D space
    random_vecs=0.0_dp
    do i=1,ldim
       call random_add_noise(random_vecs(i,:),1.0_dp)
       random_vecs(i,:)=random_vecs(i,:)/norm2(random_vecs(i,:))
!!       do j=1,hdim
!!          if(random_vecs(i,j).lt.0.0_dp) random_vecs(i,j)=-1.0_dp*random_vecs(i,j)
!!       end do
       if(i.gt.1) then
          do k=1,i-1
             coeff=dot_product(random_vecs(i-k,:),random_vecs(i,:))
             do j=1,hdim
                random_vecs(i,j)=random_vecs(i,j)-coeff*random_vecs(i-k,j)
             end do
          end do
       end if
       random_vecs(i,:)=random_vecs(i,:)/norm2(random_vecs(i,:))
    end do

    ! * Compute random projection of source date into low-D space (embedded postions scaled by average sigma and compression factor)
    sig_ave=sum(sigma)/real(reduced_point_count,dp)
    point_pos=matmul(random_vecs,data_vec(1:hdim,1:npoints))/sig_ave/projection_compression

    ! * Compute centre-of-mass of initial embedding
    com=0.0_dp
    !$omp parallel do reduction(+:com) schedule(dynamic) 
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       com=com+point_pos(:,ni)
    end do
    !$omp end parallel do
    com=com/real(reduced_point_count,dp)

    ! * Move CoM to the origin and compute the initial sphere radius (if not set)
    if(sphere_radius.le.0.0_dp) then
    
       sphere_radius=0.0_dp
       !$omp parallel do reduction(+:point_pos,sphere_radius) schedule(dynamic)
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          point_pos(:,ni)=point_pos(:,ni)-com(:)
          sphere_radius=sphere_radius+dot_product(point_pos(:,ni),point_pos(:,ni))
       end do
       !$omp end parallel do
       sphere_radius=sqrt(sphere_radius)/real(reduced_point_count)/10.0_dp*projection_compression 
       recalc_spheres=.true.

    else

       !$omp parallel do reduction(+:point_pos) schedule(dynamic)
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          point_pos(:,ni)=point_pos(:,ni)-com(:)
       end do
       !$omp end parallel do

    end if

    ! * Set initial ion radii
    point_radius=0.0_dp

    !$omp parallel do schedule(static)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       point_radius(ni)=sphere_radius*real(point_count(ni),dp)**(1.0_dp/real(ldim,dp))
    end do
    !$omp end parallel do

  end subroutine model_random_projection

  subroutine model_pca

    integer                                    :: ni,nj,point_counter,e_num,info,iwork(5*hdim),ifail(hdim)

    real(kind=dp)                              :: com(ldim)
    real(kind=dp)                              :: e_val(hdim),vl,vu,lwork
    real(kind=dp), allocatable, dimension(:)   :: work
    real(kind=dp), allocatable, dimension(:,:) :: data_red,cov,e_vec

    allocate(cov(hdim,hdim),e_vec(hdim,hdim))

    ! * Compute covariance matrix
    if(reduced_point_count.ne.npoints) then
       allocate(data_red(hdim,reduced_point_count))
       point_counter=0
       do ni=1,npoints
          if(point_count(ni).eq.0) cycle
          point_counter=point_counter+1
          data_red(:,point_counter)=data_vec(:,ni)
       end do
       cov=matmul(data_red(1:hdim,1:reduced_point_count),transpose(data_red(1:hdim,1:reduced_point_count)))
       deallocate(data_red)
    else
       cov=matmul(data_vec(1:hdim,1:reduced_point_count),transpose(data_vec(1:hdim,1:reduced_point_count)))
    end if

    ! * Eigendecomposition of covariance matrix

       ! * Determine optimal lwork and liwork
       ifail=0
       call dsyevx('V','I','U',hdim,cov,hdim,vl,vu,hdim-ldim+1,hdim,epsilon(1.0),&
               e_num,e_val,e_vec,hdim,lwork,-1,iwork,ifail,info)
       allocate(work(int(lwork)))

       ! * Do eigendecomposition
       call dsyevx('V','I','U',hdim,cov,hdim,vl,vu,hdim-ldim+1,hdim,epsilon(1.0),&
               e_num,e_val,e_vec,hdim,work,int(lwork),iwork,ifail,info)
       deallocate(cov)

    ! * Project data onto principal components
    point_pos=matmul(transpose(e_vec(1:hdim,ldim:1:-1)),data_vec(1:hdim,1:npoints))/sig_ave/projection_compression
    deallocate(e_vec)

  end subroutine model_pca

  subroutine model_resume()

    integer :: ni,stat

    character(len=4) element

    nexag=0

    open(unit=88,file=resume_file,iostat=stat,status='old')
    if(stat.ne.0) then
       write(stderr,*) 'ERROR: resume file does not exist'
       stop
    end if
    read(88,*)
    read(88,*)
    do ni=1,npoints
       if(point_count(ni).eq.0) cycle
       read(88,*,end=1004) element,point_pos(:,ni)
    end do
    1004 continue
    close(88)

  end subroutine model_resume

  subroutine model_report()
    
    ! * Header
    write (stderr,*) 'model parameters'
    write (stderr,*)

    ! * Number of data points
    if(total_point_count.ne.npoints) write (stderr,'(a25,a3,i8)') 'total_point_count',' : ', total_point_count
    write (stderr,'(a25,a3,i8)') 'npoints',' : ', npoints
    if(reduced_point_count.ne.npoints) write (stderr,'(a25,a3,i8)') 'reduced_point_count',' : ', reduced_point_count
    write (stderr,'(a25,a3,i8)') 'hdim',' : ', hdim

    ! * Minimum sphere radius
    if(recalc_spheres) write (stderr,'(a25,a3,f14.8)') 'initial sphere_radius',' : ', sphere_radius

    write (stderr,*)

  end subroutine model_report

  ! -------------------------------------------------------------------
  
end module model
