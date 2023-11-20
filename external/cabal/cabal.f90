!==================================================================================!
!                                     cabal                                        !
!==================================================================================!
!                                                                                  !
! This is a modified version of a file from the AIRSS structure prediction package.!
!                                                                                  !
! AIRSS is free software; you can redistribute it and/or modify it under the terms !
! of the GNU General Public License version 2 as published by the Free Software    !
! Foundation.                                                                      !
!                                                                                  !
! This program is distributed in the hope that it will be useful, but WITHOUT ANY  !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  !
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.        !           
!                                                                                  !
! You should have received a copy of the GNU General Public License along with this!
! program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street,!                   
! Fifth Floor, Boston, MA  02110-1301, USA.                                        !
!                                                                                  !
!----------------------------------------------------------------------------------!
! This program converts structure formats                                          !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2020                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program cabal

  use constants
  use spglib_f08

  implicit none

  integer, parameter :: max_words=120
  integer, parameter :: max_species=20

  integer :: num_lines=1000000 ! ** STDIN can only be read once, so need to estimate a maximum input length

  integer :: num_ions,num_ions_orig,num_species

  integer :: mp_kgrid(3)

  integer, allocatable, dimension(:) :: num_words,ion_index

  real(kind=dp)                                   :: lattice_abc(6)=0.0_dp,lattice_car(3,3)=0.0_dp,lattice_rec(3,3)=0.0_dp
  real(kind=dp)                                   :: kspacing=1.0_dp
  real(kind=dp)                                   :: lattice_volume=0.0_dp
  real(kind=dp)                                   :: pad=5.0_dp,tol=1e-1_Dp
  real(kind=dp),      allocatable, dimension(:)   :: ion_occ,ion_spin
  real(kind=dp),      allocatable, dimension(:,:) :: ion_fractional,ion_absolute

  character(len=10)                                       :: inext,outext,species_names(max_species)=''
  character(len=10),          allocatable, dimension(:)   :: ion_names
  character(len=max_words*20), allocatable, dimension(:)  :: buff,rembuff
  character(len=max_words*2), allocatable, dimension(:,:) :: buff_words

  character(len=40) :: structure_name=''
  real(kind=dp)     :: pressure=0.0_dp
  real(kind=dp)     :: enthalpy=0.0_dp
  real(kind=dp)     :: total_spin=0.0_dp
  real(kind=dp)     :: total_modspin=0.0_dp
  character(len=20) :: spacegroup='n/a'
  integer           :: ncopies=1
  
  logical :: x2x
  logical :: kpts=.false.

  type(SpglibDataset) :: dset

  integer :: ni
  
  real(kind=dp) :: lattice_trans(3,3)=0.0_dp
  
  ! * =======================================================================

  call get_arguments()

  call read_buff()

  call read_input()

  call consolidate()

  if(x2x) then
     
     if(abs(tol).gt.delta) then

        ! ** spglib uses the transpose of lattice_car
        
        lattice_trans=transpose(lattice_car)
        
        dset=spg_get_dataset(lattice_trans,ion_fractional,ion_index,num_ions,abs(tol))

        ! ** Get the international symbol for the space group based on the provided tolerance
        
        if(dset%spacegroup_number/=0) then
           spacegroup=deunder(trim(dset%international_symbol))
        else
           spacegroup='P1'
        end if

        ! ** Refine to conventional or primitive cell
        
        if(tol.gt.0.0_dp) then
           num_ions=spg_refine_cell(lattice_trans,ion_fractional,ion_index,num_ions,abs(tol))
        else
           num_ions=spg_find_primitive(lattice_trans,ion_fractional,ion_index,num_ions,abs(tol))
        end if

        ! ** Reset the lattice

        lattice_car=transpose(lattice_trans)
        call car2abc(lattice_car,lattice_abc)

        ! ** Reset the ion names

        do ni=1,num_ions
           ion_names(ni)=species_names(ion_index(ni))
        end do
        
     else
        
        call niggli() ! ** Replace with spglib
        
     end if
     
  end if

  call write_output()

contains

  subroutine get_arguments()

    integer :: iargc

    character(len=20) :: ctemp

    if((iargc().lt.2).or.(iargc().gt.3)) then
       
       write (stderr,'(a)') 'Usage: cabal in out [x] < seed.in > seed.out'
       write (stderr,'(a)')
       write (stderr,'(a)') '  in==out  0   : Niggli reduce'
       write (stderr,'(a)') '  in==out +tol : Refine to conventional cell'
       write (stderr,'(a)') '  in==out -tol : Refine to primitive cell'
       write (stderr,'(a)')
       write (stderr,'(a)') '  supports castep+,cell,shx,res,gulp*,cif*,psi4*,xtl,xyz(e),config+,poscar,qe,cp2k,conf'
       write (stderr,'(a)') '  *output only +input only'
       write (stderr,'(a)')
       write (stderr,'(a)') '  x : for xyz input, pad unit cell with x'
       stop
    end if

    call getarg(1,inext)
    call getarg(2,outext)

    x2x=inext.eq.outext

    if(iargc().eq.3) then
       call getarg(3,ctemp)
       if(.not.x2x) then
          read (ctemp,*) pad
       else
          read (ctemp,*) tol          
       end if
    end if

  end subroutine get_arguments

  subroutine read_buff()

    integer :: n,nw,indx

    character(len=max_words*2) :: line

    allocate(buff(num_lines))
    n=0
    do 
       n=n+1
       if(n.gt.num_lines) stop 'cabal: increase num_lines'
       read(5,'(a)',end=99,err=100) buff(n)
    end do
100 stop 'cabal: problem reading input data'
99  continue

    num_lines=n-1

    allocate(buff_words(max_words,num_lines),num_words(num_lines))
    
    do n=1,num_lines

       nw=0
       line=trim(adjustl(detab(buff(n))))
       do while (len_trim(line).gt.0)
          nw=nw+1
          indx=index(line,' ')
          buff_words(nw,n)=trim(line(:indx))
          line=adjustl(line(indx+1:))
       end do
       num_words(n)=nw

    end do

  end subroutine read_buff

  subroutine read_input()

    select case (inext)
    case('castep')
       call read_castep()
    case('cell')
       call read_cell()
    case('res','shx')
       call read_res()
    case('xyz')
       call read_xyz()
    case('xyze')
       call read_xyze()
    case('gulp')
       call read_gulp()
    case('xtl')
       call read_xtl()
    case('cif')
       call read_cif()
    case('psi4')
       call read_psi4()
    case('config')
       call read_config()
    case('poscar')
       call read_poscar()
    case('qe')
       call read_qe()
    case('cp2k')
       call read_cp2k()
    case('conf')
       call read_conf() 
    case default
       stop 'cabal: input format not recognised'
    end select

  end subroutine read_input

  subroutine write_output()

    select case (outext)
    case('cell')
       call write_cell()
    case('castep')
       call write_castep()
    case('res','shx')
       call write_res()
    case('xyz')
       call write_xyz()
    case('xyze')
       call write_xyze()
    case('gulp')
       call write_gulp()
    case('xtl')
       call write_xtl()
    case('cif')
       call write_cif()
    case('psi4')
       call write_psi4()
    case('config')
       call write_config()
    case('poscar')
       call write_poscar()
    case('qe')
       call write_qe()
    case('cp2k')
       call write_cp2k()
    case('conf')
       call write_conf()
    case('json')
       call write_json()
    case default
       stop 'cabal: output format not recognised'
    end select

  end subroutine write_output

  subroutine consolidate()

    integer :: n,ni

    if(all(lattice_car.eq.0.0_dp)) call abc2car(lattice_abc,lattice_car)

    if(all(lattice_abc.eq.0.0_dp)) call car2abc(lattice_car,lattice_abc)

    call car2rec(lattice_car,lattice_rec)

    if(all(ion_fractional.eq.0.0_dp)) then
       do ni=1,num_ions
          ion_fractional(:,ni)=matmul(lattice_rec,ion_absolute(:,ni))
       end do
    end if

    if(all(ion_absolute.eq.0.0_dp)) then
       do ni=1,num_ions
          ion_absolute(:,ni)=matmul(lattice_car,ion_fractional(:,ni))
       end do
    end if

    if(all(ion_index.eq.0)) then
       num_species=0
       species_names=''
       do ni=1,num_ions
          if(.not.any(species_names==ion_names(ni))) then
             num_species=num_species+1
             species_names(num_species)=ion_names(ni)
          end if
       end do

       do ni=1,num_ions
          do n=1,num_species
             if(ion_names(ni).eq.species_names(n)) exit
          end do
          ion_index(ni)=n
       end do

    end if

    lattice_volume=car_volume(lattice_car)

    if(kpts) then
       ! Compute k-point grid
       ! lattice_rec has reciprocal lattice vectors as rows (without 2*pi prefactor)
       mp_kgrid(1)=ceiling(sqrt(lattice_rec(1,1)**2 + lattice_rec(1,2)**2 + lattice_rec(1,3)**2)/kspacing)
       mp_kgrid(2)=ceiling(sqrt(lattice_rec(2,1)**2 + lattice_rec(2,2)**2 + lattice_rec(2,3)**2)/kspacing)
       mp_kgrid(3)=ceiling(sqrt(lattice_rec(3,1)**2 + lattice_rec(3,2)**2 + lattice_rec(3,3)**2)/kspacing)
    end if

    num_ions_orig=num_ions
    
  end subroutine consolidate

  subroutine niggli()

    real(kind=dp) :: lattice_car0(3,3)
    integer       :: Cmat(3,3)

    integer :: ni

    interface
       subroutine niggli_reduce(Lm,Cm)
         use constants
         real(kind=dp), dimension(3,3), intent(inout) :: Lm
         integer,       dimension(3,3), intent(out)   :: Cm
       end subroutine niggli_reduce
    end interface

    lattice_car0=lattice_car
    call niggli_reduce(lattice_car,Cmat)
    call car2rec(lattice_car,lattice_rec)    
    call car2abc(lattice_car,lattice_abc)

    do ni=1,num_ions
       ion_fractional(:,ni)=matmul(lattice_rec,matmul(lattice_car0,ion_fractional(:,ni)))
       ion_absolute(:,ni)=matmul(lattice_car,ion_fractional(:,ni))
    end do

  end subroutine niggli


  ! * =======================================================================

  ! ** Castep output file

  subroutine read_castep()

    integer :: n,i,ni

    do n=1,num_lines
       if(index(buff(n),'Real Lattice(A)').gt.0) then
          do i=1,3
             read(buff_words(1:3,n+i),*) lattice_car(1:3,i)
          end do
       end if
    end do

    do n=1,num_lines
       if(index(buff(n),'Total number of ions in cell =').gt.0) read(buff_words(8,n),*) num_ions
    end do

    call init_ions()

    do n=1,num_lines
       if(index(buff(n),'Number           u          v          w').gt.0) then
          do ni=1,num_ions
             read(buff_words(2,n+ni+1),*) ion_names(ni)
             ion_names(ni)=trim(strip(ion_names(ni)))
             read(buff_words(4:6,n+ni+1),*) ion_fractional(1:3,ni)
          end do
       end if
       if(index(buff(n),'Species          Ion Spin      s      p      d      f').gt.0) then
          do ni=1,num_ions
             read(buff_words(10,n+2*ni),*) ion_spin(ni)
          end do
       end if
    end do

  end subroutine read_castep

  subroutine write_castep()
    stop 'write_castep: not implemented'
  end subroutine write_castep

  ! ** Castep cell file format

  subroutine read_cell()

    integer :: n,i,ni,nstart,nfinish,indx

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_CART'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_CART'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.eq.4) then
          do i=1,3
             read(buff_words(1:3,nstart+i),*) lattice_car(1:3,i)
          end do
       else if(nfinish-nstart.eq.5) then
          if(up(buff_words(1,nstart+1)).ne.'ANG') stop 'cabal: LATTICE_CART units not recognised'
          do i=1,3
             read(buff_words(1:3,nstart+i+1),*) lattice_car(1:3,i)
          end do
       else
          stop 'cabal: problem reading LATTICE_CART in cell data'
       end if
    end if

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_ABC'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_ABC'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.ne.3) stop 'cabal: problem reading LATTICE_ABC in cell data'
       read(buff_words(1:3,nstart+1),*) lattice_abc(1:3)
       read(buff_words(1:3,nstart+2),*) lattice_abc(4:6)
    end if

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_FRAC'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_FRAC'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.lt.2) stop 'cabal: problem reading POSITIONS_FRAC in cell data'
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          ion_names(ni)=trim(buff_words(1,nstart+ni))
          read(buff_words(2:4,nstart+ni),*) ion_fractional(1:3,ni)
          indx=index(buff(nstart+ni),'SPIN=')
          if(indx.gt.0) then
             read(buff(nstart+ni)(indx+5:),*) ion_spin(ni)
          end if
       end do
    end if

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_ABS'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_ABS'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.lt.2) stop 'cabal: problem reading POSITIONS_ABS in cell data'
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          ion_names(ni)=trim(buff_words(1,nstart+ni))
          read(buff_words(2:4,nstart+ni),*) ion_absolute(1:3,ni)
          indx=index(buff(nstart+ni),'SPIN=')
          if(indx.gt.0) then
             read(buff(nstart+ni)(indx+5:),*) ion_spin(ni)
          end if
       end do
    end if

    ! K-points
    kpts=.false.
    do n=1,num_lines
       if(up(buff_words(1,n)).eq.'KPOINTS_MP_SPACING') then
          kpts=.true.
          read(buff_words(2,n),*) kspacing
          if(abs(kspacing).lt.epsilon(1.0_dp)) then
             write (stdout,'(a)') 'cabal: read_cell'
             write (stdout,'(a)') 'KPOINTS_MP_SPACING is apparently 0 (zero) which is not allowed'
             stop
          end if
          exit
       end if
    end do

    return

  end subroutine read_cell

  subroutine write_cell()

    integer :: ni

    write (stdout,'(a)') '%BLOCK LATTICE_CART'
    write (stdout,'(3f20.14)') lattice_car
    write (stdout,'(a)') '%ENDBLOCK LATTICE_CART'
    write (stdout,*)
    write (stdout,'(a)') '%BLOCK POSITIONS_FRAC'
    if(any(ion_spin.ne.0.0_dp)) then
       do ni=1,num_ions
          write (stdout,'(a,3f20.14," SPIN=",f7.3)') trim(adjustl(species_names(ion_index(ni)))),ion_fractional(:,ni),ion_spin(ni)
       end do
    else
       do ni=1,num_ions
          write (stdout,'(a,3f20.14)') trim(adjustl(species_names(ion_index(ni)))),ion_fractional(:,ni)
       end do
    end if

    write (stdout,'(a)') '%ENDBLOCK POSITIONS_FRAC'

  end subroutine write_cell

  ! ** SHLX res (result) file format

  subroutine read_res()

    integer :: n,i,ns,ni,nrem

    if(.not.(buff_words(1,1).eq.'TITL')) stop 'cabal: first line of res/shx data should start with TITL'
    if(.not.(buff_words(1,num_lines).eq.'END')) stop 'cabal: last line of res/shx data should start be END'
    
    if(num_words(1).ne.12) stop 'problem reading TITL line'
    
    structure_name=trim(buff_words(2,1))
    read(buff_words(3,1),*) pressure
    read(buff_words(4,1),*) lattice_volume
    read(buff_words(5,1),*) enthalpy
    read(buff_words(6,1),*) total_spin
    read(buff_words(7,1),*) total_modspin
    spacegroup=deunder(trim(buff_words(9,1)))
    spacegroup=trim(spacegroup(2:len_trim(spacegroup)-1))
    read(buff_words(12,1),*) ncopies

    nrem=count(buff_words(1,:).eq.'REM ')

    allocate(rembuff(nrem))

    n=0
    do i=1,num_lines
       if(buff_words(1,i).eq.'REM ') then
          n=n+1
          rembuff(n)=buff(i)
       end if
    end do
    
    do n=1,num_lines
       if(buff_words(1,n).eq.'CELL') then
          do i=1,6
             read(buff_words(i+2,n),*) lattice_abc(i)
          end do
       end if
       if(buff_words(1,n).eq.'SFAC') then
          num_species=num_words(n)-1
          do ns=1,num_species
             read(buff_words(ns+1,n),*) species_names(ns)
          end do
          exit
       end if
    end do

    num_ions=num_lines-n-1

    call init_ions()

    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,n+ni))
       read(buff_words(2,n+ni),*) ion_index(ni)
       read(buff_words(3:5,n+ni),*) ion_fractional(:,ni)
       read(buff_words(6,n+ni),*) ion_occ(ni)
       if(num_words(n+ni).gt.6) read(buff_words(7,n+ni),*) ion_spin(ni)
    end do

    return

  end subroutine read_res

  subroutine write_res()

    integer :: n

    real(kind=dp) :: factor

    character(len=40) :: fmt,ctemp,ctemp2,ctemp3

    factor=real(num_ions,dp)/real(num_ions_orig,dp)
    
    write(ctemp,'(f15.6)') lattice_volume
    write(ctemp2,*) num_ions
    write(ctemp3,*) ncopies

    if(structure_name.eq.'') then
       structure_name='cabal-'//'in-'//'out'
    end if
    
    write (stdout,'(a,a,3f20.12,2f6.2,1x,a,a,i5)') 'TITL ',trim(structure_name),pressure,lattice_volume*factor,&
         enthalpy*factor,total_spin*factor,total_modspin*factor,trim(adjustl(ctemp2)),&
         ' ('//trim(spacegroup)//') n - '//trim(adjustl(ctemp3))

    if(allocated(rembuff)) then
       do n=1,size(rembuff)
          write (stdout,'(a)') trim(rembuff(n))
       end do
    end if
    
    write (stdout,'(a,6f11.5)') 'CELL 1.54180',lattice_abc
    write (stdout,'(a)') 'LATT -1'
    write (ctemp,*) num_species
    write (fmt,*) "(a,"//trim(adjustl(ctemp))//"a3)"
    write (stdout,trim(adjustl(fmt))) 'SFAC ',species_names(1:num_species)
    if(any(ion_spin.ne.0.0_dp)) then
       do n=1,num_ions
          write (stdout,'(a4,i4,3f17.13,f4.1,f7.2)') &
               species_names(ion_index(n)),ion_index(n),ion_fractional(:,n),ion_occ(n),ion_spin(n)
       end do
    else
       do n=1,num_ions
          write (stdout,'(a4,i4,3f17.13,f4.1)') &
               species_names(ion_index(n)),ion_index(n),ion_fractional(:,n),ion_occ(n)
       end do
    end if
    write (stdout,'(a)') 'END'

  end subroutine write_res

  ! ** XYZ file format

  subroutine read_xyz()

    integer :: ni
    real(kind=dp) :: lat(3),com(3),moi(3,3),eig(3),dii

    read(buff_words(1,1),*) num_ions

    call init_ions()

    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,ni+2))
       ion_absolute(:,ni)=0.0_dp
       read(buff_words(2:2,ni+2),*,err=101) ion_absolute(1,ni)
101    read(buff_words(3:3,ni+2),*,err=102) ion_absolute(2,ni)
102    read(buff_words(4:4,ni+2),*,err=103) ion_absolute(3,ni)
103    continue
    end do

    com=sum(ion_absolute,2)/real(num_ions)

    moi=0.0_dp
    do ni=1,num_ions
       dii=dot_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
       moi(1,1)=moi(1,1)+dii
       moi(2,2)=moi(2,2)+dii
       moi(3,3)=moi(3,3)+dii
       moi=moi-outer_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
    end do

    call diag(moi,eig)

    do ni=1,num_ions
       ion_absolute(:,ni)=matmul(moi,ion_absolute(:,ni)-com)
    end do

    lat(1)=maxval(ion_absolute(1,:))-minval(ion_absolute(1,:))+pad
    lat(2)=maxval(ion_absolute(2,:))-minval(ion_absolute(2,:))+pad
    lat(3)=maxval(ion_absolute(3,:))-minval(ion_absolute(3,:))+pad

    lattice_abc(1:6) =(/lat(1),lat(2),lat(3),90.0_dp,90.0_dp,90.0_dp/)

    do ni=1,num_ions
       ion_absolute(:,ni)=ion_absolute(:,ni)+lat/2.0_dp
    end do

  end subroutine read_xyz

  subroutine write_xyz()

    integer :: ni
    character(len=10) :: ctemp
    real(kind=dp) :: com(3),moi(3,3),eig(3),dii

    com=sum(ion_absolute,2)/real(num_ions)

    moi=0.0_dp
    do ni=1,num_ions
       dii=dot_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
       moi(1,1)=moi(1,1)+dii
       moi(2,2)=moi(2,2)+dii
       moi(3,3)=moi(3,3)+dii
       moi=moi-outer_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
    end do

!!$    write (stderr,'(3f10.5)') moi
!!$    write (stderr,*)

    call diag(moi,eig)

!!$    write (stderr,'(3f10.5)') eig
!!$    write (stderr,*)
!!$    
!!$    write (stderr,'(3f10.5)') moi
!!$    write (stderr,*)

    moi=inv(moi)

    do ni=1,num_ions
       ion_absolute(:,ni)=matmul(moi,ion_absolute(:,ni)-com)
    end do

    write (ctemp,'(i10)') num_ions
    write (stdout,'(a/)') trim(adjustl(ctemp))
    do ni=1,num_ions
       write (stdout,'(a,3f20.13)') trim(adjustl(species_names(ion_index(ni)))),ion_absolute(:,ni)
    end do

  end subroutine write_xyz

  ! ** Extended XYZ file format

  subroutine read_xyze()

    integer :: ni,indx1,indx2
    real(kind=dp) :: lat(9),com(3)

    read(buff_words(1,1),*) num_ions

    indx1=index(buff(2),'"')
    indx2=index(buff(2),'"',back=.true.)

    read(buff(2)(indx1+1:indx2-1),*) lat

    lattice_car=reshape(lat,(/3,3/))

    call init_ions()

    com=0.0_dp
    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,ni+2))
       read(buff_words(2:4,ni+2),*) ion_absolute(:,ni)
       com=com+ion_absolute(:,ni)/real(num_ions,dp)
    end do

    do ni=1,num_ions
       ion_absolute(:,ni)=ion_absolute(:,ni)-com(:)
    end do

  end subroutine read_xyze

  subroutine write_xyze()

    integer :: ni
    character(len=10) :: ctemp

    write (ctemp,'(i10)') num_ions
    write (stdout,'(a)') trim(adjustl(ctemp))
    if(any(ion_spin.ne.0.0_dp)) then
       write (stdout,'(a,9f20.6,a)') 'Lattice="',reshape(lattice_car,(/9/)),'" Properties=species:S:1:pos:R:3:spin:R:1'
       do ni=1,num_ions
          write (stdout,'(a,4f20.13)') trim(adjustl(species_names(ion_index(ni)))),ion_absolute(:,ni),ion_spin(ni)
       end do
    else
       write (stdout,'(a,9f20.6,a)') 'Lattice="',reshape(lattice_car,(/9/)),'" Properties=species:S:1:pos:R:3'
       do ni=1,num_ions
          write (stdout,'(a,3f20.13)') trim(adjustl(species_names(ion_index(ni)))),ion_absolute(:,ni)
       end do
    end if
  end subroutine write_xyze

  ! ** gulp file format

  subroutine read_gulp()
    stop 'read_gulp: not implemented'
  end subroutine read_gulp

  subroutine write_gulp()

    integer :: ni    

    write (stdout,'(a/a/a/a/a)') 'opti prop','title','Output','end','cell'
    write (stdout,'(6f8.3)') lattice_abc
    write (stdout,'(a)') 'frac'
    do ni=1,num_ions
       write (stdout,'(a,a6,3f18.13)') trim(adjustl(species_names(ion_index(ni)))),'core',ion_fractional(:,ni)
    end do

  end subroutine write_gulp

  ! ** xtl file format

  subroutine read_xtl()

    integer :: n,nstart,nfinish,ni

    do n=1,num_lines
       if(up(buff_words(1,n)).eq.'CELL') then
          read(buff_words(1:6,n+1),*) lattice_abc
          exit
       end if
    end do

    nstart=0
    nfinish=0
    do n=1,num_lines
       if(up(buff_words(1,n)).eq.'EOF') then
          nfinish=n
          exit
       end if
       if((up(buff_words(1,n)).eq.'ATOMS').and.((up(buff_words(1,n+1)).eq.'NAME'))) nstart=n+1
    end do

    if(nstart.lt.nfinish) then
       if(nfinish-nstart.lt.2) stop 'cabal: problem reading positions in xtl data'
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          read(buff_words(1:4,nstart+ni),*) ion_names(ni),ion_fractional(1:3,ni)
          ion_names(ni)=trim(strip(ion_names(ni)))
       end do

    end if

  end subroutine read_xtl

  subroutine write_xtl()

    integer :: n

    write (stdout,'(4a)') 'TITLE cabal-',trim(inext),'2',trim(outext)
    write (stdout,'(a/,6f11.5)') 'CELL',lattice_abc
    write (stdout,'(a/a)') 'ATOMS','NAME         X           Y           Z'

    do n=1,num_ions
       write (stdout,'(a4,3f17.13)') &
            species_names(ion_index(n)),ion_fractional(:,n)
    end do
    write (stdout,'(a)') 'EOF'

  end subroutine write_xtl

  ! ** cif file format

  subroutine read_cif()

    stop 'read_cif: not implemented'

  end subroutine read_cif

  subroutine write_cif()

    integer :: ni

    write (stdout,'(a/)')       "data_cif"
    write (stdout,'(a/)')       "_audit_creation_method             'generated by cabal'"
    write (stdout,'(a)')        "_symmetry_space_group_name_H-M     'P1'"
    write (stdout,'(a)')        "_symmetry_Int_Tables_number        1"
    write (stdout,'(a/)')       "_symmetry_cell_setting             triclinic"
    write (stdout,'(a,f10.4)')  "_cell_length_a                     ",lattice_abc(1) 
    write (stdout,'(a,f10.4)')  "_cell_length_b                     ",lattice_abc(2) 
    write (stdout,'(a,f10.4)')  "_cell_length_c                     ",lattice_abc(3) 
    write (stdout,'(a,f10.4)')  "_cell_angle_alpha                  ",lattice_abc(4) 
    write (stdout,'(a,f10.4)')  "_cell_angle_beta                   ",lattice_abc(5) 
    write (stdout,'(a,f10.4)')  "_cell_angle_gamma                  ",lattice_abc(6) 
    write (stdout,'(a)')        "loop_"
    write (stdout,'(a)')        "_atom_site_label"
    write (stdout,'(a)')        "_atom_site_fract_x"
    write (stdout,'(a)')        "_atom_site_fract_y"
    write (stdout,'(a)')        "_atom_site_fract_z"
    write (stdout,'(a)')        "_atom_site_occupancy"
    do ni=1,num_ions
       write (stdout,'(a5,4f9.5)') trim(species_names(ion_index(ni))),ion_fractional(:,ni),ion_occ(ni)
    end do

  end subroutine write_cif

  ! ** psi4 file format

  subroutine read_psi4()

    real(kind=dp) :: com(3)

    integer :: ni,nstart,nfinish,n

    nstart=0
    nfinish=0
    do n=1,num_lines
       if(index(buff(n),'Cartesian Geometry (in Angstrom)').gt.0) nstart=n
       if(index(buff(n),'Saving final (previous) structure.').gt.0) nfinish=n
    end do

    if(nfinish-nstart.lt.2) then
       stop 'cabal: problem reading atomic positions from psi4 data'
    else
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          ion_names(ni)=trim(buff_words(1,nstart+ni))
          read(buff_words(2:4,nstart+ni),*) ion_absolute(1:3,ni)
       end do

    end if

    com=sum(ion_absolute,2)/real(num_ions,dp)

    lattice_car=ident*maxval(abs(ion_absolute))*8.0_dp

    do ni=1,num_ions
       ion_absolute(:,ni)=ion_absolute(:,ni)-com+sum(lattice_car,2)/2.0_dp
    end do

  end subroutine read_psi4

  subroutine write_psi4()

    integer :: ni

    write (stdout,'(a)') 'molecule {'
    do ni=1,num_ions
       write (stdout,'(a,3f18.13)') trim(species_names(ion_index(ni))),ion_absolute(:,ni)
    end do
    write (stdout,'(a)') '}'

  end subroutine write_psi4

  ! ** DL_POLY config file format

  subroutine read_config()

    integer :: i,n,indx

    character(len=10) :: ctemp

    read(buff_words(3,2),*) num_ions

    do i=1,3
       read(buff(2+i),*) lattice_car(:,i)
    end do

    call init_ions()

    do n=1,num_ions
       ctemp=trim(buff_words(1,4*n+2))
       indx=index(ctemp,'_')
       if(indx.gt.1) ctemp=ctemp(:indx-1)
       ctemp(1:1)=up(ctemp(1:1))
       ion_names(n)=ctemp
       read(buff(4*n+3),*) ion_absolute(:,n)
    end do

  end subroutine read_config

  subroutine write_config()

    stop 'write_config not implemented'

  end subroutine write_config

  ! ** VASP POSCAR format

  subroutine read_poscar()

    integer :: i,j,n,nisp

    real(kind=dp) :: scale    

    read(buff(2),*) scale
    do i=1,3
       read(buff(2+i),*) lattice_car(:,i)
    end do
    lattice_car=lattice_car*scale

    num_species=num_words(6)

    num_ions=0
    do i=1,num_species
       read(buff_words(i,7),*) n
       num_ions=num_ions+n
    end do

    call init_ions()

    n=0
    do i=1,num_species
       read(buff_words(i,7),*) nisp
       do j=1,nisp
          n=n+1
          ion_names(n)=trim(buff_words(i,6))
       end do
    end do

    if(trim(adjustl(buff(8))).ne."Direct") stop 'cabal: read_poscar - direct only'

    do i=1,num_ions
       read(buff(8+i),*) ion_fractional(:,i)
    end do

  end subroutine read_poscar

  subroutine write_poscar()

    integer :: n,ns,num_type(num_species)

    write (stdout,'(a)') 'cabal: POSCAR file converted at '//trim(date())
    write (stdout,*) 1.0_dp
    write (stdout,'(3f15.10)') lattice_car
    write (stdout,*) species_names(1:num_species)
    do n=1,num_species
       num_type(n)=count(ion_names.eq.species_names(n))
    end do
    write (stdout,*) num_type
    write (stdout,'(a)') ' Direct'
    do ns=1,num_species
       do n=1,num_ions
          if(ion_names(n).eq.species_names(ns)) write (stdout,'(3f15.10)') ion_fractional(:,n)
       end do
    end do

  end subroutine write_poscar

  ! ** Quantum Espresso input file format
  ! ** This implementation requires CELL_PARAMETERS and ATOMIC_POSITIONS to be present
  !    (i.e. ibrav =0)

  subroutine read_qe()

    integer :: i,n

    ! Loop through and find CELL_PARAMETERS
    do n=1,num_lines

       if(buff_words(1,n).eq.'CELL_PARAMETERS') then
          ! Decide whether the units are angstrom or bohr
          if (buff_words(2,n).eq.'angstrom') then
             do i=1,3
                read(buff_words(1:3,n+i),*) lattice_car(1:3,i)
             end do
          else if(buff_words(2,n).eq.'bohr') then
             do i=1,3
                read(buff_words(1:3,n+i),*) lattice_car(1:3,i)
             end do
             lattice_car = lattice_car*bohr2ang
          else
             write (stderr,'(a)') 'cabal: read_qe'
             write (stderr,'(a)') 'CELL_PARAMETERS units not recognised - expected angstrom or bohr'
             stop 
          end if
          exit
       end if

       if(n.eq.num_lines) then
          write (stderr,'(a)') 'cabal: read_qe'
          write (stderr,'(a)') 'Reached end of input file and CELL_PARAMETERS not found.'
          stop 
       end if

    end do

    ! Do ATOMIC_POSITIONS    
    do n=1,num_lines
       if(buff_words(1,n).eq.'nat') then
          read(buff_words(3,n),*) num_ions
          exit
       end if
       if(n.eq.num_lines) then
          write (stderr,'(a)') 'cabal: read_qe'
          write (stderr,'(a)') 'Reached end of input file and number of atoms (nat = ...) not found'
          write (stderr,'(a)') '[ nat should appear on its own line, in the form nat = xxx , ]'
          stop 
       end if
    end do

    call init_ions()

    do n=1,num_lines
       if(buff_words(1,n).eq.'ATOMIC_POSITIONS') then
          ! Decide whether the positions are bohr, angstrom, or crystal
          if(buff_words(2,n).eq.'bohr') then
             do i=1,num_ions
                ion_names(i)=trim(buff_words(1,n+i))
                read(buff_words(2:4,n+i),*) ion_absolute(1:3,i)
                ion_absolute = ion_absolute*bohr2ang
             end do
          else if(buff_words(2,n).eq.'angstrom') then
             do i=1,num_ions
                ion_names(i)=trim(buff_words(1,n+i))
                read(buff_words(2:4,n+i),*) ion_absolute(1:3,i)
             end do
          else if(buff_words(2,n).eq.'crystal') then
             do i=1,num_ions
                ion_names(i)=trim(buff_words(1,n+i))
                read(buff_words(2:4,n+i),*) ion_fractional(1:3,i)
             end do
          else
             write (stderr,'(a)') 'cabal: read_qe'
             write (stderr,'(a)') 'ATOMIC_POSITIONS units not recognised - must be bohr, angstrom or crystal'
             write (stderr,'(a)') '[ alat and crystal_sg are not supported yet ]'
             stop
          end if
          exit
       end if
       if(n.eq.num_lines) then
          write (stderr,'(a)') 'cabal: read_qe'
          write (stderr,'(a)') 'Reached end of input file and ATOMIC_POSITIONS not found.'
          stop 
       end if
    end do

  end subroutine read_qe

  subroutine write_qe()

    integer :: n

    write (stdout,'(a)') '&SYSTEM'
    write (stdout,'(a)') '  ibrav = 0 ,'
    write (stdout,'(a,i5,a)') '  nat   = ',num_ions,' ,'
    write (stdout,'(a,i5,a)') '  ntyp  = ',num_species,' ,'
    write (stdout,*)
    write (stdout,'(a)') 'CELL_PARAMETERS angstrom'
    write (stdout,'(3f20.14)') lattice_car
    write (stdout,*)
    write (stdout,'(a)') 'ATOMIC_POSITIONS crystal'
    do n=1,num_ions
       write (stdout,'(a,3f20.14)') trim(adjustl(species_names(ion_index(n)))),ion_fractional(:,n)
    end do
    write (stdout,*)
    if (kpts) then
       write (stdout,'(a)') 'K_POINTS automatic'
       write (stdout,'(6i4)') mp_kgrid(1),mp_kgrid(2),mp_kgrid(3),0,0,0
       write (stdout,*)
    end if

  end subroutine write_qe

  ! ** CP2K file format

  subroutine read_cp2k()

    integer   ::  i,n,nstart,nfinish
    logical   ::  cart,abc,abc_noang,scaled

    cart=.false.
    abc=.false.
    abc_noang=.false.
    ! We look for *consecutive* lines beginning with A,B,C, or with ABC,ALPHA_BETA_GAMMA
    do n=1,num_lines
       if((buff_words(1,n).eq.'A').and.(buff_words(1,n+1).eq.'B').and.(buff_words(1,n+2).eq.'C')) then
          cart=.true.
       end if
       if((buff_words(1,n).eq.'ABC').and.(buff_words(1,n+1).eq.'ALPHA_BETA_GAMMA')) then
          abc=.true.
       end if
       if((buff_words(1,n).eq.'ABC').and.(buff_words(1,n+1).ne.'ALPHA_BETA_GAMMA')) then
          abc_noang=.true.
       end if
    end do
    ! Check that only one of cart, abc, and abc_noang, is true
    if((cart.and.abc).or.(cart.and.abc_noang).or.(abc.and.abc_noang).or.(cart.and.abc.and.abc_noang)) then
       write (stderr,'(a)') 'cabal: read_cp2k'
       write (stderr,'(a)') 'Lattice multiply defined: only one of A,B,C or'
       write (stderr,'(a)') 'ABC/ALPHA_BETA_GAMMA should be present'
       stop
    end if
    if((.not.cart).and.(.not.abc).and.(.not.abc_noang)) then
       write (stderr,'(a)') 'cabal: read_cp2k'
       write (stderr,'(a)') 'Lattice information (A,B,C or ABC/ALPHA_BETA_GAMMA) not found'
       stop
    end if
    ! Caution about angles (not a stopping condition)
    if(abc_noang) then
       write (stderr,'(a)') 'cabal: read_cp2k'
       write (stderr,'(a)') 'Remark: lattice angles not explicitly given'
       write (stderr,'(a)') 'Defaulting to alpha = beta = gamma = 90 degrees'
    end if
    ! Now read and store lattice information
    do n=1,num_lines
       if((buff_words(1,n).eq.'A').and.(buff_words(1,n+1).eq.'B').and.(buff_words(1,n+2).eq.'C')) then
          do i=0,2
             read(buff_words(2:4,n+i),*) lattice_car(1:3,i+1)
          end do
       else if((buff_words(1,n).eq.'ABC').and.(buff_words(1,n+1).eq.'ALPHA_BETA_GAMMA')) then
          read(buff_words(2:4,n),*)   lattice_abc(1:3)
          read(buff_words(2:4,n+1),*) lattice_abc(4:6)
       else if((buff_words(1,n).eq.'ABC').and.(buff_words(1,n+1).ne.'ALPHA_BETA_GAMMA')) then
          read(buff_words(2:4,n),*)   lattice_abc(1:3)
          lattice_abc(4) = 90.0_dp
          lattice_abc(5) = 90.0_dp
          lattice_abc(6) = 90.0_dp
       end if
    end do

    ! Do atomic positions
    num_ions = 0 
    nstart = 0
    nfinish = 0
    scaled=.false.
    do n=1,num_lines
       if(buff_words(1,n).eq.'&COORD') then
          nstart = n+1
          do i=nstart,num_lines
             if((buff_words(1,i).ne.'UNIT').and.(buff_words(1,i).ne.'SCALED').and.(buff_words(1,i).ne.'&END')) then
                num_ions = num_ions + 1
             end if
             if((buff_words(1,i).eq.'SCALED').and.(buff_words(2,i).eq.'T')) then
                scaled=.true.
             end if
             if(buff_words(1,i).eq.'&END') then
                nfinish = i-1
                exit
             end if
             if(i.eq.num_lines) then
                write (stderr,'(a)') 'cabal: read_cp2k'
                write (stderr,'(a)') '&COORD found, but end of file reached before closing &END'
                stop
             end if
          end do
          exit   
       end if
       if(n.eq.num_lines) then
          write (stderr,'(a)') 'cabal: read_cp2k'
          write (stderr,'(a)') 'No atom coordinate information found (&COORD not present)'
          stop
       end if
    end do
    if(nstart.gt.nfinish) then
       write (stderr,'(a)') 'cabal: read_cp2k'
       write (stderr,'(a)') 'Error assigning nstart and nfinish: nstart > nfinish'
       stop
    end if

    call init_ions()

    i=0
    do n=nstart,nfinish
       if((buff_words(1,n).ne.'UNIT').and.(buff_words(1,n).ne.'SCALED')) then
          i = i+1
          if(scaled) then
             ! Coordinates are fractional
             ion_names(i)=trim(buff_words(1,n))
             read(buff_words(2:4,n),*) ion_fractional(1:3,i)
          else if(.not.scaled) then
             ! Coordinates are absolute
             ion_names(i)=trim(buff_words(1,n))
             read(buff_words(2:4,n),*) ion_absolute(1:3,i)
          end if
       end if
    end do
    if(i.ne.num_ions) then
       write (stderr,'(a)') 'cabal: read_cp2k'
       write (stderr,'(a)') 'Inconsistency assigning ion arrays on read-in'
       stop
    end if

    ! Finally, we look for UNITS in the case of non-fractional coordinates, and convert if need be
    if(.not.scaled) then
       do n=nstart,nfinish
          if(buff_words(1,n).eq.'UNIT') then
             if(buff_words(2,n).eq.'bohr') then
                ion_absolute = ion_absolute*bohr2ang
             else if(buff_words(2,n).eq.'m') then
                ion_absolute = ion_absolute*1e10_dp
             else if(buff_words(2,n).eq.'pm') then
                ion_absolute = ion_absolute*0.01_dp
             else if(buff_words(2,n).eq.'nm') then
                ion_absolute = ion_absolute*10_dp
             else if(buff_words(2,n).eq.'angstrom') then
                ! Nothing needs doing
             else
                write (stderr,'(a)') 'cabal: read_cp2k'
                write (stderr,'(a)') 'UNITS for atomic coordinates (&COORD) not recognised'
                write (stderr,'(a)') 'Must be bohr, m, pm, nm, or angstrom'
                stop            
             end if
             exit
          end if
       end do
    end if

  end subroutine read_cp2k

  subroutine write_cp2k()

    integer :: n

    write (stdout,'(a)') '&SYSTEM'
    write (stdout,'(a)') '  ibrav = 0 ,'
    write (stdout,'(a,i5,a)') '  nat   = ',num_ions,' ,'
    write (stdout,'(a,i5,a)') '  ntyp  = ',num_species,' ,'
    write (stdout,*)
    write (stdout,'(a)') 'CELL_PARAMETERS angstrom'
    write (stdout,'(3f20.14)') lattice_car
    write (stdout,*)
    write (stdout,'(a)') 'ATOMIC_POSITIONS crystal'
    do n=1,num_ions
       write (stdout,'(a,3f20.14)') trim(adjustl(species_names(ion_index(n)))),ion_fractional(:,n)
    end do
    write (stdout,*)

  end subroutine write_cp2k

  subroutine read_conf()

    integer :: n,is

    read(buff_words(1,3),*) num_ions
    read(buff_words(1,4),*) num_species

    call init_ions()

    do is=1,num_species
       species_names(is) = trim(buff_words(is,1))
    end do

    lattice_car=0.0_dp
    do n=1,3
       read(buff_words(2,5+n),*) lattice_car(n,n)
    end do
    read(buff_words(1:2,9),*) lattice_car(1,2:3)
    read(buff_words(3,9),*) lattice_car(2,3)

    do n=1,num_ions
       read(buff_words(2,15+num_species+n),*) is
       ion_names(n)=species_names(is)
       read(buff_words(4:6,15+num_species+n),*) ion_absolute(1:3,n)
    end do

  end subroutine read_conf

  subroutine write_conf()

    integer :: n,ns

    real(kind=dp) :: lattice_lammps(3,3),ax,bx,by,cx,cy,cz

    character(len=30) :: ctemp

    do ns=1,num_species
       write (stdout,'(a3)',advance='no') species_names(ns)
    end do
    write (stdout,*)
    write (stdout,*)
    write (ctemp,*) num_ions,'atoms'
    write (stdout,'(a)') trim(adjustl(ctemp))
    write (ctemp,*) num_species,'atom types'
    write (stdout,'(a)') trim(adjustl(ctemp))
    write (stdout,*)


    ax=norm2(lattice_car(:,1))
    bx=dot_product(lattice_car(:,2),lattice_car(:,1))/ax
    by=sqrt(norm2(lattice_car(:,2))**2-bx**2)
    cx=dot_product(lattice_car(:,3),lattice_car(:,1))/ax
    cy=(dot_product(lattice_car(:,2),lattice_car(:,3))-bx*cx)/by
    cz=sqrt(norm2(lattice_car(:,3))**2-cx**2-cy**2)

    lattice_lammps=0.0_dp
    lattice_lammps(1,1)=ax
    lattice_lammps(1,2)=bx    
    lattice_lammps(2,2)=by
    lattice_lammps(1,3)=cx
    lattice_lammps(2,3)=cy
    lattice_lammps(3,3)=cz   

    write (stdout,'(f7.4,f8.4,a)') 0.0_dp,lattice_lammps(1,1),' xlo xhi'
    write (stdout,'(f7.4,f8.4,a)') 0.0_dp,lattice_lammps(2,2),' ylo yhi'
    write (stdout,'(f7.4,f8.4,a)') 0.0_dp,lattice_lammps(3,3),' zlo zhi'
    write (stdout,'(f7.4,2f8.4,a)') lattice_lammps(1,2:3),lattice_lammps(2,3),' xy xz yz'
    write (stdout,*)

    write (stdout,'(a)') 'Masses'
    write (stdout,*)
    do ns=1,num_species
       write (ctemp,'(i2,a4)') ns," 1.0"
       write (stdout,'(a)') trim(adjustl(ctemp))
    end do
    write (stdout,*)
    write (stdout,'(a)') 'Atoms'
    write (stdout,*)

    do ns=1,num_species
       do n=1,num_ions
          if(ion_names(n).eq.species_names(ns))  write (stdout,'(i5,i6,4f15.8)') n,ns,0.0_dp,&
               matmul(lattice_lammps,ion_fractional(:,n))
       end do
    end do

  end subroutine write_conf

  subroutine write_json()

    integer :: ni,nj
    character(len=15) :: str

    write (stdout,'(a,i5,a)')  '            "size": ',num_ions,','
    write (stdout,'(a)',advance='no') '            "names": ['
    do ni=1,num_ions-1
       write (stdout,'(3a)',advance='no') '"',trim(adjustl(species_names(ion_index(ni)))),'", '
    end do
    write (stdout,'(3a)') '"',trim(adjustl(species_names(ion_index(num_ions)))),'"],'
    write (stdout,'(a)',advance='no') '            "x": ['
    do ni=1,num_ions-1
       write (str,'(f15.8)') ion_absolute(1,ni)
       write (stdout,'(a,a)',advance='no') trim(adjustl(str)),', '
    end do
    write (str,'(f15.8)') ion_absolute(1,ni)
    write (stdout,'(a,a)') trim(adjustl(str)),'],'
    write (stdout,'(a)',advance='no') '            "y": ['
    do ni=1,num_ions-1
       write (str,'(f15.8)') ion_absolute(2,ni)
       write (stdout,'(a,a)',advance='no') trim(adjustl(str)),', '
    end do
    write (str,'(f15.8)') ion_absolute(2,ni)
    write (stdout,'(a,a)') trim(adjustl(str)),'],'
    write (stdout,'(a)',advance='no') '            "z": ['
    do ni=1,num_ions-1
       write (str,'(f15.8)') ion_absolute(3,ni)
       write (stdout,'(a,a)',advance='no') trim(adjustl(str)),', '
    end do
    write (str,'(f15.8)') ion_absolute(3,ni)
    write (stdout,'(a,a)') trim(adjustl(str)),'],'
    write (stdout,'(a)',advance='no') '            "cell": ['
    do ni=1,3
       do nj=1,3
          if(ni.eq.3.and.nj.eq.3) then
             write (str,'(f15.8)') lattice_car(nj,ni)
             write (stdout,'(a,a)') trim(adjustl(str)),']'
          else
             write (str,'(f15.8)') lattice_car(nj,ni)
             write (stdout,'(a,a)',advance='no') trim(adjustl(str)),', '
          end if
       end do
    end do

  end subroutine write_json

  !* =======================================================================

  subroutine init_ions()

    if(allocated(ion_names)) stop 'cabal: multiple ion definitions'

    allocate(ion_names(num_ions*4),ion_index(num_ions*4),ion_occ(num_ions*4),ion_spin(num_ions*4))
    allocate(ion_fractional(3,num_ions*4),ion_absolute(3,num_ions*4))

    ion_names=''
    ion_index=0
    ion_occ=1.0_dp
    ion_spin=0.0_dp
    ion_fractional=0.0_dp
    ion_absolute=0.0_dp

  end subroutine init_ions

  subroutine car2abc(car,abc)

    real(kind=dp), intent(in)  :: car(3,3)
    real(kind=dp), intent(out) :: abc(6)

    abc(1) = sqrt(car(1,1)**2+car(2,1)**2+car(3,1)**2)
    abc(2) = sqrt(car(1,2)**2+car(2,2)**2+car(3,2)**2)
    abc(3) = sqrt(car(1,3)**2+car(2,3)**2+car(3,3)**2)
    abc(4) = acos(dot_product(car(:,2),car(:,3))/abc(2)/abc(3))/dgrd
    abc(5) = acos(dot_product(car(:,1),car(:,3))/abc(1)/abc(3))/dgrd
    abc(6) = acos(dot_product(car(:,1),car(:,2))/abc(1)/abc(2))/dgrd

  end subroutine car2abc

  subroutine abc2car(abc,car)

    real(kind=dp), intent(in) :: abc(6)
    real(kind=dp), intent(out):: car(3,3)

    car(:,1) = (/abc(1),0.0_dp,0.0_dp/)
    car(:,2) = (/abc(2)*cos(dgrd*abc(6)),abc(2)*sin(dgrd*abc(6)),0.0_dp/)
    car(1,3) = abc(3)*cos(dgrd*abc(5))
    car(2,3) = abc(3)*(cos(dgrd*abc(4))-cos(dgrd*abc(5))*cos(dgrd*abc(6)))/sin(dgrd*abc(6))
    car(3,3) = sqrt(abc(3)**2-car(1,3)**2-car(2,3)**2)

  end subroutine abc2car

  subroutine car2rec(car,rec)

    real(kind=dp), intent(in)  :: car(3,3)
    real(kind=dp), intent(out) :: rec(3,3)

    real(kind=dp) :: volume

    volume = car(1,1)*(car(2,2)*car(3,3)-car(3,2)*car(2,3))+&
         car(2,1)*(car(3,2)*car(1,3)-car(1,2)*car(3,3))+&
         car(3,1)*(car(1,2)*car(2,3)-car(2,2)*car(1,3))

    if(abs(volume).lt.epsilon(1.0_dp)) stop 'cabal: zero volume cell detected'

    ! ** Calculate the reciprocal lattice vectors

    rec(1,1)=car(2,2)*car(3,3)-car(3,2)*car(2,3)
    rec(2,1)=car(2,3)*car(3,1)-car(3,3)*car(2,1)
    rec(3,1)=car(2,1)*car(3,2)-car(3,1)*car(2,2)
    rec(1,2)=car(3,2)*car(1,3)-car(1,2)*car(3,3)
    rec(2,2)=car(3,3)*car(1,1)-car(1,3)*car(3,1)
    rec(3,2)=car(3,1)*car(1,2)-car(1,1)*car(3,2)
    rec(1,3)=car(1,2)*car(2,3)-car(2,2)*car(1,3)
    rec(2,3)=car(1,3)*car(2,1)-car(2,3)*car(1,1)
    rec(3,3)=car(1,1)*car(2,2)-car(2,1)*car(1,2)

    rec(:,:)=rec(:,:)/volume

    ! ** Remark: 'rec' is the actual inverse of 'car'
    !    With 'car' having column lattice vectors, an additional transpose
    !    would be needed to get 'rec' having column reciprocal lattice vectors

  end subroutine car2rec

  function car_volume(car)

    real(kind=dp), dimension(3,3), intent(in) :: car

    real(kind=dp) :: car_volume

    car_volume = car(1,1)*(car(2,2)*car(3,3)-car(3,2)*car(2,3))+&
         car(2,1)*(car(3,2)*car(1,3)-car(1,2)*car(3,3))+&
         car(3,1)*(car(1,2)*car(2,3)-car(2,2)*car(1,3))

    if(abs(car_volume).lt.epsilon(1.0_dp)) stop 'cabal: zero volume cell detected'

  end function car_volume

  function up(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=240)            :: up

    integer :: i,ic,ia,iz,ishift

    ! ** Initialise system dependant bounds in collating sequence

    ia     = ichar('a')
    iz     = ichar('z') 
    ishift = ichar('A')-ia

    ! ** Raise all lowercase to upper case

    up=string

    do i=1,200
       ic = ichar(up(i:i))
       if((ic.ge.ia).and.(ic.le.iz)) up(i:i) = char(ishift+ic) 
    end do

  end function up

  function strip(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: strip

    integer                      :: i,ia,iz,iAA,iZZ,itest

    ! ** Initialise system dependant bounds in collating sequence

    ia=ichar('a')
    iz=ichar('z')
    iAA=ichar('A')
    iZZ=ichar('Z')

    ! ** Strip all characters after the first non-alphabetical character

    do i=1,len(string)
       itest=ichar(string(i:i))
       if(((itest.lt.ia).or.(itest.gt.iz)).and.((itest.lt.iAA).or.(itest.gt.iZZ))) exit
    end do

    strip=string(1:i-1)

  end function strip

  function detab(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: detab

    integer                      :: i

    detab=string
    do while (scan(detab,char(9)).gt.0)
       i=scan(detab,char(9))
       detab=trim(detab(1:i-1))//'      '//trim(detab(i+1:))
    end do

  end function detab

  function deunder(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: deunder

    integer                      :: i

    deunder=string
    do while (scan(deunder,'_').gt.0)
       i=scan(deunder,'_')
       deunder=trim(deunder(1:i-1))//trim(deunder(i+1:))
    end do

  end function deunder
  
  function date()

    implicit none

    character(len=80) :: date
    character(len=80) :: time,month,zone
    character(len=2)  :: suffix  

    call date_and_time(date,time,zone)

    ! ** Reformat into human readable form

    ! * First, the time (including zone information)

    if(zone(2:2).eq.'0') then
       zone = '(GMT'//zone(1:1)//zone(3:3)//'.'//zone(4:4)//')'
    else
       zone = '(GMT'//zone(1:1)//zone(2:3)//'.'//zone(4:4)//')'
    end if

    time = time(1:2)//':'//time(3:4)//':'//time(5:6)//' '//trim(zone)

    ! * Now, the date

    ! Find the month

    select case(date(5:6))
    case ('01') ; month='January'
    case ('02') ; month='February'
    case ('03') ; month='March'
    case ('04') ; month='April'
    case ('05') ; month='May'
    case ('06') ; month='June'
    case ('07') ; month='July'
    case ('08') ; month='August'
    case ('09') ; month='September'
    case ('10') ; month='October'
    case ('11') ; month='November'
    case ('12') ; month='December'
    case default
       month = 'Month Unknown'
    end select

    ! Date suffix

    select case(date(7:8))
    case ('01') ; suffix='st'
    case ('02') ; suffix='nd'
    case ('03') ; suffix='rd'
    case ('21') ; suffix='st'
    case ('22') ; suffix='nd'
    case ('23') ; suffix='rd'
    case ('31') ; suffix='st'
    case default
       suffix='th'
    end select

    ! * Finally, combine into a single string

    if(date(7:7).eq.'0') then
       date = trim(time)//' '//date(8:8)//suffix//' '//&
            trim(month)//' '//date(1:4)
    else
       date = trim(time)//' '//date(7:8)//suffix//' '//&
            trim(month)//' '//date(1:4)
    end if

  end function date

  function outer_product(a,b) result(c)

    real(kind=dp), intent(in) :: a(3),b(3)

    real(kind=dp) :: c(3,3)

    c=spread(a,dim=2,ncopies=3)*spread(b,dim=1,ncopies=3)

  end function outer_product

  subroutine diag(a,b) 

    real(kind=dp), intent(inout) :: a(3,3)

    real(kind=dp), intent(out) :: b(3)

    integer :: info

    real(kind=dp) :: work(12)

    call dsyev('V','U',3,a,3,b,work,12,info)

    if(info.ne.0) stop 'error in diag'

  end subroutine diag

  function inv(A) result(Ainv)

    real(kind=dp), dimension(:,:), intent(in) :: A
    real(kind=dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(kind=dp), dimension(size(A,1)) :: work 
    integer,       dimension(size(A,1)) :: ipiv
    integer :: n, info

    external dgetrf
    external dgetri

    Ainv=A
    n=size(A,1)

    call dgetrf(n,n,Ainv,n,ipiv,info)

    if (info.ne.0)  stop 'inv: matrix is numerically singular!'

    call dgetri(n,Ainv,n,ipiv,work,n,info)

    if (info.ne.0) stop 'inv: matrix inversion failed!'

  end function inv

end program cabal
