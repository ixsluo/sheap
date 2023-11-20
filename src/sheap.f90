!=================================================================================!
!                                      SHEAP                                      !
!=================================================================================!
!                                                                                 !
! This is the main program of SHEAP --- a dimensionality reduction tool           !
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

program sheap

  ! * Prerequisite modules
  use constants
  use control
  use model
  use optim
  use rng

  ! * All variables must be declared
  implicit none

  ! * Write out SHEAP banner
  if(.not.quiet) call banner()

  ! * Initialise random number generator 
  call init_pseudorandom()

  ! * Read in command line arguments
  call control_arguments()

  ! * Read in source-data vector file
  if(verbose) write(stderr,*) 'reading in vector file'
  if(file_type.eq."vec") then
     call model_read_vec()
  else if(file_type.eq."old") then
     call model_read_old()
  else if(file_type.eq."xyz") then
     call model_read_xyz()
  else
     write(stderr,*) 'invalid file type selected'
     stop
  end if

  ! * Write out control module report
  if(.not.quiet) call control_report()

  ! * Compute neighbourhood graph of source-data
  call model_init()

  ! * Initialise embedding - place N data points in low-dimensional space
  if(verbose) write(stderr,'(a)',advance='NO') ' initialising embedding with:'
  if(uniform_packing) then
     if(verbose) write(stderr,*) 'uniform packing'
     call model_uniform_packing()
  else if(random_projection) then
     if(verbose) write(stderr,*) 'random projection'
     call model_random_projection()
  else if(pca) then
     if(verbose) write(stderr,*) 'pca'
     call model_pca()
  else if(resume) then
     if(verbose) write(stderr,*) 'resuming from previous projection'
     call model_resume()
  end if
  if(verbose.and.init_only) write(stderr,*) 'stopping - initialise only set'

  ! * Write out model module report
  if(.not.quiet) call model_report()

  ! * Optimisation
  if(init_only) then

     ! * Write out result of initialisation
     call model_write(final_cost)

  else

     ! * Test analytical gradient against numerical scheme - does not carry out optimisation
     if(grad_test) then
        write(stderr,*) 'computing numerical gradient'
        call optim_test();stop
     end if

     ! * Write out choice of cost function to be optimised
     if(verbose.and.use_kl) write(stderr,*) 'using KL-divergence cost function'
     if(verbose.and.use_ce) write(stderr,*) 'using cross-entropy cost function'

     ! * Write out choice of optimisation scheme
     if(.not.quiet.and.use_tpsd) write(stderr,*) 'doing tpsd'
     if(.not.quiet.and.use_sgd) write(stderr,*) 'doing sgd'

     ! * Minimise cost function to optimie layout of map
     if(.not.init_only) call optim_ise(tolerance,niter)        

     ! * Write map points to output
     if(.not.trace) call model_write(final_cost)

     ! * Write out final values from optimisation
     write (stderr,*)
     write (stderr,'(a25,a3,f18.10)') 'Final cost',' : ',final_cost
     write (stderr,'(a25,a3,i13)') 'Number of TPSD steps',' : ',steps
     if(steps.ne.total_steps) then
        write (stderr,'(a25,a3,i13)') 'Total number of optimisation steps',' : ',total_steps   
     end if
     if(grow) write (stderr,'(a25,a3,i13)') 'HS growth from step',' : ',start_growth

  end if

contains

  subroutine banner()

    write (stderr,*) '                                              '
    write (stderr,*) '                            8888              ' 
    write (stderr,*) '               8888888  8888    88            '
    write (stderr,*) '              88      88          888         '
    write (stderr,*) '             8                      8         '
    write (stderr,*) '             88                     88        '
    write (stderr,*) '              88 888         88   888         '
    write (stderr,*) '         .8888888---888888888--88888888.      '
    write (stderr,*) '       8888  888-----------------88   888     '
    write (stderr,*) '      88     888--8888-----8888--888    88    '
    write (stderr,*) '     88     888---88 8-----88 8--888      8   '
    write (stderr,*) '    88      8 8---8888-----8888---88       8  '
    write (stderr,*) '    8       8 8-------------------88       8  '
    write (stderr,*) '    88      8 8-------------------8 88    88  '
    write (stderr,*) '     88    88 8------------------8  88   88   '
    write (stderr,*) '       8888   88----88-----88----8    8888    '
    write (stderr,*) '               88-----8---8----88             '
    write (stderr,*) '                88------------88              '
    write (stderr,*) '                  88--------88                '
    write (stderr,*) '                    "888888"                  '
    write (stderr,*) '                                              '
    write (stderr,*) '           888                                '
    write (stderr,*) '           888                                '
    write (stderr,*) '           888                                '
    write (stderr,*) '  .d8888b  88888b.   .d88b.   8888b.  88888b. '
    write (stderr,*) '  888      888 "88b d88  88b     "88b 888 "88b'
    write (stderr,*) '  "888888. 888  888 88888888 .d888888 888  888'
    write (stderr,*) '       888 888  888 88b.     888  888 888 d88P'
    write (stderr,*) '   888888" 888  888  "88888  "8888888 88888P" '
    write (stderr,*) '                                      888     '
    write (stderr,*) '                                      888     '
    write (stderr,*) '                                      888     '
    write (stderr,*) '                                              '
    write (stderr,*) '   Authors: Ben Shires (bs511@cam.ac.uk)      '
    write (stderr,*) '            Chris Pickard (cjp20@cam.ac.uk)   '
    write (stderr,*) '                                              '
    write (stderr,*) '                 2019-2022 (c)                '
    write (stderr,*) '                                              '

  end subroutine banner

end program sheap
