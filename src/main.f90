!------------------------------------------------------------------------------
!
! `MAIN` Source File
!
! Main source function file. Contains top level calls to functions for
! initialization of data structures, parallel routines, and input/output.
! Also contains the top level function call for the entire MD operation, as
! dictated by the input file's parameters. Ends by closing out all parallel
! and input/output operations.
!
!------------------------------------------------------------------------------

PROGRAM silicon

  USE prms
  USE data
  USE stillweb
  USE timeint
  USE out
  USE par
  USE parallel

  IMPLICIT none

  call init
  call tint
  call close

CONTAINS

  SUBROUTINE init

    call initparallel
	print *, "Initparallel done"
    call readprms
	print *, "readprms done"
    call initran
	print *, "initran done"
    call initdata
	print *, "initdata done"
    call initparaops(X) ! previously Xi
	print *, "initparaops done"
    call init_nlist
	print *, "init_nlist done"
    call initio
	print *, "initio done"
    call mktables
	print *, "mktables done"
    call initmeans
	print *, "initmeans done"
    call si_nlist(X) ! previously Xi
  print *,"si_nlist done"
    if (myid.eq.0) then
       write(out_unit,*) "DONE INITIALIZING"
       write(out_unit,*) "--------------------------------"
    end if

  END SUBROUTINE init

  SUBROUTINE initran

    real    :: r
    real    :: ran1

    !initranions generates fwhm.dat which dictates the random gaussian distribution
    !of initial ion locations, based on parameters from siga.in, prms.f90

    call initranions

  END SUBROUTINE initran

  SUBROUTINE close

    call closeio
    call closeparaops
    call closeparallel

  END SUBROUTINE close

END PROGRAM silicon
