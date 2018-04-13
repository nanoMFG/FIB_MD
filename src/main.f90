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
	if (myid.eq.0) then
    print *, "Initparallel done"
  end if
    call readprms
  if (myid.eq.0) then
    print *, "readprms done"
  end if
    call initran
  if (myid.eq.0) then
    print *, "initran done"
  end if
    call initdata
  if (myid.eq.0) then
    print *, "initdata done"
  end if
    call initparaops(X) ! previously Xi
  if (myid.eq.0) then
    print *, "initparaops done"
  end if
    call init_nlist
  if (myid.eq.0) then
    print *, "init_nlist done"
  end if
    call initio
  if (myid.eq.0) then
    print *, "initio done"
  end if
    call mktables
  if (myid.eq.0) then
    print *, "mktables done"
  end if
    call si_nlist(X) ! previously Xi
  if (myid.eq.0) then
    print *,"si_nlist done"
  end if
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
