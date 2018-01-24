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
    !call writeatoms(X,0)
    call mktables
	print *, "mktables done"
    !call TempList(X)
    call initmeans
	print *, "initmeans done"
    call si_nlist(X,0) ! previously Xi
    print *,"si_nlist done"
    if (myid.eq.0) then
       write(out_unit,*) "DONE INITIALIZING"
       write(out_unit,*) "--------------------------------"
    end if

  END SUBROUTINE init

  SUBROUTINE initran
    
    real    :: r
    real    :: ran1

    !if (myid.eq.0) write(out_unit,*) "INTIALIZING RANDOM"
    r = ran1(ranseed)

  END SUBROUTINE initran

  SUBROUTINE close 

    call closeio
    call closeparaops
    call closeparallel

  END SUBROUTINE close

END PROGRAM silicon
