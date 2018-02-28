!------------------------------------------------------------------------------
!
! `knock` Source File
!
! knock.f90 source function file. This file contains the subroutine for
! placing the ion currently being fired at the surface within its random
! position distribution. It also assigns the ion's velocity vector, according
! to the assigned beam energy and angle of incidence.
!
!------------------------------------------------------------------------------


MODULE knockmod

  USE prms
  USE data

CONTAINS

  SUBROUTINE knock(knc,lt)

    integer            ::knc, ierr,lt
    real            :: knock_vel
    real            :: ran1
    integer         :: ranindex
    real            :: dr, dx, dy

    knock_vel = 1663.6334*SQRT(eV)  !this is for Ga, prev. 2198.312*SQRT(eV)    ! 1KeV  Knock on velocity (of Ar atom) in m/s

    dr = TAN(phiz*Pi/180)*(knockz*Lb(3)/10 - Lb(3)/10)
    dx = dr*COS(phixy*Pi/180)*SIN(phiz*Pi/180)
    dy = dr*SIN(phixy*Pi/180)*SIN(phiz*Pi/180)


    if (myid .eq. 0) then
        ranindex = NINT(ran1(ranseed)*Nrand)
        print*,'             |'
        print*,'           .\:/.'
        print*,'KNOCK    --=>*<=--    RANINDEX', ranindex, 'ion ID', Nsg+knc
        print*,'           //:\\'
    end if

    call MPI_BCAST(ranindex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(dx,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(dy,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    V(Nsg+knc,1) = knock_vel*COS(phixy*Pi/180)*SIN(phiz*Pi/180)
    V(Nsg+knc,2) = knock_vel*SIN(phixy*Pi/180)*SIN(phiz*Pi/180)
    V(Nsg+knc,3) = -1.0*knock_vel*COS(phiz*Pi/180)

    X(Nsg+knc,1) = Xrand(ranindex,1)+Lb(1)/2.0-dx
    X(Nsg+knc,2) = Xrand(ranindex,2)+Lb(2)/2.0-dy
    X(Nsg+knc,3) = knockz*Lb(3)/10.0

    if(myid .eq. 0) then
        write(*,"(I10,4E20.10,I10)")lt,X(Nsg+knc,1:3),V(Nsg+knc,3),Nsg+knc
        !print*, " KNOCK...the ion is",Nsg+knc, "z velocity is",V(Nsg+knc,3)
        !print*, " The X position is ", X(Nsg+knc,1:3)
    end if
    !print *, "theta equals", theta

  END SUBROUTINE knock

END MODULE knockmod
