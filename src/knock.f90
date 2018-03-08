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

!knock(knc,lt) assigns the 'knc'th ion fired at timestep lt to its firing location
  SUBROUTINE knock(knc,lt)

    integer            ::knc, ierr,lt
    real            :: knock_vel
    real            :: ran1
    integer         :: ranindex
    real            :: dr, dx, dy

    knock_vel = 1663.6334*SQRT(eV)  !this is for Ga, in m/s

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

    V(Nsg+knc,1) = knock_vel*COS(phixy*Pi/180)*SIN(phiz*Pi/180)
    V(Nsg+knc,2) = knock_vel*SIN(phixy*Pi/180)*SIN(phiz*Pi/180)
    V(Nsg+knc,3) = -1.0*knock_vel*COS(phiz*Pi/180)

    X(Nsg+knc,1) = Xrand(ranindex,1)+Lb(1)/2.0-dx
    X(Nsg+knc,2) = Xrand(ranindex,2)+Lb(2)/2.0-dy
    X(Nsg+knc,3) = knockz*Lb(3)/10.0

    if(myid .eq. 0) then
        write(*,"(I10,4E20.10,I10)")lt,X(Nsg+knc,1:3),V(Nsg+knc,3),Nsg+knc
    end if

  END SUBROUTINE knock

END MODULE knockmod
