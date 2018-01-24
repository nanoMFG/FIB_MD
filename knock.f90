MODULE knockmod

  USE prms
  USE data

CONTAINS

  SUBROUTINE knock(knc,lt)

    integer            ::knc, ierr,lt
    real            ::theta, knock_vel
    real            :: ran1
    integer         :: ranindex

    knock_vel = 1663.6334*SQRT(eV)  !this is for Ga, prev. 2198.312*SQRT(eV)    ! 1KeV  Knock on velocity (of Ar atom) in m/s
    theta=0.0

    if (myid .eq. 0) then
        ranindex = NINT(ran1(ranseed)*Nrand)
        print*,'             |'
        print*,'           .\:/.'
        print*,'KNOCK    --=>*<=--    RANINDEX', ranindex, 'ion ID', Nsg+knc
        print*,'           //:\\'
    end if

    call MPI_BCAST(ranindex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    !V(Nsg+knc,1) = SIND(theta)*knock_vel
    !V(Nsg+knc,2) = 0.0
    !V(Nsg+knc,3) = -COSD(theta)*knock_vel
    V(Nsg+knc,1) = 0.0
    V(Nsg+knc,2) = 0.0
    V(Nsg+knc,3) = -1.0*knock_vel
    
    X(Nsg+knc,1) = Xrand(ranindex,1)+Lb(1)/2.0
    X(Nsg+knc,2) = Xrand(ranindex,2)+Lb(2)/2.0
    X(Nsg+knc,3) = knockz*Lb(3)/10.0
    
    !write(*,*)X(Nsg+knc,1),X(Nsg+knc,2),X(Nsg+knc,3)

    if(myid .eq. 0) then
        write(*,"(I10,4E20.10,I10)")lt,X(Nsg+knc,1:3),V(Nsg+knc,3),Nsg+knc
        !print*, " KNOCK...the ion is",Nsg+knc, "z velocity is",V(Nsg+knc,3)
        !print*, " The X position is ", X(Nsg+knc,1:3)
    end if
    !print *, "theta equals", theta

  END SUBROUTINE knock

END MODULE knockmod
