!------------------------------------------------------------------------------
!
! `tools` Source File
!
! tools.f90 source function file. This file contains subroutines for reboxing
! individual atoms according to periodic boundary conditions, as well as
! freezing atoms that are defined as part of the simulation edge, or boundary.
!
!------------------------------------------------------------------------------

MODULE tools

  USE prms
  USE parallel
  USE par

  IMPLICIT none

  integer,dimension(3)                  ::  Nlat ! number of lattice rungs
  real, allocatable, dimension(:)       ::  Xlat ! x1 of each lattice rung
  real, allocatable, dimension(:)       ::  Ylat ! y1 of each lattice rung
  real, allocatable, dimension(:)       ::  Zlat ! z1 of each lattice rung
  integer,allocatable,dimension(:,:) ::  Lijk ! x lattice index of each atom

CONTAINS

  SUBROUTINE rebox(X)
    real, dimension(Natm,3) :: X

    integer :: l

! Put wandering atoms back in box
    do l = 1,3
       X(:,l) = MOD(X(:,l) + 2.*Lb(l),Lb(l))
    end do

  END SUBROUTINE rebox

  SUBROUTINE frzrow1(X,Xi,V,mass)
    real, dimension(Natm,3) :: X,Xi,V
    real, dimension(Natm)   :: mass
    integer    :: i,ii

    if (freeze.eq.1) then
       do ii = 1,Nl
          i = il(ii)
          if (Lijk(i,1).eq.1.or.Lijk(i,1).eq.Nlat(1)) then
             X(i,:) = Xi(i,:)
             V(i,:) = 0.  !SQRT(kB*0.5*(Ttar1+Ttar2)/mass(i))
          end if
       end do
    end if

  END SUBROUTINE frzrow1

  SUBROUTINE freezbottom (X, Xi, V)
    real, dimension(Natm,3) :: X,Xi,V
    integer    :: i,ii

    do ii = 1,Nl
        i = il(ii)
        if (Lijk(i,3).eq.1) then
            X(i,:) = Xi(i,:)
            V(i,:) = 0.
        end if
    end do

    !do i = 1,Natm
    !    if (Lijk(i,3).eq.1) then
    !        !write(31,*)X(i,:),Xi(i,:),i
    !        X(i,:) = Xi(i,:)
    !        V(i,:) = 0.
    !    end if
    !end do

  END SUBROUTINE freezbottom


  SUBROUTINE freezsides (X, Xi, V)

    real, dimension(Natm,3) :: X,Xi,V
    integer    :: i,ii, xmax, ymax, ierr

    if(myid .eq. 0) then
        xmax = MAXVAL(Lijk(:,1))
        ymax = MAXVal(Lijk(:,2))
    end if
    call MPI_BCAST(xmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ymax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    !print*,"xmax and ymax, myid", xmax, ymax,myid
    do ii = 1,Nl
        i = il(ii)
        if (Lijk(i,1).eq.1 .or. Lijk(i,2).eq.1  .or.  Lijk(i,1).eq.xmax .or. Lijk(i,2).eq.ymax) then
            X(i,:) = Xi(i,:)
            V(i,:) = 0.
        end if

    end do

  END SUBROUTINE FREEZSIDES


  SUBROUTINE freezsidesN (X, Xi, V)

    real, dimension(Natm,3) :: X,Xi,V
    integer    :: i,ii,j, xmax, ymax, ierr

    do ii = 1,Nl
        i = il(ii)
        if(xsides(i).gt.0) then
            X(i,:) = Xi(i,:)
            V(i,:) = 0.
        end if
    end do

  END SUBROUTINE freezsidesN

  SUBROUTINE flatrow1(X,V)
    real, dimension(Natm,3) :: X,V
    real, dimension(6)      :: xvm,xvm_g
    integer                 :: N1,NN
    integer    :: i,ii
    integer    :: ierr

    xvm = 0.

    if (freeze.eq.1) then
       do ii = 1,Nl
          i = il(ii)
          if (Lijk(i,1).eq.1) then
             xvm(1) = xvm(1) + X(i,1)
             xvm(2) = xvm(2) + V(i,1)
             xvm(5) = xvm(5) + 1.
          else if (Lijk(i,1).eq.Nlat(1)) then
             xvm(3) = xvm(3) + X(i,1)
             xvm(4) = xvm(4) + V(i,1)
             xvm(6) = xvm(6) + 1.
          end if
       end do

       call MPI_ALLREDUCE(xvm,xvm_g,6, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

       xvm_g(1:2) = xvm_g(1:2)/xvm_g(5)
       xvm_g(3:4) = xvm_g(3:4)/xvm_g(6)

       do ii = 1,Nl
          i = il(ii)
          if (Lijk(i,1).eq.1) then
             X(i,1) = xvm_g(1)
             V(i,1) = xvm_g(2)
          else if (Lijk(i,1).eq.Nlat(1)) then
             X(i,1) = xvm_g(3)
             V(i,1) = xvm_g(4)
          end if
       end do
    end if

  END SUBROUTINE flatrow1

  SUBROUTINE constrainmom(V,mass)
    real, dimension(Natm,3) :: V
    real, dimension(Natm)   :: mass

    real, dimension(3)      :: sw,sw_l
    integer                 :: i,ii
    integer                 :: ierr

    sw_l = 0.
    do ii = 1,Nl
       i = il(ii)
       sw_l = sw_l + mass(i)*V(i,:)
    end do
    call MPI_ALLREDUCE(sw_l,sw,3, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

! enforce zero net momentum
    V(:,1) = V(:,1) - sw(1)/REAL(Natm)/2./mass
    V(:,2) = V(:,2) - sw(2)/REAL(Natm)/2./mass
    V(:,3) = V(:,3) - sw(3)/REAL(Natm)/2./mass

  END SUBROUTINE constrainmom

  SUBROUTINE lattice(X,Xi)
    real, dimension(Natm,3)  :: X,Xi
    !real, allocatable, dimension(:,:) :: X,Xi,Xli_loc
    real, dimension(Natm,3)  :: Xli_loc
    integer                  :: i,j,l,k
    logical                  :: flag

    real                     :: xlt
    real, dimension(3)       :: Xt
    integer, dimension(3)    :: Lt

    real, parameter          :: eps=1.E-13

    if (myid .eq. -1) then
        do i=1,Natm
            write(32,"(6E30.15e3)")X(i,:),Xi(i,:)
        end do
    end if

    allocate (Lijk(Natm,3))

    if (myid.eq.0) write(out_unit,*) "COMPUTING LATTICE"
    do l = 1,3
       Nlat(l) = 1
       Xli_loc(1,l) = X(1,l)
       Lijk(1,l) = 1
       do i = 1,Nsg
          flag = .true.
          do j = 1,i-1
             if (ABS(X(i,l)-Xli_loc(Lijk(j,l),l)).lt.eps) then
                Lijk(i,l) = Lijk(j,l)
                flag = .false.
                exit
             end if
          end do
          if (flag.AND.ABS(X(i,l)-Xli_loc(Nlat(l),l)).gt.eps) then
             Nlat(l) = Nlat(l) + 1
             Xli_loc(Nlat(l),l) = X(i,l)
             Lijk(i,l) = Nlat(l)
          end if
       end do
    end do

    if(allocated(Xlat)) deallocate(Xlat,Ylat,Zlat)
    allocate (Xlat(Nlat(1)),Ylat(Nlat(2)),Zlat(Nlat(3)))

    Xlat = Xli_loc(1:Nlat(1),1)
    Ylat = Xli_loc(1:Nlat(2),2)
    Zlat = Xli_loc(1:Nlat(3),3)

!!$    print *,"*** WARNING: NOT SORTING LATTICE ***"
!!$
    if (myid.eq.0) write(out_unit,*) "SORTING LATTICE"

    do k = 1,Nlat(1)
       do j = 1,Nlat(1)-1
          if (Xlat(j).gt.Xlat(j+1)) then
             xlt = Xlat(j)
             Xlat(j) = Xlat(j+1)
             Xlat(j+1) = xlt
             do i = 1,Natm
                if (Lijk(i,1).eq.j) then
                   Lijk(i,1) = j+1
                else if (Lijk(i,1).eq.j+1) then
                   Lijk(i,1) = j
                end if
             end do
          end if
       end do
    end do

    if(myid.eq. -1) then !taking this portion out. remove if statement to activate.
        do i = 1,Nsg
            do j = 1,Nsg-1
                if (X(j,1).gt.X(j+1,1)) then

                    Xt = X(j,:)
                    X(j,:) = X(j+1,:)
                    X(j+1,:) = Xt

                    Xt = Xi(j,:)
                    Xi(j,:) = Xi(j+1,:)
                    Xi(j+1,:) = Xt

                    Lt = Lijk(j,:)
                    Lijk(j,:) = Lijk(j+1,:)
                    Lijk(j+1,:) = Lt
                end if
            end do
        end do
    end if
    !print*,'sorting done for myid = ', myid
    if (myid.eq.-1) then
       open(lat_unit,file='D/xlat.out')
        do i = 1,Nsg
          write(lat_unit,"(2I7,2E20.10)")i,Lijk(i,1),Xlat(Lijk(i,1)),X(i,1)
       end do
       close(lat_unit)
       open(lat_unit,file='D/ylat.out')
       do i = 1,Nsg
          write(lat_unit,"(2I7,2E20.10)")i,Lijk(i,2),Ylat(Lijk(i,2)),X(i,2)
       end do
       close(lat_unit)
       open(lat_unit,file='D/zlat.out')
       do i = 1,Nsg
          write(lat_unit,"(2I7,2E20.10)")i,Lijk(i,3),Zlat(Lijk(i,3)),X(i,3)
       end do
       close(lat_unit)
   end if

   if (myid .eq. -1) then
       do i=1,Natm
           write(33,"(6E30.15e3)")X(i,:),Xi(i,:)
       end do
   end if


  END SUBROUTINE lattice

END MODULE tools
