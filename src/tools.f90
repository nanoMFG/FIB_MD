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

CONTAINS

  SUBROUTINE rebox(X)
    real, dimension(Natm,3) :: X

    integer :: l

! Put wandering atoms back in box
    do l = 1,3
       X(:,l) = MOD(X(:,l) + 2.*Lb(l),Lb(l))
    end do

  END SUBROUTINE rebox

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

END MODULE tools
