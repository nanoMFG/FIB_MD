!------------------------------------------------------------------------------
!
! `nlist` Source File
!
! nlist.f90 source function file. This contains subroutines for the for defining
! pointers for the neighborlist, and for determining the spatial `cell` each
! atom is located in.
!
!------------------------------------------------------------------------------

MODULE nlistmod

  USE prms
  USE data

  IMPLICIT none

CONTAINS

  SUBROUTINE atomcell(X,Icell,hspace,Nmax)
    real,dimension(Natm,3)       :: X        ! particle positions
    integer,dimension(Natm,3)    :: Icell    ! atom cell list
    real,dimension(3)            :: hspace   ! mesh cell spacing for distribution
    integer,dimension(3)         :: Nmax     ! max cell
    integer                      :: l
    integer i

    do l = 1,3
       Icell(:,l) = INT( ( X(:,l) - FLOOR( X(:,l)*iLb(l) )*Lb(l) )/hspace(l) ) + 1
    end do

  END SUBROUTINE atomcell

  SUBROUTINE chainlist(LL,HOC,N,X)
    integer,dimension(Natm)                          :: LL    ! next in list pointer
    integer,dimension(3)                             :: N     ! either Ncsr or Nclj (# cells in each dimension)
    integer,dimension(0:N(1)+1,0:N(2)+1,0:N(3)+1)    :: HOC   ! Head-of-Chain pointer
                                                              ! with perioidic continuation
    real   ,dimension(Natm,3)                        :: X     ! particle positions
    integer,dimension(Natm,3)                        :: Icell ! atom cell list
    integer                                          :: i
    real,dimension(3)                                :: HMf   ! cell spacing

    HMf = Lb/REAL(N)
    call atomcell(X,Icell,HMf,N)

    HOC = 0
    do i = 1,Natm
       LL(i) = HOC(Icell(i,1),Icell(i,2),Icell(i,3))
       HOC(Icell(i,1),Icell(i,2),Icell(i,3)) = i
   end do
   !print*, LL(2222)
   ! added by Kallo, checking value of LL
   !open(29,file = 'll.dat')
   !do i = 1,Natm
   !    write(29,"(I10)")LL(i)
   !end do
   !close(29)

    if (N(1).gt.2) then
       HOC(0,:,:)        = HOC(N(1),:,:)
       HOC(N(1)+1,:,:)   = HOC(1,:,:)
    end if
    if (N(2).gt.2) then
       HOC(:,0,:)        = HOC(:,N(2),:)
       HOC(:,N(2)+1,:)   = HOC(:,1,:)
    end if
    if (N(3).gt.2) then
       HOC(:,:,0)        = HOC(:,:,N(3))
       HOC(:,:,N(3)+1)   = HOC(:,:,1)
    end if

  END SUBROUTINE chainlist

END MODULE nlistmod
