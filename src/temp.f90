!------------------------------------------------------------------------------
!
! `temp` Source File
!
! temp.f90 source function file. This file contains subroutines for calculating
! temperatures in the system, as well as subroutines for controlling the
! temperature using thermostats applied to various locations. Also contains
! the subroutine responsible for relocating and freezing atoms that were
! sputtered from the target.
!
!------------------------------------------------------------------------------


MODULE temp

  USE prms
  USE data
  USE tools
  USE stats
  USE parallel
  USE par

  IMPLICIT none

  integer   ::  Ntemp1,Ntemp2,NtempTot1,NtempTot2
  integer,allocatable,dimension(:)   ::  T1list,T2list

CONTAINS

!function for applying a berendsen thermostat to atoms along the side of the domain
! X,V are position,velocity and Ttar is target temperature
  SUBROUTINE controlsidetemp (X,V, Ttar)
    real, dimension(Natm,3)     :: X,V
    real                        :: tempside, tempside_l, s, Ttar, nm
    integer                     :: kk_l, ierr, kk, i, ii
    integer, dimension(Nl)      :: tempside_list

    tempside_list = 0
    kk_l = 0
    tempside_l = 0.0
    tempside = 0.0
    nm = 1.0e-9
    do ii=1,Nl
        i = il(ii)
        if(   ((X(i,1).lt.sidewidth*nm) .or. (X(i,1).gt.(Lb(1)-sidewidth*nm)))  &
        .or. ((X(i,2).lt.sidewidth*nm) .or. (X(i,2).gt.(Lb(2)-sidewidth*nm)))  &
        .and.  ((X(i,3).lt.outz(1)*Lb(3)/10.0) .or. (X(i,3).gt.outz(4)*Lb(3)/10.0))) then
          if( atype(i) < 3 ) then
              tempside_l = tempside_l + mass(i)*(V(i,1)**2   + V(i,2)**2   + V(i,3)**2  )
              kk_l = kk_l+1
              tempside_list(kk_l) = i
          end if
        end if
    end do

    if(kk_l .gt. 0) then
        tempside_l = tempside_l/3./REAL(kk_l)/kB

        s = SQRT(1. + Ts/Tau*(Ttar/tempside_l-1.))

        do ii = 1,kk_l
            i = tempside_list(ii)
            V(i,:) = s*V(i,:)
        end do
    end if

  END SUBROUTINE controlsidetemp

!function for counting sputtered atoms, currently unused
  SUBROUTINE cntsputter(X,lt,impact)

    real, dimension(Natm,3)     :: X
    integer                     :: i,ii,ierr,lt,impact, sputtercnt, sputtercnt_l
    real                        :: buffer

    sputtercnt_l = 0
    sputtercnt = 0
    buffer = 0.7

    do ii=1,Nl
        i = il(ii)
        if(X(i,3) .gt. ((atomz(1)-buffer)*Lb(3)/10.0) .and. X(i,3) .lt. ((atomz(1)+buffer)*Lb(3)/10.0) .and. atype(i).lt.3) then
            sputtercnt_l = sputtercnt_l + 1
        end if
    end do
    call MPI_REDUCE(sputtercnt_l,sputtercnt,1,MPI_INTEGER,MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if(myid .eq. 0) then
        open(91, file='sputtercnt.dat', POSITION='APPEND')
        write(91,"(2I9)") impact, sputtercnt
        close(91)
    end if

  END SUBROUTINE CNTSPUTTER


!function for freezing and relocating atoms that have just been sputtered,
!and refreezing atoms that are currently in the sputter containment layers
  SUBROUTINE frzsptrdatms (X,V,F)
    real, dimension(Natm,3)     :: X,V,F
    integer                     :: i,ii,ierr
    real                        :: buffer

    buffer = 0.3
    do ii=1,Nl
        i = il(ii)
        if(   ((X(i,3).gt.outz(1)*Lb(3)/10.0) .and. (X(i,3).lt.outz(2)*Lb(3)/10.0))   ) then
            X(i,3) = atomz(1)*Lb(3)/10.0
            V(i,:) = 0.0
            F(i,:) = 0.0
        elseif (  ((X(i,3).gt.outz(3)*Lb(3)/10.0) .and. (X(i,3).lt.outz(4)*Lb(3)/10.0))    ) then
            X(i,3) = atomz(2)*Lb(3)/10.0
            V(i,:) = 0.0
            F(i,:) = 0.0
        elseif( ((X(i,3).gt.(atomz(1)-buffer)*Lb(3)/10.0) .and. (X(i,3).lt.(atomz(1)+buffer)*Lb(3)/10.0))  ) then
            V(i,:) = 0.0
            F(i,:) = 0.0
        elseif( ((X(i,3).gt.(atomz(2)-buffer)*Lb(3)/10.0) .and. (X(i,3).lt.(atomz(2)+buffer)*Lb(3)/10.0))  ) then
            V(i,:) = 0.0
            F(i,:) = 0.0
        elseif( ((X(i,3).gt.(ionz-buffer)*Lb(3)/10.0) .and. (X(i,3).lt.(ionz+buffer)*Lb(3)/10.0))  ) then
            V(i,:) = 0.0
            F(i,:) = 0.0
        end if
    end do

  END SUBROUTINE frzsptrdatms

!adjust temperature according to andersen thermostat, currently unused / untested
  SUBROUTINE adjusttempAnd(V,Ntemp,Tlist,Temp,Ttar)
    real, dimension(Natm,3)  :: V
    integer                  :: Ntemp
    integer,dimension(Ntemp) :: Tlist
    real                     :: Temp
    real                     :: Ttar
    real                     :: s
    integer                  :: i,ii
    integer                  :: l
    real, dimension(3)       :: Vb
    real                     :: ran1
    real   :: PE,Bltz

    do i = 1,Ntemp
       if (ran1(ranseed).lt.nu) then
          ii = Tlist(i)
          do l = 1,3
             do
                Vb(l) = -Vmax + 2.*Vmax*ran1(ranseed)
                Bltz =EXP(-mass(ii)*Vb(l)*Vb(l)/(2.*kB*Ttar))
                if (ran1(ranseed).lt.Bltz) exit
             end do
          end do
          V(ii,:) = Vb
       end if
    end do
  END SUBROUTINE adjusttempAnd


!function for berendsen thermostat to entire domain, for initialization
  SUBROUTINE adjusttemp(V,Temp,Ttar)
    real, dimension(Natm,3)  :: V

    real                     :: Temp
    real                     :: Ttar

    real                     :: s

    integer                  :: i,ii

    s = SQRT(1. + Ts/Tau*(Ttar/Temp-1.))

    !print *,myid,": ",Ttar,Temp,s
    !if (myid .eq. 0) print*,s, Temp
    do ii = 1,Nl
        i = il(ii)
        if(atype(i).lt.3)  V(i,:) = s*V(i,:)
    end do

  END SUBROUTINE adjusttemp


!  Compute the mean temperature for the system, not including Ga ions
  SUBROUTINE temperature(V,Temp)
    real, dimension(Natm,3)  :: V
    real                     :: Temp,Temp_l
    integer                  :: i,ii
    integer                  :: ierr

!                       3
!       3/2 kB T   =   SUM ( 1/2 m v_i^2 )
!                      i=1

    Temp_l = 0.
    do ii = 1,Nl
        i = il(ii)
        if(atype(i).lt.3) then
            Temp_l = Temp_l + mass(i)*(V(i,1)**2   + V(i,2)**2   + V(i,3)**2  )
        end if
    end do
    Temp_l = Temp_l/3./REAL(Nsg)/kB
    !Temp_l = Temp_l/3./REAL(Nl)/kB  !Kallol

    call MPI_ALLREDUCE(Temp_l, Temp, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  END SUBROUTINE temperature

END MODULE temp
