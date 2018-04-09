!------------------------------------------------------------------------------
!
! `data` Source File
!
! data.f90 source file. This contains subroutines for initializing the
! positions, velocities, and masses for atoms in the simulation.
!
!------------------------------------------------------------------------------


MODULE data

  USE prms
  USE tools
  USE parallel

  IMPLICIT none

  real, allocatable, dimension(:,:)   :: X, X_old, X_new     ! atomic positions
  real, allocatable, dimension(:,:)   :: Xrand ! random atomic positions of the ions
  integer, allocatable, dimension(:)  :: sputter_index,sputter_index_l,moved_index,moved_index_l ! index of the sputtered atoms
  real, allocatable, dimension(:,:)   :: Xi    ! initial atomic positions
  real, allocatable, dimension(:,:)   :: V     ! atomic velocities
  real, allocatable, dimension(:)     :: mass,ken  ! atomic mass
  integer, allocatable, dimension(:)  :: atype,pp  ! atom type
  integer                             :: Nt0

CONTAINS

!This function initializes array storage for data stored for each atom
!X: position, V: velocities, mass:mass, atype:integer flag for type of atom
!Xrand: randomized starting ion positions
  SUBROUTINE initdata

    integer   :: i,l,j,jo
    integer   :: ix,iy,iz
    integer   :: Ncs

    real      :: ran1

    allocate (X(Natm,3),Xi(Natm,3),V(Natm,3),mass(Natm),atype(Natm),Xrand(Nrand,2), sputter_index(Natm),moved_index(Natm),sputter_index_l(Natm),moved_index_l(Natm), pp(Natm),ken(Natm))
    !allocate (X_old(Natm,3), X_new(Natm,3))
    sputter_index = 0
    sputter_index_l = 0
    moved_index = 0
    moved_index_l = 0
    call ic

  END SUBROUTINE initdata

!determines which initial conditions to start the simulation with
!4 starts from scratch, 0 starts from a previously generated restart file
!-2 starts from a .dat file containing atomic positions and atom types
  SUBROUTINE ic

    integer   :: i,ierr, reason
    real      :: ran1
    character(100) :: cmd, fn


!for a from-scratch run, set timestep, time, and ion count to 0
    if (whichic.gt.0) then
        Nt0 = 0
        ions = 0
        restart_time = 0.0
    end if

    select case (whichic)
    case (4) ! added by Kallol, from-scratch run
        call initlat(X,Xi,atype)
        call readrandimp(Xrand)
        call assignmass
        if(amorcrys.eq.0) call randvel
        !call constrainmom(V,mass)

    case(0)
        !from a .restart file generated in a previous run
        !call initlat(X,Xi,atype)
        call readlat(Xi,atype)
        call MPI_BCAST(Xi,Natm*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        if(startimpact .eq. 0) then
            open(res_unit,file='D/restart.in_old',form='UNFORMATTED',status='OLD')
        else
            write(cmd,"('ls D/restart.in',I6.6, '* > F/restart_list_',I3.3,'.dat')") startimpact, myid
            call system (cmd)
            write(fn,"('F/restart_list_',I3.3,'.dat')")myid
            open(33, file=fn)
            read(33,"(A28)",IOSTAT=reason)fn
            close (33)
            open(res_unit,file=fn,form='UNFORMATTED',status='OLD')
            if (myid .eq. 0) print*, "-----Restarting from ", startimpact, "th impact with restart file ", trim(fn)," -----"
        end if
        read(res_unit)Nt0,restart_time
        read(res_unit)X,V,mass,atype
        read(res_unit)ions
        close(res_unit)
        Xi(Nsg+1:Natm,:) = X(Nsg+1:Natm,:)
        V(Nsg+ions+1:Natm,:) = 0.0
        call readrandimp(Xrand)
        !call lattice(X,Xi)
        !write(*,"('new_atom 7567',3E20.10)")V(7567,:)

    case(-2)
        call readlat(Xi,atype)
        call MPI_BCAST(Xi,Natm*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        Nt0 = 4684000; ions = 400; restart_time = 0.3482100000E-09;
        open(1,file='D/atoms004684000.dat')
        do i = 1,Natm
            read(1,*)X(i,:),atype(i),pp(i),ken(i), V(i,:)
        end do
        call assignmass
        Xi(Nsg+1:Natm,:) = X(Nsg+1:Natm,:)
        V(Nsg+ions+1:Natm,:) = 0.0
        call readrandimp(Xrand)

    case default
        stop 'bad whichic'
    end select

  END SUBROUTINE ic

!initlat(X,Xi,atyp) reads in initial lattice conditions as X, assigns to Xi
  SUBROUTINE initlat(X,Xi,atype)
    integer, dimension(Natm) :: atype
    real, dimension(Natm,3)     :: X,Xi

    call readlat(X,atype)
    call initIons(X,atype)
    Xi = X
    !call lattice(X,Xi)

  END SUBROUTINE initlat

!readlat(X,atype) reads positions and types from si.dat to X and atype
  SUBROUTINE readlat(X,atype)
    integer, dimension(Natm) :: atype
    real, dimension(Natm,3)     :: X
    integer                  :: ierr
    integer                  :: i
    real                     :: ran1

    if (myid.eq.0) then
        if (amorcrys .eq. 0) then
            open(1,file='si.dat')
            do i = 1,Nsg
                read(1,*)X(i,:),atype(i)
            end do
        else
            open(1,file='5nm_relaxed_done.dat')
            do i = 1,Nsg
                read(1,*)X(i,:),atype(i), V(i,:)
            end do
        end if
        close(1)
    end if
    !call MPI_BCAST(X,Natm*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !call MPI_BCAST(atype,Natm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE readlat

!readrandimp(Xrand) reads randomized ion positions from fwhm.dat
  SUBROUTINE readrandimp(Xrand)
    real, dimension(Nrand,2)     :: Xrand
    integer                  :: ierr
    integer                  :: i
    real                     :: ran1

    if (myid.eq.0) then
        open(1,file='fwhm.dat')
        do i = 1,Nrand
            read(1,*)Xrand(i,:)
        end do
        close(1)
    end if
    call MPI_BCAST(Xrand,Nrand*2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE readrandimp

!randvel assigns initially random velocities to atoms, sets net momentum to zero
  SUBROUTINE randvel
    integer :: i, ierr
    real    :: ran1
    real, dimension(3)	:: mmtm

    !added by Kallol
    if (myid.eq.0) then
        mmtm = 0.0
        if (amorcrys .eq.0) then
            do i = 1,Nsg
                V(i,1) = Uo*(ran1(ranseed)-0.5)
                V(i,2) = Uo*(ran1(ranseed)-0.5)
                V(i,3) = Uo*(ran1(ranseed)-0.5)
                mmtm(:) = mmtm(:)+mass(i)*V(i,:)
            end do
      			mmtm = mmtm/Nsg

      			!subtract average momentum per particle from each particle (to make net zero)
      			do i = 1,Nsg
      				V(i,:) = V(i,:) - mmtm/mass(i)
      			end do

            do i = Nsg+1, Natm
                V(i,:) = 0.0
            end do

        else
            do i = Nsg+1, Natm
                V(i,:) = 0.0
            end do
        end if

    end if
    call MPI_BCAST(V,Natm*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE randvel

!assignmass initializes mass to each atom, according to its atype value
  SUBROUTINE assignmass
    integer :: i

    do i = 1,Natm
        if (atype(i).eq.1) then
            mass(i) = massSi
        else if (atype(i).eq.2) then
            mass(i) = massGe
        else if (atype(i).eq.3) then
            mass(i) = massGa
        end if
    end do

  END SUBROUTINE assignmass

!initIons(X,atype) initializes ion positions in the pre-firing layer, into X
  SUBROUTINE initIons(X,atype)
    integer ,dimension(Natm)	:: atype
    integer                    	:: ii,i
    integer                     :: j, ierr
    real,dimension(Natm,3)      :: X
    real,dimension(3)           :: xxij
    real                        :: ran1
    real                        :: rij

    if(myid .eq. 0) then
        do i=Nsg+1,Natm
            X(i,1) = Lb(1)*ran1(ranseed)
            X(i,2) = Lb(2)*ran1(ranseed)
            X(i,3) = ionz*Lb(3)/10.0
            !V(i,:) = 0.
            atype(i) = 3
        end do

        !  Check Ion positions and reassign if necessary
        do i=Nsg+1,Natm
            do j=Nsg+1,Natm
                if (i .gt. j) then
79                  xxij=X(i,:)-X(j,:)
                    rij = sqrt(xxij(1)**2 + xxij(2)**2)
                    if (rij .lt. sigmaSi) then
                        X(i,1) = Lb(1)*ran1(ranseed)
                        X(i,2) = Lb(2)*ran1(ranseed)
                        goto 79
                    end if
                end if
            end do
        end do
    end if

    call MPI_BCAST(X,Natm*3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atype,Natm,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  END SUBROUTINE initIons

!initranions creates fwhm.dat, a random 2D gaussian distribution for ion impact locations
  SUBROUTINE initranions
    real          :: rion1, rion2, rion3, rion4, ran1
    integer       :: i
    character(30) :: fn

    write(fn,"('fwhm.dat')")
    open(ran_unit,file=fn)

    do i=1,Nrand
      rion1 = ran1(ranseed)
      rion2 = ran1(ranseed)

      !convert to a normal distribution
      !also scale to a fwhm given in nm
      rion3 = sqrt(-2*LOG(rion1))*COS(2*Pi*rion2)*fwhm/2.35*1.E-9
      rion4 = sqrt(-2*LOG(rion2))*SIN(2*Pi*rion1)*fwhm/2.35*1.E-9

      write(ran_unit,"(E12.6,A1,E12.6)")rion3,' ',rion4
    end do

    close(ran_unit)

  END SUBROUTINE initranions


END MODULE data
