!------------------------------------------------------------------------------
!
! `prms` Source File
!
! prms.f90 source function file. This file contains all of the subroutines
! responsible for defining physical constants, reading input parameters from
! sige.in, and broadcasting all of these parameters to the parallel
! processes. Also contains subroutines for creating the initial diamond cubic
! crystal structure for Si and for swapping between the fast/slow timestepping
! and thermostat regimes based on current ion energy.
!
!------------------------------------------------------------------------------

MODULE prms

  USE parallel

  IMPLICIT none

  integer              :: Natm                ! number of atoms
  integer              :: Nlj 				  ! number of lj atoms
  integer			         :: Nsg				  ! number of si-ge atoms
  integer              :: Nrand                 ! number of random points
  integer              :: impact               ! impact number
  integer              :: startimpact          ! starting impact number, only valid if which=0, if 0, start from last, else strat from that impact using restart file.
  integer              :: ions                ! write impact number to restart file
  integer, dimension(3):: Npc                 ! number of pieces in the x,y,z direction for atom division
  integer              :: Nt                  ! number timesteps
  real                 :: fwhm              !full width half max of ion distribution
  real                 :: phiz, phixy         !angles of incidence for ion
  real                 :: Ts                  ! timestep
  real                 :: Ts_i, Ts_r          ! initial and secondary timestep
  real,dimension(3)    :: Lb                  ! box sixe
  real,dimension(3)    :: iLb                 ! inverse box sixe
  real, dimension(4)   :: outz                ! places to collect the out of simulation atoms, in multiples of Lb(3)/10
  real, dimension(2)   :: atomz               ! keep out of sim atoms
  real                 :: ionz, knockz ! place to keep out-of-sim atoms, repository of ions before sim, starting position of ions
  real                 :: sidewidth           ! sidewidth to control temperature in nanometers
  real                 :: ke_limit_2            ! kinetic energy threshhold, for which the atoms are considered slow
  real                 :: ke_limit_1          ! threshold to increase timestep
  logical              :: atom_is_slow, atom_is_fast
  real                 :: restart_time
  integer, allocatable, dimension(:) :: xsides

  real                 :: A                   ! S-W prm
  real                 :: B                   ! S-W prm
  real                 :: psi                 ! S-W prm
  real                 :: qsi                 ! S-W prm
  real                 :: al                  ! S-W prm
  real                 :: lambdaSi            ! S-W prm
  real                 :: lambdaGe            ! S-W prm
  real                 :: gamma               ! S-W prm
  real                 :: epsSi               ! S-W prm
  real                 :: epsGe               ! S-W prm
  real                 :: sigmaSi             ! S-W prm
  real                 :: sigmaGe             ! S-W prm
  real                 :: massSi
  real                 :: massGe
  real                 :: massGa              ! mass of Ga atoms

  real                 :: ZSi                 ! Atomic number Si
  real                 :: ZGe                 ! Atomic number Ge
  real                 :: ZGa                 ! Atomic number Ar
  real                 :: Amol, Bmol, Cmol    ! Moliere prms
  real                 :: alpha_mol           ! Moliere prm
  real                 :: beta_mol            ! Moliere prm
  real                 :: gamma_mol           ! Moliere prm
  real,dimension(3,3)  :: aamol               ! A Moliere
  real,dimension(3)    :: Zmol                ! Z Moliere
  real                 :: ie0                 ! one over eps0 (constant for unit matching)



  real,dimension(2,2)  :: isigma              ! 1./sigma
  real,dimension(2,2)  :: eps                 ! eps
  real,dimension(2,2,2):: leps                ! lambda * eps
  real,dimension(2,2)  :: epssig              ! eps * sigma

  real                 :: Ttar1               ! Target Temperature 1
  real                 :: Qheat               ! Heat flux for #3
  real                 :: Tau                 ! Termostat prm
  real                 :: Tau_i, Tau_r
  integer              :: Tinit               ! # Berendsen on all

  integer              :: whichthrm           ! 1 - Anderson; 2 - Berendsen; 3 - Energy
  integer              :: whichtemp           ! 1 - ends; 2 - left,middle
  integer              :: freeze              ! 1 - freeze ends, 2 - no freeze ends

  integer              :: Nltemp1             ! Temp control # lat left
  integer              :: Nltemp2             ! Temp control # lat right

  real                 :: nu                  ! probability of collision for And Therm.
  real                 :: Vmax                ! Max boltzman velocity for And Therm.

  real                 :: Uo                  ! Initial velocity

  real                 :: Xint                ! location of interface
  integer              :: Nsp                ! N of monolayers in superlattice layer
  real                 :: Amp                 ! forcing amplitude
  real                 :: alpha               ! wave angle
  real                 :: Xofrac              ! postion (fraction of Lb(1))

  integer              :: atm_out             ! atom position output interval
  integer              :: atm_out_i, atm_out_r
  integer              :: abn_out             ! atom pos-vel output interval
  integer              :: res_out             ! restart file output interval
  integer              :: zen_out             ! zone energy change
  integer              :: eng_out             ! energy output
  integer              :: tmp_out             ! mean temperature output

  real                 :: listfac             ! factor rl/rc for neighborlists
  integer              :: ntlist, ntlist_i, ntlist_r              ! how often to reconstruct the neighborlist
  integer              :: atlist, atlist_i, atlist_r              ! how often updates global X
  integer              :: dslist, dslist_i, dslist_r              ! displaced atom update intercval

  integer              :: Ntable              ! number of values is Fsr table

  integer              :: ranseed             ! random number seed
  integer              :: whichic             ! which initial conditions
  integer              :: amorcrys             ! amorphous or crystalline target? 0-cryst, 1-amorphous
  real 				   :: Teps               ! epsilon of temperature



  integer, parameter :: prm_unit  =   4  ! read unit for parameters
  integer, parameter :: out_unit  =   6  ! stdout
  integer, parameter :: atm_unit  =   7  ! atom locations
  integer, parameter :: rdf_unit  =   8  ! radial distribution function
  integer, parameter :: abn_unit  =  10  ! atoms+velocities in binary
  integer, parameter :: lat_unit  =  11  ! lattice definition
  integer, parameter :: res_unit  =  12  ! restart
  integer, parameter :: wav_unit  =  13  ! perturbed position unit
  integer, parameter :: zen_unit  =  14  ! zone energy change
  integer, parameter :: eng_unit  =  15  ! energy diagnostics
  integer, parameter :: tmp_unit  =  16  ! temperature diagnostics
  integer, parameter :: max_unit  =  17  ! maxwell (?) distribution
  integer, parameter :: egk_unit  =  18  ! heat current
  integer, parameter :: fly_unit  =  19  ! y dissipation
  integer, parameter :: yen_unit  =  20  ! y energy output
  integer, parameter :: pbx_unit  =  21  ! x histories
  integer, parameter :: pby_unit  =  22  ! y histories
  integer, parameter :: bce_unit  =  23  ! bc energy
  integer, parameter :: ion_unit  =  24  ! ion trajectory
  integer, parameter :: ds_unit   =  25  ! atoms trajectory
  integer, parameter :: ran_unit  =  26  ! generated random ion positions
  integer, parameter :: snap_unit =  27  ! snapshot xyz file

  real               :: Pi                         ! Pi (duh)
  real               :: srPi                       ! SQRT(Pi)
  real               :: thrd                       ! 1./3.
  real   , parameter :: kB        =   1.381E-23    ! k-Boltzmann (J/K)
  real   , parameter :: eps0      =   8.854187E-12 ! permittivity constant (F/m = C^2/J m)
  real   , parameter :: ec        =   1.60217646E-19 ! elementary charge (C)
  real   , parameter :: ao        =   0.053E-09    ! Bohr radius

  integer, parameter :: Nmaxg=200      ! RDF num bins
  integer, parameter :: Nmaxw=200      ! M-B num bins
  real               :: eV
  real               :: dti            ! time between subsequent impacts

CONTAINS

!readprms() defines material parameters and reads input from siga.in
  SUBROUTINE readprms ()

    integer   :: ierr
    integer :: Nc1, Nc2, Nc3
    integer, dimension (8) :: seed

    if (myid.eq.0) write(out_unit,*) "READING PARAMETERS"

    Pi = 4.*ATAN(1.0)
    srPi = SQRT(Pi)
    thrd = 1./3.
    ZSi = 14.0
    ZGe = 32.0
    ZGa = 31.0
    massSi = 46.637063E-27
    massGe = 120.6049350E-27
    massGa = 1.157777417E-25

    Amol     = 0.35                     ! Moliere prm
    Bmol     = 0.55                     ! Moliere prm
    Cmol     = 0.10                     ! Moliere prm
    alpha_mol= 0.3                      ! Moliere prm
    beta_mol = 1.2                      ! Moliere prm
    gamma_mol= 6.0                      ! Moliere prm

    A     = 7.049556277                 ! S-W prm
    B     = 0.6022245584                ! S-W prm
    psi   = 4.0                         ! S-W prm
    qsi   = 0.0                         ! S-W prm
    al    = 1.8                         ! S-W prm

    lambdaSi = 21.0                     ! S-W prm
    lambdaGe = 31.0						! S-W prm
    gamma = 1.2                     	! S-W prm
    epsSi = 3.473928E-19                ! S-W prm
    epsGe = 3.085E-19                   ! S-W prm
    sigmaSi = 2.0951E-10                ! S-W prm
    sigmaGe = 2.0951E-10                ! S-W prm, made same as Silicon, since we do not have any germanium in our system, 2.1810E-10 for Ge

    nu = 0.01
    Vmax = 3500.0
    Tinit = 2000
    whichthrm = -2
    whichtemp = -1
    freeze = 0
    Qheat = 8.E-23
    Nltemp1 = 15
    Nltemp2 = 15
    listfac = 1.1   ! this is extra padding of neighborlist. prev 1.05, higher value lists more atoms as neighbors, while they will be given 0 force. Kallol

    abn_out = -1
    zen_out = -1
    Ntable = 12000  ! number of division in table, prev 10000 kallol

    if (myid.eq.0) then

        call DATE_AND_TIME(values=seed)
        open(prm_unit,file='siga.in')

        read(prm_unit,*) amorcrys
        read(prm_unit,*) whichic
        read(prm_unit,*) startimpact
        read(prm_unit,*) eV
        read(prm_unit,*) Nc1,Nc2,Nc3
        read(prm_unit,*) Npc(1), Npc(2), Npc(3)
        read(prm_unit,*) Nlj

        read(prm_unit,*) Nt
        read(prm_unit,*) Ts_i
        read(prm_unit,*) Ts_r

        read(prm_unit,*) Ttar1
        read(prm_unit,*) ke_limit_1
        read(prm_unit,*) ke_limit_2

        read(prm_unit,*) Tau_i
        read(prm_unit,*) Tau_r
        read(prm_unit,*) Teps

        read(prm_unit,*) Uo

        read(prm_unit,*) ntlist_i
        read(prm_unit,*) ntlist_r
        read(prm_unit,*) atlist_i
        read(prm_unit,*) atlist_r
        read(prm_unit,*) dslist_i
        read(prm_unit,*) dslist_r

        read(prm_unit,*) atm_out_i
        read(prm_unit,*) atm_out_r
        read(prm_unit,*) res_out
        read(prm_unit,*) eng_out
        read(prm_unit,*) tmp_out

        read(prm_unit,*) Nrand

        read(prm_unit,*) outz(1)
        read(prm_unit,*) outz(2)
        read(prm_unit,*) outz(3)
        read(prm_unit,*) outz(4)

        read(prm_unit,*) atomz(1)
        read(prm_unit,*) atomz(2)
        read(prm_unit,*) ionz
        read(prm_unit,*) knockz
        read(prm_unit,*) sidewidth
        read(prm_unit,*) dti

        read(prm_unit,*) fwhm

        read(prm_unit,*) phiz
        read(prm_unit,*) phixy

        close(prm_unit)

        if (amorcrys .eq. 0) then
            call pureepi(Nc1, Nc2,Nc3) !calculates Lb, Natm, Nsg of crystalline target
            ! disabling prev line adding new lines
            ! also enable xsides!!!!!!!!!!!!!!!!!#####################
            !Lb = 5.431073E-10
            !Lb(1) = Lb(1)*Nc1
            !Lb(2) = Lb(2)*Nc2
            !Lb(3) = Lb(3)*Nc3*10

        else
            Lb = 5.431073E-10
            Lb(1) = Lb(1)*Nc1
            Lb(2) = Lb(2)*Nc2
            Lb(3) = Lb(3)*Nc3*10
            Natm = Nsg+Nlj
        end if

        ranseed = -seed(7)*seed(8)

    end if

    call MPI_BCAST(amorcrys,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(eV,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Natm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nsg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nc1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nc2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nc3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Npc,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nlj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Ts_i,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Ts_r,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Lb,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Ttar1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ke_limit_1,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ke_limit_2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Tau_i,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Tau_r,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Teps,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Uo,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    call MPI_BCAST(ntlist_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ntlist_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atlist_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atlist_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dslist_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dslist_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    call MPI_BCAST(atm_out_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atm_out_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(res_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(eng_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tmp_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    call MPI_BCAST(ranseed,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(whichic,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(startimpact,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nrand,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(outz,4,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atomz,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ionz,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(knockz,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(sidewidth,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dti,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fwhm,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phiz,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phixy,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    call MPI_BCAST(Lb,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    !add these two lines in production run
    if(myid .ne. 0) allocate (xsides(Natm))
    call MPI_BCAST(xsides,Nsg,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    call setupSiGe
    call setupMol
    call ionrun(1)
    iLb = 1./Lb

    !if (myid.eq.Np-1) then
    if(myid .eq. -1) then

       write(out_unit,*)
       write(out_unit,*) " ---- Numerical parameters ----- "
       write(out_unit,*) "Nt      = ",Nt
       write(out_unit,*) "Ts      = ",Ts
       write(out_unit,*) "ntlist  = ",ntlist
       write(out_unit,*) "atlist  = ",atlist
       write(out_unit,*) "listfac = ",listfac
       write(out_unit,*) "ranseed = ",ranseed
       write(out_unit,*) "whichic = ",whichic
       write(out_unit,*) "outz 1-2= ",outz(1), outz(2)
       write(out_unit,*) "outz 3-4= ",outz(3), outz(4)
       write(out_unit,*) "atomz   = ", atomz(1), atomz(2)
       write(out_unit,*) "ionz    = ", ionz
       write(out_unit,*) "knockz  = ", knockz
       write(out_unit,*) "Npc =     ", Npc
       write(out_unit,*) "phiz =     ", phiz
       write(out_unit,*) "phixy =     ", phixy
       write(out_unit,*)
       write(out_unit,*) " ----- Physical parameters ----- "
       write(out_unit,*) "Natm    = ",Natm
       write(out_unit,*) "Nsg     = ",Nsg
       write(out_unit,*) "Nlj     = ",Nlj
       write(out_unit,*) "Nrand   = ",Nrand
       write(out_unit,*) "Uo      = ",Uo
       write(out_unit,*) "Lb      = ",Lb
       write(out_unit,*) "A       = ",A
       write(out_unit,*) "B       = ",B
       write(out_unit,*) "psi     = ",psi
       write(out_unit,*) "qsi     = ",qsi
       write(out_unit,*) "al      = ",al
       write(out_unit,*) "lambdaSi= ",lambdaSi
       write(out_unit,*) "lambdaGe= ",lambdaGe
       write(out_unit,*) "gamma   = ",gamma
       write(out_unit,*) "epsSi   = ",epsSi
       write(out_unit,*) "epsGe   = ",epsGe
       write(out_unit,*) "sigmaSi = ",sigmaSi
       write(out_unit,*) "sigmaGe = ",sigmaGe
       write(out_unit,*) "massSi  = ",massSi
       write(out_unit,*) "massGe  = ",massGe
       write(out_unit,*)
       write(out_unit,*) " --- Temperature parameters ---- "
       write(out_unit,*) "Ttar1   = ",Ttar1
       write(out_unit,*) "Qheat   = ",Qheat
       write(out_unit,*) "Tau     = ",Tau
       write(out_unit,*) "Tinit   = ",Tinit
       write(out_unit,*) "whichthrm = ",whichthrm
       write(out_unit,*) "whichtemp = ",whichtemp
       write(out_unit,*) "freeze  = ",freeze
       write(out_unit,*) "Nltemp1 = ",Nltemp1
       write(out_unit,*) "Nltemp2 = ",Nltemp2
       write(out_unit,*) "nu      = ",nu
       write(out_unit,*) "Vmax    = ",Vmax
       write(out_unit,*)
       write(out_unit,*) " ------ Output parameters ----- "
       write(out_unit,*) "atm_out   = ",atm_out
       write(out_unit,*) "abn_out   = ",atm_out
       write(out_unit,*) "res_out   = ",res_out
       write(out_unit,*) "zen_out   = ",zen_out
       write(out_unit,*) "eng_out   = ",eng_out
       write(out_unit,*) "tmp_out   = ",tmp_out
       write(out_unit,*)
       write(out_unit,*) " ----- Ranging parameters ----- "
       write(out_unit,*) "listfac   = ",listfac
       write(out_unit,*) "Ntlist    = ",Ntlist
       write(out_unit,*) "Ntable    = ",Ntable

    end if

  END SUBROUTINE readprms

!setupMol defines cross interaction coefficients for moliere potential
  SUBROUTINE setupMol

    Zmol(1) = ZSi
    Zmol(2) = ZGe                                                   ! previously    0. zubaer
    Zmol(3) = ZGa

    aamol = 0.
    aamol(1,3) = 0.885*ao*(SQRT(ZSi) + SQRT(ZGa))**(-2./3.)
    aamol(2,3) = 0.885*ao*(SQRT(ZGe) + SQRT(ZGa))**(-2./3.)         ! previously set  aamol(1,3)        zubaer
    aamol(3,1) = aamol(1,3)
    aamol(3,2) = aamol(2,3)
    aamol(3,3) = 0.885*ao*(SQRT(ZGa) + SQRT(ZGa))**(-2./3.)

    ie0 = 1/(4*Pi*eps0)    !  Coulomb force constant
    !  where eps0 is the permittivity of free space

  END SUBROUTINE setupMol

!setupSiGe sets up cross interaction constants for Si-Ge Stillinger-Weber potential
  SUBROUTINE setupSiGe

    eps(1,1) =  epsSi
    eps(1,2) = 0.5*(epsSi+epsGe)
    eps(2,2) =  epsGe
    eps(2,1) = 0.5*(epsSi+epsGe)

    isigma(1,1) = 1./sigmaSi
    isigma(1,2) = 2./(sigmaSi+sigmaGe)
    isigma(2,2) = 1./sigmaGe
    isigma(2,1) = 2./(sigmaSi+sigmaGe)

    leps(1,1,1) = lambdaSi*epsSi
    leps(1,1,2) = 2.*lambdaSi*epsSi/3. + lambdaGe*epsGe/3.
    leps(1,2,1) = 2.*lambdaSi*epsSi/3. + lambdaGe*epsGe/3.
    leps(2,1,1) = 2.*lambdaSi*epsSi/3. + lambdaGe*epsGe/3.
    leps(1,2,2) = lambdaSi*epsSi/3. + 2.*lambdaGe*epsGe/3.
    leps(2,2,1) = lambdaSi*epsSi/3. + 2.*lambdaGe*epsGe/3.
    leps(2,1,2) = lambdaSi*epsSi/3. + 2.*lambdaGe*epsGe/3.
    leps(2,2,2) = lambdaGe*epsGe

    epssig(1,1) = epsSi/sigmaSi
    epssig(1,2) = (epsSi+epsGe)/(sigmaSi+sigmaGe)
    epssig(2,1) = (epsSi+epsGe)/(sigmaSi+sigmaGe)
!    epssig(1,2) = 0.5*(epsSi*sigmaSi+epsGe*sigmaGe)
!    epssig(2,1) = 0.5*(epsSi*sigmaSi+epsGe*sigmaGe)
    epssig(2,2) = epsGe/sigmaGe

  END SUBROUTINE setupSiGe

!pureepi(Nc1,Nc2,Nc3) defines a perfect starting diamond cubic lattice,
!according to unit cell counts Nc1,Nc2,Nc3 in the x,y,z directions respectively
  SUBROUTINE pureepi(Nc1,Nc2,Nc3)
    ! GROWTH DIRECTION is Z
    real,dimension(8,3)                :: uc
    integer                            :: Nc1, Nc2, Nc3, m, mnew
    real, allocatable, dimension(:,:)  :: XX
    integer                            ::  l,i,j,k,cnt, ns
    real                               :: shift, maxx, maxy

    allocate(XX(8*Nc1*Nc2*Nc3,3))

    Lb = 5.431073E-10;
    !Lb(3) = 10E-10;

    ! VALUES r fr MATERIAL ONE with X MLs in [001] DIRECTION
    uc(1,1) = 0.00;         uc(1,2) = 0.00;            uc(1,3) = 0.00
    uc(2,1) = 0.00;         uc(2,2) = 0.50*Lb(2) ;     uc(2,3) = 0.50*Lb(3)
    uc(3,1) = 0.250*Lb(1);  uc(3,2) = 0.250*Lb(2) ;    uc(3,3) = 0.750*Lb(3)
    uc(4,1) = 0.250*Lb(1);  uc(4,2) = 0.750*Lb(2) ;    uc(4,3) = 0.250*Lb(3)
    uc(5,1) = 0.50*Lb(1);   uc(5,2) = 0.0*Lb(2) ;      uc(5,3) = 0.50*Lb(3)
    uc(6,1) = 0.50*Lb(1);   uc(6,2) = 0.50*Lb(2) ;     uc(6,3) = 0.0*Lb(3)
    uc(7,1) = 0.750*Lb(1);  uc(7,2) = 0.250*Lb(2) ;    uc(7,3) = 0.250*Lb(3)
    uc(8,1) = 0.750*Lb(1);  uc(8,2) = 0.750*Lb(2) ;    uc(8,3) = 0.750*Lb(3)

    m = 1
    do i = 1,Nc1
        do j = 1,Nc2
            do k = 1,Nc3
                do l = 1,8
                    XX(m,1) = uc(l,1) + Lb(1)*REAL(i-1)
                    XX(m,2) = uc(l,2) + Lb(2)*REAL(j-1)
                    XX(m,3) = uc(l,3) + Lb(3)*REAL(k-1)
                    m = m+1
                end do
            end do
        end do
    end do

    mnew = m-1
    cnt = 0
    open (2,file = 'si.dat')
    do i = 1,m-1
        write(2,"(3E20.10,I10)")XX(i,:),1
        cnt =cnt+1
    end do
    close(2)

    print *,"Si atoms",m-1
    print*, "Domain size (nm): x,y,z"
    print*,Lb(1)*Nc1*1e9
    print*,Lb(2)*Nc2*1e9
    print*,Lb(3)*Nc3*1e9

    Lb(1) = Lb(1)*Nc1
    Lb(2) = Lb(2)*Nc2
    Lb(3) = Lb(3)*Nc3*10

    Nsg = cnt
    Natm = Nsg+Nlj

    maxx = MAXVAL(XX(:,1))
    maxy = MAXVAL(XX(:,2))

!    open(2, file='si_sides.dat')
    allocate(xsides(Natm))
    xsides = 0
    do i=1,Nsg
        if(XX(i,1).lt.1.0e-10 .or. XX(i,2).lt.1.0e-10 .or. XX(i,1).gt.maxx-1.0e-10 .or. XX(i,2).gt.maxy-1.0e-10) then
            xsides(i) = i
!            write(2,"(I10)")i
        end if
    end do
!    close(2)

  END SUBROUTINE pureepi

!ionrun(tic) adjusts frequency of operations based on whether simulations is
!in fast or slow regimes (by ion energy theshold), tic = 1 or 2 respectively
!controls file write frequency, neighborlist generation frequency, core assignment frequency
!thermostat coefficient, dslist
  SUBROUTINE ionrun(tic)
    integer  :: tic

    if(myid .eq. 0) print*, 'called ionrun', tic

    if(tic .eq. 1) then
        Ts = Ts_i
        atm_out = atm_out_i
        ntlist = ntlist_i
        atlist = atlist_i
        Tau = Tau_i
        dslist = dslist_i
    elseif (tic .eq. 2) then
        Ts = Ts_r
        atm_out = atm_out_r
        ntlist = ntlist_r
        atlist = atlist_r
        Tau = Tau_r
        dslist = dslist_r
    end if

  END SUBROUTINE ionrun

END MODULE prms
