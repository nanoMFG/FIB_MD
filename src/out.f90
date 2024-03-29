!------------------------------------------------------------------------------
!
! `OUT` Source File
!
! out.f90 source function file. This contains all of the subroutines used for
! initializing data outputs and then actually writing output data when called
! during the simulation.
!
!------------------------------------------------------------------------------

MODULE out

  USE prms
  USE stats
  USE temp
  USE data
  USE parallel
  USE par


  IMPLICIT none

  integer                         :: NYen
  real,allocatable,dimension(:,:) :: YenS

  integer                         :: NmaxS
  real,allocatable,dimension(:,:) :: MaxwS

CONTAINS

!writerestart(X,V,lt) writes a restart file containing atom positions X,
!velocities V, and current timestep lt
!requires output folder 'D/'
  SUBROUTINE writerestart(X,V,lt)
    real, dimension(Natm,3) :: X,V
    integer                 :: lt
    real                    :: time
    character(30)           :: fn

    call updateroot(X)
    call updateroot(V)

    if (myid.eq.0) then
        write(fn,"('D/restart.in_old')")
        open(res_unit,file=fn,form='UNFORMATTED')
        write(res_unit)lt,restart_time
        write(res_unit)X,V,mass,atype
        write(res_unit)impact
        close(res_unit)

        write(fn,"('D/restart.in',I6.6,'_',I9.9)")impact,lt
        open(res_unit,file=fn,form='UNFORMATTED')
        write(res_unit)lt,restart_time
        write(res_unit)X,V,mass,atype
        write(res_unit)impact
        close(res_unit)
    end if

  END SUBROUTINE writerestart

!writeatoms(X,V,F,lt,time,temp) writes current atomic values to .xyz file, also
!writes basic tracking info to command line
  SUBROUTINE writeatoms(X,V,F,lt,time,temp)
    real, dimension(Natm,3) :: X,V,F
    real, dimension(Natm) :: kin_eng
    integer  :: lt, fon
    integer  :: i,k,l,j,ii
    integer :: ierr
    real    :: temp, TE, time, coef
    character(30)           :: fn
    character(2)            :: tempatype

      if (myid.eq.0) then		!added by Josh, output in xyz format appended to atm_unit
        write(atm_unit,"(I9)"),Natm !heads the section with # of atoms
!        write(atm_unit,"(1A12,F12.5,3A12,2A8,3A16,A8,F14.5)")'Time (ps) = ',time*1e12,'x (A)','y (A)','z (A)','type','proc#','Vx','Vy','Vz'
        write(atm_unit,"(1A12,F12.5,3A12,2A8,3A16,A8,F14.5)")'Time (ps) = ',time*1e12,'x (A)','y (A)','z (A)','type','proc#','KE (eV)'
        open(snap_unit,file='data/3dmol2.xyz',STATUS="REPLACE")
        write(snap_unit,"(I9)"),Natm
        write(snap_unit,"(5A12,F12.5)")'type','x (A)','y (A)','z (A)','Time (ps) = ',time*1e12

        coef = 6.252e18 !conversion from J to eV
        do i = 1,Natm
      	    if(atype(i).eq.1) then
			          tempatype = "Si"
            elseif(atype(i).eq.3) then
		            tempatype = "Ga"
            end if

            TE = 0.5*mass(i)*(V(i,1)**2 + V(i,2)**2 + V(i,3)**2)*coef
        	  write(atm_unit,"(3F12.4,A8,I8,1ES14.3)")X(i,:)*1e10,tempatype,P(i),TE
            write(snap_unit,"(A8,3F10.4)")tempatype,X(i,:)*1e10

        end do
        close(snap_unit)
        !write a line tracking the current system temperature and step number
    		write(*,"(I10,F10.4,F10.2)")lt,time*1e12,temp

      end if

  END SUBROUTINE writeatoms

!initio initializes definitions for intput/output files
  SUBROUTINE initio
    character(10)  fn

!Nt0 is the starting timestep, 0 or set by restart file
    write(fn,"(I10.10)")Nt0

    if (myid.eq.0) then
       if (zen_out.gt.0) open(zen_unit,file='G/zen.out'//fn)
       if (tmp_out.gt.0) open(tmp_unit,file='G/tmp.out'//fn)
       if (eng_out.gt.0) open(eng_unit,file='G/eng.out'//fn)
       if (atm_out.gt.0) open(atm_unit,file='data/mdrun2.xyz')
    end if

  END SUBROUTINE initio

!closeio closes any output files that were opened
  SUBROUTINE closeio
    if (myid.eq.0) then
       if (tmp_out.gt.0) close(tmp_unit)
       if (eng_out.gt.0) close(eng_unit)
       if (zen_out.gt.0) close(zen_unit)
       if (atm_out.gt.0) close(atm_unit)
    end if

  END SUBROUTINE closeio

!iostuff(X,V,F,lt,time,temp) is general function call to check if data should
!be output at any given timestep
  SUBROUTINE iostuff(X,V,F,lt,time,temp)
    real, dimension(Natm,3)  :: X,V,F
    real, dimension(Natm)    :: q
    real                     :: time,temp
    integer                  :: lt

    if (res_out.gt.0.and.MOD(lt,res_out).eq.0.or.   &
         atm_out.gt.0.and.MOD(lt,atm_out).eq.0)  then
          call updateroot(X)
          call updateroot(V)
          call temperature(V,temp)
          if (res_out.gt.0.and.MOD(lt,res_out).eq.0) then
              call writerestart(X,V,lt)
          end if
          if (atm_out.gt.0.and.MOD(lt,atm_out).eq.0) then
              call writeatoms(X,V,F,lt,time,temp)
          end if
    end if

  END SUBROUTINE iostuff

!diagnostics(X,V,lt,time) outputs energies
  SUBROUTINE diagnostics(X,V,lt,time)
    real, dimension(Natm,3)  :: X,V
    real, dimension(Natm)    :: q
    real                     :: time
    integer                     :: lt

    real,dimension(Nlat(1),5)    :: Yen

    real                      :: T,KE,PE           ! Temperatures + Energies
    real,dimension(3)         :: Pm                ! Momentum
    real,dimension(0:Nmaxw,3) :: Maxw              ! Maxwellian distributions
    integer,dimension(Natm)   :: Lall              ! list of all atoms
    real                      :: e1,e2
    real,dimension(0:Nmaxg)   :: g, rg             ! Radial distribution fct
    character(30)             :: fn
    integer                   :: ib,i,j

!instantaneous temperature output
    if (tmp_out.gt.0.and.MOD(lt,tmp_out).eq.0) then
       call temperature(V,T)
       if (myid.eq.0) write(tmp_unit,"(I11, E20.10,F30.10)") lt, time,T
    end if

!potential and kinetic energy output
    if (eng_out.gt.0.and.MOD(lt,eng_out).eq.0) then
       call kineticenergy(V,KE)
       call potentialenergy(X,PE)
       call momentum(V,Pm)
       KE = KE/epsSi/REAL(Natm)
       PE = PE/epsSi/REAL(Natm)
       write(eng_unit,"(E20.10,3F25.15,3E10.2)") time,KE,PE,KE+PE,Pm(:)
    end if

  END SUBROUTINE diagnostics

END MODULE out
