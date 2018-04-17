!------------------------------------------------------------------------------
!
! `timeint` Source File
!
! timeint.f90 source function file. This file contains all the subroutines for
! calling necessary force calculations at each timestep, as well as time
! integrating the results to update atomic positions and velocities. Can
! generally be considered the next step down from the top level, main.f90.
!
!------------------------------------------------------------------------------

MODULE timeint

  USE prms
  USE data
  USE out
  USE stats
  USE temp
  USE stillweb
  USE tools
  USE par
  USE knockmod

  IMPLICIT none
  real(8) :: nt_time, at_time, ds_time, nt_time1, at_time1, ds_time1, wt_time, rs_time, f_time1, f_time, time_per_step, ion_time, wtime, wr_time
  integer :: at_count, nt_count, ds_count, f_count, wr_count
  integer :: nrmax, nzmax, vel_size, velstart, globcnt, globcnti, veltscnt
  real    :: density, dr, dz
  real, allocatable, dimension(:,:) :: vel, vel_all, veli, vel_alli, vel_ts, vel_ts_all, vel_tsi, vel_ts_alli
  real, dimension(100000,12) :: ion_data, ion_data_l
  integer :: ion_tstep, ion_tstep_l

CONTAINS

!top level time integration routine
  SUBROUTINE tint

    integer    				:: lt
    real       				:: time,temp, tmprt
    real,dimension(3) 		:: Mom
    real,dimension(Natm,3)	:: F,Fp, W, masslong
    integer 				:: i,j, ierr, whichctrl  !whichctrl=0-vel rescaling, 1-knocking
    logical 				:: test, freerun
    character(30)           :: ion_traj,cmd

!initialize things
    at_time = 0.0; nt_time = 0.0; ds_time = 0.0; f_time=0.0
    at_count = 0; nt_count = 0; ds_count = 0; f_count = 0;
    F = 0. ; Fp = 0.;
    lt = Nt0;    time = restart_time
    masslong(:,1) = mass
    masslong(:,2) = mass
    masslong(:,3) = mass
    test = .false.
    freerun = .false.

    if(myid.eq.1) then
      print *,"Timestep #, Time (ps), Temperature (K)"
    end if

    !if starting fresh, write out initial configuration
    call temperature(V,temp)
    if(lt.eq.0) then
        wr_time = MPI_WTIME()
        call writeatoms(X,V,F,lt,time,temp) !X = X + Ts*V  !call sendposN(X)
        wr_time = MPI_WTIME()-wr_time
    end if

    call si_force(X,F,lt)
    wtime = MPI_WTIME()

    ! testing time it takes for each iteration
    if(test) then
        call TI(X,V,F,masslong,2,Ttar1,Nt,time,lt,temp)
        call calc_time(lt)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr); call MPI_FINALIZE(ierr); stop
    end if

!     initial temperature increase
    if (Nt0 .eq. 0) then
        !initially increase the temperature
        if (myid.eq.0) then
          print *, "Initializing target temperature. . ."
        end if
        Tau = 1.E-14 !accelerated from Tau_i
        call TI(X,V,F,masslong,0,Ttar1,Nt,time,lt,temp)
        call calc_time(lt)
        Tau = Tau_i
    end if

    !freerun enables time integration without ion bombardment
    if(freerun) then
        call ionrun(2)
        call TI(X,V,F,masslong,2,Ttar1,Nt,time,lt,temp)

    !doing impacts, otherwise
    else
        do impact = ions+1,Nlj
            ion_time = MPI_WTIME()

            atom_is_slow = .false.
            atom_is_fast = .true.
            call ionrun(1)
            call knock(impact,lt)
            ion_tstep = 0
            ion_tstep_l = 0
            ion_data = 0.0
            ion_data_l = 0.0

            call TI(X,V,F,masslong,1,Ttar1,Nt,time,lt,temp)
            call writeatoms(X,V,F,lt,time,temp)
            call calc_time(lt)
        end do
    end if

    if (myid.eq.0) print *,"TIME PER STEP: ", (MPI_WTIME()-wtime)/REAL(lt-Nt0)
    if (myid.eq.0) print *,"Total Time in simulation: ", (MPI_WTIME()-wtime)

  END SUBROUTINE tint

!second level time integraiton function; called once per fired ion
!whichctrl determined by objective.
!0: fresh run, initialize system
!1: normal running, ion fired
!2: free running, no check on ion timer
  SUBROUTINE TI(X,V,F,masslong,whichctrl,Ttar,Tstep,time,lt,temp)

    real, dimension(Natm,3)  :: X,V,F,masslong
    real                     :: temp, time, Ttar, T
    integer                  :: whichctrl,Tstep,lt,i


    if(whichctrl .eq. 0) then
        call temperature(V,temp)
        do while (ABS(temp-Ttar) .gt. Teps)
            lt = lt+1
            call onestep(X,V,F,masslong,whichctrl,Ttar,lt,time,temp)
        end do
        if(myid.eq.0) print*, "Equalized to temperature ", temp

    elseif (whichctrl .eq. 1) then
        do while (.not. atom_is_slow)
            lt = lt+1
            if(time .gt. real(impact)*dti) then
                if(myid.eq.0) print*, 'Time reached', time/dti, ' dti. Exiting and firing next one...'
                atom_is_slow = .true.
                exit
            end if
            call onestep(X,V,F,masslong,whichctrl,Ttar,lt,time,temp)
        end do
    elseif (whichctrl .eq. 2) then
        if(myid.eq.0) print*, 'Running the system aplying thermostat ...'
        do i=1,Nt
            lt = lt+1
            call onestep(X,V,F,masslong,whichctrl,Ttar,lt,time,temp)
        end do
    end if
  END SUBROUTINE TI

!bottom level time step, called every step for integration / input output
  SUBROUTINE onestep(X,V,F,masslong,whichctrl,Ttar,lt,time,temp)

    real, dimension(Natm,3)  :: X,V,F,Fp,masslong
    real                     :: temp, time,Ttar
    integer                  :: Tstep,lt,i,whichctrl,ii
    character(30)            :: fn

    time = time+Ts
    restart_time = time

!if initializing, adjust temperature of whole domain
    if(whichctrl .eq. 0) then
        call adjusttemp(V,temp,Ttar)
        call temperature(V,temp)
        !progress report
        if (atm_out.gt.0.and.MOD(lt,atm_out).eq.0) then
          if (myid.eq.0) write(*,"(I10,F10.4, F10.2)")lt,time*1e12,temp
        end if
    end if

!if not initializing, adjust temperature of boundaries
    if(whichctrl .gt. 0) then
        call controlsidetemp (X,V, Ttar)
    end if

    !update positions with velocities and forces
    !X = X + V*Ts + 0.5*F*Ts*Ts/masslong
    do ii=1,Nl
        i=il(ii)
        X(i,:) = X(i,:) + V(i,:)*Ts + 0.5*F(i,:)*Ts*Ts/masslong(i,1)
    end do

    call freezsidesN (X,Xi,V)
    call rebox(X)
    call sendposN(X)

    ! geometrically divides atoms and updates X
    if (mod(lt,atlist).eq.0) then
        at_time1 = MPI_WTIME()
        at_count = at_count+1
        call frzsptrdatms (X,V,F)
        call sendposN(X)
        call initparaops_l(X,V,F)
        call init_nlist
        at_time = at_time + (MPI_WTIME()-at_time1)
    end if

    ! neighborlist redone
    if (mod(lt,ntlist).eq.0) then
        nt_time1 = MPI_WTIME()
        nt_count = nt_count+1
        call si_nlist(X)
        nt_time = nt_time + (MPI_WTIME()-nt_time1)
    end if

    ! determine if moving from fast to slow regime
    if (mod(lt,dslist).eq.0) then
        if (whichctrl .eq. 1) then
            ds_time1 = MPI_WTIME()
            ds_count = ds_count+1
            call calc_atm (lt,time)
            ds_time = ds_time + (MPI_WTIME()-ds_time1)
        end if
    end if

    f_time1 = MPI_WTIME()
    f_count = f_count+1
    !calculate forces at updated positions
    call si_force(X,Fp,lt)
    f_time = f_time + (MPI_WTIME()-f_time1)

    !V = V + 0.5*(F + Fp)*Ts/masslong
    !update velocities by forces
    do ii=1,Nl
        i=il(ii)
        V(i,:) = V(i,:) + 0.5*(F(i,:) + Fp(i,:))*Ts/masslong(i,1)
    end do

    F = Fp

    if (whichctrl .gt. 0) then
      call iostuff(X,V,F,lt,time,temp)
    end if

  END SUBROUTINE onestep


  SUBROUTINE calc_ion (lt,time, F)

    integer :: i,lt
    real :: time
    real, dimension(Natm,3) :: F

    i = Nsg+impact

    if (myid .eq. P(i)) then
        ion_tstep_l = ion_tstep_l+1
        ion_data_l(ion_tstep_l,:) = (/ X(i,:),V(i,:),F(i,:),real(lt),time,real(i)  /)
        !write(ion_unit, "(3E20.10E3, 2I4, 3E20.10E3, 2I9, E15.5)")X(i,:),atype(i),P(i),V(i,:),i,lt, time
    end if

  END SUBROUTINE calc_ion

  SUBROUTINE write_ion

    integer :: i, j, k
    character(30) :: ionfn

    if(ion_tstep_l .gt. 1) then
        write(ionfn, "('G/ion_traj_',I4.4,'_',I4.4,'.dat')") impact, myid
        open(ion_unit, file=ionfn)
        do i=1,ion_tstep_l
            write(ion_unit,"(12E15.6)")ion_data_l(i,:)
        end do
    end if

  END SUBROUTINE write_ion

!function for checking kinetic energy of system, to change timestep size
  SUBROUTINE calc_atm(lt,time)

    integer :: lt,i,ii,ierr
    real :: time
    real, dimension(Natm) :: ke_l,ke

    !calculate kinetic energies of atoms
    ke_l = 0.0
    ke = 0.0
    do ii=1,Nl
        i = il(ii)
        ke_l(i) = 0.5*mass(i)*SUM(V(i,:)**2)/ec
    end do
    call MPI_ALLREDUCE(ke_l, ke, Natm, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    !if maximum kinetic energy is below threshold, switch to slow regime
    if(MAXVAL(ke).lt.ke_limit_1 .and. atom_is_fast) then
        call ionrun(2)
        atom_is_fast = .false.
    end if

  END SUBROUTINE calc_atm


!calculate and print out total time spent on each task
  SUBROUTINE calc_time(lt)

    integer :: lt

    if(myid.eq.0) write(*,"('Total time for ',I6,'th impact is ',F10.3,' seconds')")impact,MPI_WTIME()-ion_time
    if(myid.eq.0) then
        write(*, "('-----------------------------------------------------------------------------')")
        write(*, "('| at_time', E12.4, ' at_count', I9, ' average', E12.4,' |')"),at_time,at_count,at_time/at_count
        write(*, "('| nt_time', E12.4, ' nt_count', I9, ' average', E12.4,' |')"),nt_time,nt_count,nt_time/nt_count
        write(*, "('| ds_time', E12.4, ' ds_count', I9, ' average', E12.4,' |')"),ds_time,ds_count,ds_time/ds_count
        write(*, "('|  f_time', E12.4, '  f_count', I9, ' average', E12.4,' |')"), f_time, f_count, f_time/f_count
        write(*, "('| wr_time = ', E12.4)"),wr_time
        write(*, "('-----------------------------------------------------------------------------')")
    end if

    time_per_step = ((MPI_WTIME()-wtime-at_time-ds_time-nt_time-f_time)/REAL(lt-Nt0))+at_time/at_count+nt_time/nt_count+ds_time/ds_count+f_time/f_count
    if(ds_count .eq. 0) time_per_step = ((MPI_WTIME()-wtime-at_time-nt_time-f_time)/REAL(lt-Nt0))+at_time/at_count+nt_time/nt_count+f_time/f_count

    if (myid.eq.0) write(*, "('at_time = ',F6.2,'% of whole')")(at_time/(at_time+nt_time+ds_time+f_time))*100
    if (myid.eq.0) write(*, "('nt_time = ',F6.2,'% of whole')")(nt_time/(at_time+nt_time+ds_time+f_time))*100
    if (myid.eq.0) write(*, "('ds_time = ',F6.2,'% of whole')")(ds_time/(at_time+nt_time+ds_time+f_time))*100
    if (myid.eq.0) write(*, "(' f_time = ',F6.2,'% of whole')")(f_time/(at_time+nt_time+ds_time+f_time))*100
    if (myid.eq.0) write(*, "('TIME PER STEP:',E12.4)")(MPI_WTIME()-wtime)/REAL(lt-Nt0)
    if (myid.eq.0) write(*, "('Total Time in simulation:',E12.4)")(MPI_WTIME()-wtime)
!    if (myid.eq.0) print *,"TIME PER STEP: ", (MPI_WTIME()-wtime)/REAL(lt-Nt0)
!    if (myid.eq.0) print *,"Total Time in simulation: ", (MPI_WTIME()-wtime)
    at_time = 0.0; nt_time = 0.0; ds_time = 0.0; f_time=0.0
    at_count = 0; nt_count = 0; ds_count = 0; f_count = 0;

  END SUBROUTINE calc_time


END MODULE timeint
