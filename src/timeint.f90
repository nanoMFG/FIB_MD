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

  SUBROUTINE tint

    integer    				:: lt       
    real       				:: time,temp, tmprt
    real,dimension(3) 		:: Mom
    real,dimension(Natm,3)	:: F,Fp, W, masslong
    integer 				:: i,j, ierr, whichctrl  !whichctrl=0-vel rescaling, 1-knocking
    logical 				:: test, freerun
    character(30)           :: ion_traj,cmd

    at_time = 0.0; nt_time = 0.0; ds_time = 0.0; f_time=0.0
    at_count = 0; nt_count = 0; ds_count = 0; f_count = 0;
    F = 0. ; Fp = 0.;
    lt = Nt0;    time = restart_time
    masslong(:,1) = mass
    masslong(:,2) = mass
    masslong(:,3) = mass
    test = .false.
    freerun = .false.
  
    call temperature(V,temp)
    if(lt.eq.0) then
        wr_time = MPI_WTIME()
        call writeatoms(X,V,F,lt,temp) !X = X + Ts*V  !call sendposN(X)
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
    
    ! temperature increase
    !if (Nt0 .eq. 0) then
        !initially increase the temperature
    !    call TI(X,V,F,masslong,0,Ttar1,Nt,time,lt,temp)
    !    call writerestart (X,V,lt)
    !    call calc_time(lt)
    !end if

    ! doing impacts
    if(freerun) then
        call ionrun(2)
        call TI(X,V,F,masslong,2,Ttar1,Nt,time,lt,temp)
        
    else
        do impact = ions+1,Nlj
            write(cmd,"('mkdir -p ',I4.4,'_det')")impact
            call system(cmd)
            write(cmd,"('mkdir -p ',I4.4,'_dim')")impact
            call system(cmd)
            ion_time = MPI_WTIME()
            ! not writing now, will write at the end of impact
            !write(ion_traj,"('G/ion_traj_',I6.6,'.dat')")impact
            !open(ion_unit, file=ion_traj, POSITION='APPEND')
            
            atom_is_slow = .false.
            atom_is_fast = .true.
            moved_index_l = 0
            moved_index = 0
            call ionrun(1)
            call knock(impact,lt)
            ion_tstep = 0
            ion_tstep_l = 0
            ion_data = 0.0
            ion_data_l = 0.0

            call initvel
            call TI(X,V,F,masslong,1,Ttar1,Nt,time,lt,temp)
            call finalizevel
            !call finalizevelts(lt)
            call write_ion
            call cntsputter(X,V,F,lt,impact)
            call writeatoms(X,V,F,lt,temp)
            call writerestart (X,V,lt)
            !close(ion_unit)
            call calc_time(lt)
            ! *************IMPORTANT******************
            ! apply this only for legolas
            !if(mod(impact,2).eq.0)then  
                !call MPI_BARRIER(MPI_COMM_WORLD,ierr); call MPI_FINALIZE(ierr); stop
            !end if
        end do
    end if
    
    if (myid.eq.0) print *,"TIME PER STEP: ", (MPI_WTIME()-wtime)/REAL(lt-Nt0)
    if (myid.eq.0) print *,"Total Time in simulation: ", (MPI_WTIME()-wtime)

  END SUBROUTINE tint


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
            !call onestep(X,V,F,masslong,whichctrl,Ttar,lt,time,temp)
            !if(atom_is_slow .and. myid.eq.0) then
            !    print*,'All the atoms are slowed down, exiting and firing next one ...'
            !end if
            !if(velstart .eq. 100) then !just for testing the velocity profile.
            !    if(myid .eq. 0) print*, 'testing velocity profile, exiting now.... at TI, velstart = ', velstart
            !    atom_is_slow = .true. ! using this condition just to break the loop and fire next one.
            !    exit
            !end if
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


  SUBROUTINE onestep(X,V,F,masslong,whichctrl,Ttar,lt,time,temp)

    real, dimension(Natm,3)  :: X,V,F,Fp,masslong
    real                     :: temp, time,Ttar
    integer                  :: Tstep,lt,i,whichctrl,ii
    character(30)            :: fn

    time = time+Ts
    restart_time = time
    !X_old = X

    if(whichctrl .eq. 0) then
        call adjusttemp(V,temp,Ttar)
        call temperature(V,temp)
    end if
    if(whichctrl .gt. 0) call controlsidetemp (X,V,F, Ttar)
    !call temperature(V,temp)
    
    !X = X + V*Ts + 0.5*F*Ts*Ts/masslong
    do ii=1,Nl
        i=il(ii)
        X(i,:) = X(i,:) + V(i,:)*Ts + 0.5*F(i,:)*Ts*Ts/masslong(i,1)
    end do
    
    !call freezbottom(X,Xi,V)
    call freezsidesN (X,Xi,V)
    call rebox(X)
    call sendposN(X)
    
    ! geometrically divides atoms and updates X
    if (mod(lt,atlist).eq.0) then
        at_time1 = MPI_WTIME()
        at_count = at_count+1
        call frzsptrdatms (X,V,F,lt)
        call sendposN(X)
        !call updateglobal(X)
        !call initparaops(X)
        call initparaops_l(X,V,F)
        call init_nlist
        at_time = at_time + (MPI_WTIME()-at_time1)
        !call frzsptrdatms (X,V,F,lt)
    end if
    
    ! neighborlist redone
    if (mod(lt,ntlist).eq.0) then
        nt_time1 = MPI_WTIME()
        nt_count = nt_count+1
        call si_nlist(X,lt)
        nt_time = nt_time + (MPI_WTIME()-nt_time1)
    end if

    ! write displaced atoms and decide if its time to fire next one
    if (mod(lt,dslist).eq.0) then
        if (whichctrl .eq. 1) then
            ds_time1 = MPI_WTIME()
            ds_count = ds_count+1
            !write(fn, "('E/ds_',I9.9,'.dat')")lt
            !open(ds_unit, file=fn)
            call calc_atm (lt,time)
            ds_time = ds_time + (MPI_WTIME()-ds_time1)
        end if
    end if

    f_time1 = MPI_WTIME()
    f_count = f_count+1
    call si_force(X,Fp,lt)
    f_time = f_time + (MPI_WTIME()-f_time1)

    
    !V = V + 0.5*(F + Fp)*Ts/masslong
    do ii=1,Nl
        i=il(ii)
        V(i,:) = V(i,:) + 0.5*(F(i,:) + Fp(i,:))*Ts/masslong(i,1)
    end do
    
    F = Fp
    call iostuff(X,V,F,lt,time,temp)
    !call diagnostics(X,V,lt,time)
    !X_new = X
    if (whichctrl .eq. 1) then
        call calc_ion (lt,time, F)

        velstart = velstart+1
        veltscnt = veltscnt+1
        call profvel
        if(mod(lt,1000) .eq. 0) then
            call finalizevelts(lt)
            !call midsec(lt)
        end if
        if(Ts.eq.Ts_i .and. mod(lt,500) .eq. 0) then
            call bub_dim (X,V,lt,impact,time)
        elseif (Ts.eq.Ts_r .and. mod(lt,25) .eq. 0) then
            call bub_dim (X,V,lt,impact,time)
        end if
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

  SUBROUTINE calc_atm(lt,time)

    integer :: lt,i,ii,ierr
    real :: time
    real, dimension(Natm) :: ke_l,ke

    ke_l = 0.0
    ke = 0.0
    do ii=1,Nl
        i = il(ii)
        ke_l(i) = 0.5*mass(i)*SUM(V(i,:)**2)/ec
        !if (ke_l(i).gt.ke_limit_2 .or. moved_index(i).gt.0) then
        !    moved_index_l(i) = 1
        !    write(ds_unit,"(3E20.10E3,2I3,4E20.10E3,I8)")X(i,:),atype(i),P(i),V(i,:),ke_l(i)
        !end if
    end do
    call MPI_ALLREDUCE(ke_l, ke, Natm, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    !call MPI_ALLREDUCE(moved_index_l, moved_index, Natm, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    
    !if(myid.eq.0) print*,MAXVAL(ke)

    if(MAXVAL(ke).lt.ke_limit_1 .and. atom_is_fast) then
        call ionrun(2)
        atom_is_fast = .false.
    end if
    ! disabling this for now, because using constant time between impacts
    !if(MAXVAL(ke).lt.ke_limit_2) then
    !    atom_is_slow = .true.
    !end if
    
    !len(1) = Lb(1); len(2) = Lb(2); len(3) = Lb(3)/10.0
    !minlength = MIN(len(1),len(2),len(3))
    !minlength = minlength/2.0
    !ds = X_new-X_old
    !do ii = 1,Nl
    !    i = il(ii)
    !    do l=1,3
    !        if(ds(i,l).gt.minlength) ds(i,l) = len(l)-abs(ds(i,l))
    !    end do
    !    sqds_l(i) = SQRT(SUM(ds(i,:)**2))
    !    if( (sqds_l(i).gt.tol1) .or. (moved_index(i).gt.0) ) then
    !        write(ds_unit,"(3E20.10E3,2I3,3E20.10E3)")X(i,:),atype(i),P(i),V(i,:)
        

  END SUBROUTINE calc_atm


  SUBROUTINE calc_time(lt)

    integer :: lt
    
    if(myid.eq.0) write(*,"('Total time for ',I6,'th impact is ',F10.3,' seconds')")impact,MPI_WTIME()-ion_time
    if(myid.eq.0) then
        write(*, "('-----------------------------------------------------------------------------')")
        write(*, "('| at_time', E20.10, ' at_count', I9, ' average', E20.10,' |')"),at_time,at_count,at_time/at_count
        write(*, "('| nt_time', E20.10, ' nt_count', I9, ' average', E20.10,' |')"),nt_time,nt_count,nt_time/nt_count
        write(*, "('| ds_time', E20.10, ' ds_count', I9, ' average', E20.10,' |')"),ds_time,ds_count,ds_time/ds_count
        write(*, "('|  f_time', E20.10, '  f_count', I9, ' average', E20.10,' |')"), f_time, f_count, f_time/f_count
        write(*, "('| wr_time = ', E20.10)"),wr_time
        write(*, "('-----------------------------------------------------------------------------')")
    end if

    time_per_step = ((MPI_WTIME()-wtime-at_time-ds_time-nt_time-f_time)/REAL(lt-Nt0))+at_time/at_count+nt_time/nt_count+ds_time/ds_count+f_time/f_count
    if(ds_count .eq. 0) time_per_step = ((MPI_WTIME()-wtime-at_time-nt_time-f_time)/REAL(lt-Nt0))+at_time/at_count+nt_time/nt_count+f_time/f_count
    
    if (myid.eq.0) write(*, "('at_time = ',F6.2,'% of whole')")at_time/at_count / time_per_step*100
    if (myid.eq.0) write(*, "('nt_time = ',F6.2,'% of whole')")nt_time/nt_count / time_per_step*100
    if (myid.eq.0) write(*, "('ds_time = ',F6.2,'% of whole')")ds_time/ds_count / time_per_step*100
    if (myid.eq.0) write(*, "(' f_time = ',F6.2,'% of whole')") f_time/f_count / time_per_step*100
    if (myid.eq.0) print *,"TIME PER STEP: ", (MPI_WTIME()-wtime)/REAL(lt-Nt0)
    if (myid.eq.0) print *,"Total Time in simulation: ", (MPI_WTIME()-wtime)
    at_time = 0.0; nt_time = 0.0; ds_time = 0.0; f_time=0.0
    at_count = 0; nt_count = 0; ds_count = 0; f_count = 0;
    
  END SUBROUTINE calc_time


  SUBROUTINE initvel
    
    integer :: i, j, k, cnt

    velstart = 0
    veltscnt = 0
    globcnt = 0
    globcnti = 0
    density = (5.431073E-10**3)/8.0
    dr = 5.431073E-10
    dz = dr*10.0
    nzmax = NINT((Lb(3)/10.0+2.0*dz)/dr) + 1
    nrmax = NINT(Lb(1)/2.0/dr) + 1
    vel_size = nrmax*nzmax
    cnt = 0
    if(.not. allocated(vel)) then
        allocate(vel(vel_size,8), vel_all(vel_size,8))
        allocate(vel_ts(vel_size,8), vel_ts_all(vel_size,8))
        vel = 0.0
        vel_all = 0.0
        vel_ts = 0.0
        vel_ts_all = 0.0
        cnt = 0
        do i=1,nrmax
            do j=1,nzmax
                cnt = cnt+1
                vel(cnt,1) = i*dr
                vel(cnt,2) = j*dr
                vel_ts(cnt,1) = i*dr
                vel_ts(cnt,2) = j*dr
                
            end do
        end do
    else
        vel (:,3:8) = 0.0
        vel_all (:,3:8) = 0.0
        vel_ts (:,3:8) = 0.0
        vel_ts_all (:,3:8) = 0.0
        
    end if
    

  END SUBROUTINE initvel


  SUBROUTINE finalizevel

    integer :: cnt, i, j, k, ierr
    real    :: vol, den, r2, r1
    character(30) :: fn

    vel_all = 0.0
    call MPI_REDUCE(vel,vel_all,vel_size*8,MPI_REAL8,MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if(myid .eq. 0) then
        cnt = 0
        do i=1,nrmax
            do j=1,nzmax
                cnt = cnt+1
                r2 = i*dr
                r1 = (i-1)*dr
                vol = Pi*(r2**2-r1**2)*dr
                !print*, velstart, vel_all(cnt,3), vol/density
                vel_all(cnt,8) = vel_all(cnt,3)/velstart/(vol/density);
                if(vel_all(cnt,3) .gt. 0) then  ! just to make sure that it has some atoms
                    vel_all(cnt,4) = vel_all(cnt,4)/vel_all(cnt,3)*vel_all(cnt,8)
                    vel_all(cnt,5) = vel_all(cnt,5)/vel_all(cnt,3)*vel_all(cnt,8)
                    vel_all(cnt,6) = vel_all(cnt,6)/vel_all(cnt,3)*vel_all(cnt,8)
                    vel_all(cnt,7) = vel_all(cnt,7)/vel_all(cnt,3)
                end if
                
            end do
        end do
        
        write(fn,"('E/velp_all_',I4.4,'.dat')")impact
        open(30, file=fn)
        do i=1,vel_size
            write(30,"(8E20.10E3)") vel(i,1:2), vel_all(i,3:8)
        end do
        close(30)
        
    end if

  END SUBROUTINE finalizevel


  SUBROUTINE profvel

    integer :: i, ii, j, k, nr, nz, dum
    real    :: sx, sy, sr, sz, vr, vt, temp

    do ii=1,Nl
        i=il(ii)
        
        if(X(i,3).gt. (Lb(3)-dz)) then
            sz = X(i,3)-Lb(3)
            sz = sz + dz
        else
            sz = X(i,3) + dz
        end if

        !if(sz .lt. (Lb(3)/10.0*outz(1)+dz) .and. atype(i).lt.2) then
        sx = X(i,1) - Lb(1)/2.0
        sy = X(i,2) - Lb(2)/2.0
        sr = sqrt(sx**2+sy**2)
        nr = FLOOR(sr/dr)
        nz = FLOOR(sz/dr)
        dum = nz + nr*nzmax + 1
        vr = (V(i,1)*sx + V(i,2)*sy)/sr
        vt = (sx*V(i,2) - sy*V(i,1))/sr
        if(sr .eq. 0) then
            vr = 0.0; vt = 0.0
        end if
        temp = mass(i)*sum(V(i,:)**2)/3.0/kb

        !if(sz .lt. (Lb(3)/10.0*outz(1)+dz) .and. atype(i).lt.2) then
            if(nr.lt.nrmax .and. nz.lt.nzmax .and. atype(i).lt.2) then
                !globcnt=globcnt+1
                vel(dum,3) = vel(dum,3)+1.0                         !count of atoms
                vel(dum,4) = vel(dum,4)+vr                          !radial
                vel(dum,5) = vel(dum,5)+vt                          !normal to radial
                vel(dum,6) = vel(dum,6)+V(i,3)                      !z component
                vel(dum,7) = vel(dum,7)+temp                        !temperature

                vel_ts(dum,3) = vel_ts(dum,3)+1.0
                vel_ts(dum,4) = vel_ts(dum,4)+vr
                vel_ts(dum,5) = vel_ts(dum,5)+vt
                vel_ts(dum,6) = vel_ts(dum,6)+V(i,3)
                vel_ts(dum,7) = vel_ts(dum,7)+temp

            end if
            !end if
            
        end do
        
    !if(myid .eq. 0) print *, globcnt
    
  END SUBROUTINE profvel
  
  SUBROUTINE finalizevelts(lt)

    integer :: cnt, i, j, k, ierr, lt
    real    :: vol, den, r2, r1
    character(30) :: fn

    vel_ts_all = 0.0
    call MPI_REDUCE(vel_ts,vel_ts_all,vel_size*8,MPI_REAL8,MPI_SUM, 0, MPI_COMM_WORLD, ierr)


    if(myid .eq. 0) then
        cnt = 0
        do i=1,nrmax
            do j=1,nzmax
                cnt = cnt+1
                r2 = i*dr
                r1 = (i-1)*dr
                vol = Pi*(r2**2-r1**2)*dr
                !print*, veltscnt, vel_ts_all(cnt,3), vol/density
                vel_ts_all(cnt,8) = vel_ts_all(cnt,3)/veltscnt/(vol/density);
                if(vel_ts_all(cnt,3) .gt. 0) then  ! just to make sure that it has some atoms
                    vel_ts_all(cnt,4) = vel_ts_all(cnt,4)/vel_ts_all(cnt,3)*vel_ts_all(cnt,8)
                    vel_ts_all(cnt,5) = vel_ts_all(cnt,5)/vel_ts_all(cnt,3)*vel_ts_all(cnt,8)
                    vel_ts_all(cnt,6) = vel_ts_all(cnt,6)/vel_ts_all(cnt,3)*vel_ts_all(cnt,8)
                    vel_ts_all(cnt,7) = vel_ts_all(cnt,7)/vel_ts_all(cnt,3)
                end if
            end do
        end do

        write(fn,"('F/vel_ts_',I9.9,'.dat')")lt
        open(30, file=fn)
        do i=1,vel_size
            write(30,"(8E20.10E3)") vel_ts(i,1:2), vel_ts_all(i,3:8)
        end do
        close(30)

    end if

    veltscnt = 0
    vel_ts (:,3:8) = 0.0
    vel_ts_all (:,3:8) = 0.0

  END SUBROUTINE finalizevelts

  SUBROUTINE midsec (lt)

    integer :: i, ii, lt, ierr, j, cnt, cntr, cnts, k
    real :: lb2, buff, lb1, lb3
    real, allocatable, dimension(:,:) :: outx, outx_send, outx_recv
    character(30) :: fn
    integer, dimension(MPI_STATUS_SIZE) :: status

    if (allocated(outx)) then
    else
        allocate(outx(Natm,8), outx_send(8,Natm), outx_recv(8,Natm))
        outx = 0.0
        outx_send = 0.0
        outx_recv = 0.0
    end if

    lb1 = Lb(1)
    lb2 = Lb(2)
    lb3 = Lb(3)/10.0
    buff = 5.431073E-10
    cnt = 0
    cntr = 0
    cnts = 0

    do ii=1,Nl
        i=il(ii)
        if(X(i,2).lt.(lb2/2.0+buff) .and. X(i,2).gt.(lb2/2.0-buff)) then
            if(X(i,3).gt.(Lb(3)-lb3) .or. X(i,3).lt.(2.0*lb3)) then
                if(X(i,3).gt.(Lb(3)-lb3)) X(i,3) = X(i,3)-Lb(3)
                cnt = cnt+1
                if(myid .eq. 0) then
                    outx(cnt,1:3) = X(i,:)
                    outx(cnt,4:6) = V(i,:)
                    outx(cnt,7) = massSi*(sum(V(i,:)**2))/3.0/kb
                    outx(cnt,8) = real(i)
                else
                    outx_send(1:3,cnt) = X(i,:)
                    outx_send(4:6,cnt) = V(i,:)
                    outx_send(7,cnt) = massSi*(sum(V(i,:)**2))/3.0/kb
                    outx_send(8,cnt) = real(i)
                end if
            end if
        end if
    end do

    if(myid .ne. 0) then
        call MPI_SEND(cnt,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
        if (cnt .gt. 0) then
            call MPI_SEND(outx_send,cnt*8,MPI_REAL8,0,2,MPI_COMM_WORLD,ierr)
        end if
    else
        do j=1,Np-1
            call MPI_RECV(cntr,1,MPI_INTEGER,j,0,MPI_COMM_WORLD,status,ierr)
            if(cntr .gt. 0) then
                call MPI_RECV(outx_recv,8*cntr,MPI_REAL8,j,2,MPI_COMM_WORLD,status,ierr)
                do k = 1,cntr
                    cnt = cnt+1
                    outx(cnt,:) = outx_recv(:,k)
                end do
            end if
        end do
    end if
    
    
    if(myid.eq.0) then
        
        write(fn,"('out/mid_',I9.9,'.dat')")lt
        open(21, file=fn)
        
        do i=1,cnt
            write(21,"(8E20.10)")outx(i,:)
        end do
        close(21)

    end if
    
    
  END SUBROUTINE midsec
  
END MODULE timeint



