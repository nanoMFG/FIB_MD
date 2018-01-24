program upsidedown

  include 'mpif.h'
  
  integer :: i, j, k, natm, fi, fo, df, nzmax, nrmax, cnt, nr, nz, dum, ff, ii, io, reason, fcount, fcount_all
  integer :: Np, myid, ierr, vel_size, fic
  real, allocatable, dimension(:,:) :: X, V, F, vel, vel_all
  real, allocatable, dimension(:) :: Lb, kin_eng
  integer, allocatable, dimension(:) :: atype, p
  real :: dr, r, dz, vr, kb, ev, vt, r2, r1, z, vol, n, density, den, pi
  character (100) :: cmd
  character (30)  :: fn

  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, Np, ierr )
  !print *, 'Process ', myid, ' of ', Np, ' is alive'
  
  natm = 401000
  allocate(X(natm,3), V(natm,3), F(natm,3), kin_eng(natm), Lb(3))
  allocate(atype(natm),p(natm))
  Lb = (/ 2.7155365E-8, 2.7155365E-8, 1.0862146E-8 /)
  kb = 1.3806503E-23
  ev = 1.60217646E-19
  density = (5.431073E-10**3)/8.0
  pi = atan(1.0)*4.0

  ii = 1
  io = ii+399
  
  write(cmd,"('ls D/restart.in',I6.6, '* > restart_list',I3.3,'.dat')")ii, myid
  call system (cmd)
  write(fn,"('restart_list',I3.3,'.dat')")myid
  open(33, file=fn)
  read(33,"(A28)",IOSTAT=reason)fn
  close (33)
  read(fn(20:28),"(I9)")fi

  write(cmd,"('ls D/restart.in',I6.6, '* > restart_list',I3.3,'.dat')")io, myid
  call system (cmd)
  write(fn,"('restart_list',I3.3,'.dat')")myid
  open(33, file=fn)
  read(33,"(A28)",IOSTAT=reason)fn
  close (33)
  read(fn(20:28),"(I9)")fo

  df = 4000
  fi = NINT(real(fi)/real(df))*df
  fo = NINT(real(fo)/real(df))*df
  ! putting the values of fi fo below to calculate the velocity of flow after bombardment
  fi = 8000
  fo = 4304000
  fic = (fo-fi)/df
  !ii = 999
  fic = NINT(real(fic)/real(Np))
  
  if (myid .eq. 0) then
      fi = fi
      fo = fi + fic*df*(myid+1)
  elseif (myid .eq. (Np-1)) then
      fi = fi + fic*df*myid + df
      fo = fo
  else
      fo = fi + fic*df*(myid+1)
      fi = fi + fic*df*myid + df
  end if

  print*, "Proc = ", myid, "from", fi, "to", fo

  !call MPI_FINALIZE(ierr)
  !stop
  !print*, fi, fo
  
  dr = 0.5E-9
  dz = 5.0E-9
  nzmax = NINT((2.0E-8+dz)/dr) + 1
  nrmax = NINT(Lb(1)/2.0/dr) + 1
  vel_size = nrmax*nzmax
  allocate(vel(vel_size,8), vel_all(vel_size,8))
  vel = 0.0
  vel_all = 0.0

  cnt = 0
  do i=1,nrmax
      do j=1,nzmax
          cnt = cnt+1
          vel(cnt,1) = i*dr
          vel(cnt,2) = j*dr
      end do
  end do

  fcount = 0
  do ff=fi,fo,df
      fcount = fcount+1
      write(fn,"('D/atoms',I9.9,'.dat')")ff
      print*, fn
      open(29, file=fn)
      do i=1,natm
          read(29,*) X(i,:),atype(i),P(i),kin_eng(i), V(i,:)

          if(X(i,3).gt. 8.0E-8) X(i,3)= X(i,3)-5.431073E-9*20
          if(X(i,3).lt. 2.0E-8 .and. atype(i).lt.2) then
              X(i,3) = X(i,3) + dz
              X(i,1) = X(i,1) - Lb(1)/2.0
              X(i,2) = X(i,2) - Lb(2)/2.0
              r = sqrt(X(i,1)**2+X(i,2)**2)
              nr = FLOOR(r/dr)
              nz = FLOOR(X(i,3)/dr)
              dum = nz + nr*nzmax + 1
              vr = sum(V(i,1:2)*X(i,1:2))/r
              vt = (X(i,1)*V(i,2)-X(i,2)*V(i,1))/r
              if(r .eq. 0) then
                  vr = 0.0; vt = 0.0
              end if
                 
              if(nr.lt.nrmax .and. nz.lt.nzmax) then
                  vel(dum,3) = vel(dum,3)+1.0                         !count of atoms
                  vel(dum,4) = vel(dum,4)+vr                          !radial
                  vel(dum,5) = vel(dum,5)+vt                          !normal to radial
                  vel(dum,6) = vel(dum,6)+V(i,3)                      !z component
                  vel(dum,7) = vel(dum,7)+kin_eng(i)*ev/kb*2.0/3.0    !temperature 
              end if

          end if
      end do
      close(29)
  end do

  !print*, "reached 1"
  !call MPI_REDUCE(vel,vel_all,vel_size*8,MPI_REAL,MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  !call MPI_REDUCE(fcount, fcount_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  !call MPI_ALLREDUCE(vel,vel_all,vel_size*8, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  !print*, "reached 2"
  
  cnt = 0
  do i=1,nrmax
      do j=1,nzmax
          cnt = cnt+1
          r2 = i*dr
          r1 = (i-1)*dr
          z = dr
          vol = pi*(r2**2-r1**2)*z
          n = vol/density
          den = vel(cnt,3)/fcount/n;
          vel(cnt,8) = den
          !if (myid.eq.0) den = vel_all(cnt,3)/fcount_all/n
          !if (myid.eq.0) vel_all(cnt,8) = den
          if(vel(cnt,3) .ge. 0.9) then  ! just to make sure that it has some atoms
              vel(cnt,4) = vel(cnt,4)/vel(cnt,3)*vel(cnt,8)
              vel(cnt,5) = vel(cnt,5)/vel(cnt,3)*vel(cnt,8)
              vel(cnt,6) = vel(cnt,6)/vel(cnt,3)*vel(cnt,8)
              vel(cnt,7) = vel(cnt,7)/vel(cnt,3)*vel(cnt,8)
              !if (myid.eq.0) then
              !    vel_all(cnt,4) = vel_all(cnt,4)/vel_all(cnt,3)*vel_all(cnt,8)
              !    vel_all(cnt,5) = vel_all(cnt,5)/vel_all(cnt,3)*vel_all(cnt,8)
              !    vel_all(cnt,6) = vel_all(cnt,6)/vel_all(cnt,3)*vel_all(cnt,8)
              !    vel_all(cnt,7) = vel_all(cnt,7)/vel_all(cnt,3)*vel_all(cnt,8)
              !end if
          end if
      end do
  end do

  call MPI_REDUCE(vel,vel_all,vel_size*8,MPI_REAL,MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if(myid .eq. 0) then
      vel_all = vel_all/real(Np)
  end if
  
  write(fn,"('velp_',I9.9,'.dat')")fi
  open(29, file=fn)
  if(myid .eq. 0)  open(30, file='vel_all.dat')
  do i=1,cnt
      write(29,"(8E20.10)") vel(i,:)
      if (myid .eq. 0) write(30,"(8E20.10)") vel(i,1:2), vel_all(i,3:8)
  end do
  close(29)
  if(myid .eq. 0) close(30)

  call MPI_FINALIZE(ierr)
  
end program upsidedown
