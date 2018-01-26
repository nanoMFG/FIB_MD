program upsidedown

  integer :: i, j, k, natm, fi, fo, df, nzmax, nrmax, cnt, nr, nz, dum, ff, ii, io, reason, fcount
  real, allocatable, dimension(:,:) :: X, V, F, vel
  real, allocatable, dimension(:) :: Lb, kin_eng
  integer, allocatable, dimension(:) :: atype, p
  real :: dr, r, dz, vr, vt, ev, kb, r2, r1, z, vol, n, density, den, pi
  character (100) :: cmd
  character (30)  :: fn

  natm = 401000
  allocate(X(natm,3), V(natm,3), F(natm,3), kin_eng(natm), Lb(3))
  allocate(atype(natm),p(natm))
  Lb = (/ 2.7155365E-8, 2.7155365E-8, 1.0862146E-8 /)
  kb = 1.3806503E-23
  ev = 1.60217646E-19
  density = (5.431073E-10**3)/8.0
  pi = atan(1.0)*4.0
  
  ii = 299
  io = ii+10

  write(cmd,"('ls D/restart.in',I6.6, '* > restart_list.dat')")ii
  call system (cmd)
  open(33, file='restart_list.dat')
  read(33,"(A28)",IOSTAT=reason)fn
  close (33)
  read(fn(20:28),"(I9)")fi

  write(cmd,"('ls D/restart.in',I6.6, '* > restart_list.dat')")io
  call system (cmd)
  open(33, file='restart_list.dat')
  read(33,"(A28)",IOSTAT=reason)fn
  close (33)
  read(fn(20:28),"(I9)")fo
  
  df = 4000
  fi = NINT(real(fi)/real(df))*df
  fo = NINT(real(fo)/real(df))*df
  ! putting the values of fi fo below to calculate the velocity of flow after bombardment
  !fi = 4412000
  !fo = 5868000
  !ii = 999
  
  print*, fi, fo
  
  dr = 0.5E-9
  dz = 5.0E-9
  nzmax = NINT((2.0E-8+dz)/dr) + 1
  nrmax = NINT(Lb(1)/2.0/dr) + 1
  cnt = nrmax*nzmax

  print*, cnt
  allocate(vel(cnt,8))
  vel = 0.0

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
      fcount = fcount + 1
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

              if(nr.lt.nrmax .and. nz.lt.nzmax) then
                  vel(dum,3) = vel(dum,3)+1
                  vel(dum,4) = vel(dum,4)+vr
                  vel(dum,5) = vel(dum,5)+vt
                  vel(dum,6) = vel(dum,6)+V(i,3)
                  vel(dum,7) = vel(dum,7)+kin_eng(i)*ev/kb*2.0/3.0
              end if

          end if
      end do
      close(29)
  end do

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
          if(vel(cnt,3) .ne. 0) then
              vel(cnt,4) = vel(cnt,4)/vel(cnt,3)*vel(cnt,8)
              vel(cnt,5) = vel(cnt,5)/vel(cnt,3)*vel(cnt,8)
              vel(cnt,6) = vel(cnt,6)/vel(cnt,3)*vel(cnt,8)
              vel(cnt,7) = vel(cnt,7)/vel(cnt,3)*vel(cnt,8)
          end if
      end do
  end do

  write(fn,"('velf_',I3.3,'.dat')")ii
  open(29, file=fn)
  do i=1,cnt
      write(29,"(8E20.10)") vel(i,:)
  end do
  close(29)
  
end program upsidedown
