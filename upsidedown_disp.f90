program upsidedown_disp

  integer :: i, j, k, natm, fi, fo, df, nzmax, nrmax, cnt, nr, nz, dum, ff, ii, io, reason
  real, allocatable, dimension(:,:) :: X, V, F, vel, Xi
  real, allocatable, dimension(:) :: Lb, kin_eng, time
  integer, allocatable, dimension(:) :: atype, p, dx
  real :: dr, r, dz, vr, dt
  character (100) :: cmd
  character (30)  :: fn

  natm = 401000
  allocate(X(natm,3), V(natm,3), F(natm,3), kin_eng(natm), Lb(3), Xi(natm,3), dx(2), time(43103))
  allocate(atype(natm),p(natm))
  Lb = (/ 2.7155365E-8, 2.7155365E-8, 1.0862146E-8 /)

  ii = 300
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
  
  dr = 1.0E-9
  dz = 5.0E-9
  nzmax = NINT((2.0E-8+dz)/dr) + 1
  nrmax = NINT(Lb(1)/2.0/dr) + 1
  cnt = nrmax*nzmax

  print*, cnt
  allocate(vel(cnt,5))
  vel = 0.0

  cnt = 0
  do i=1,nrmax
      do j=1,nzmax
          cnt = cnt+1
          vel(cnt,1) = i*dr
          vel(cnt,2) = j*dr
      end do
  end do

  write(fn,"('D/atoms',I9.9,'.dat')") fi-df
  open(29, file=fn)
  do i=1,natm
      read(29,*) Xi(i,:),atype(i),P(i),kin_eng(i), V(i,:)
      if(Xi(i,3).gt. 8.0E-8) Xi(i,3)= Xi(i,3)-5.431073E-9*20
      Xi(i,1) = Xi(i,1) - Lb(1)/2.0
      Xi(i,2) = Xi(i,2) - Lb(2)/2.0
      Xi(i,3) = Xi(i,3) + dz
      
  end do
  close(29)

  open(29, file= 'temperature.dat')
  do i=1,43103
      read(29,*) time(i)
  end do
  close(29)
      
  do ff=fi,fo,df
      write(fn,"('D/atoms',I9.9,'.dat')")ff
      dt = time(NINT(ff/100.0)) - time(NINT((ff-df)/100.0))
      print*, fn, dt

      open(29, file=fn)
      do i=1,natm
          read(29,*) X(i,:),atype(i),P(i),kin_eng(i), V(i,:)
          if(X(i,3).gt. 8.0E-8) X(i,3)= X(i,3)-5.431073E-9*20
          X(i,1) = X(i,1) - Lb(1)/2.0
          X(i,2) = X(i,2) - Lb(2)/2.0
          X(i,3) = X(i,3) + dz
          
          if(X(i,3).lt. 2.0E-8+dz .and. atype(i).lt.2) then
              r = sqrt(X(i,1)**2+X(i,2)**2)
              nr = FLOOR(r/dr)
              nz = FLOOR(X(i,3)/dr)
              dum = nz + nr*nzmax + 1
              dx = (X(i,1:2)-Xi(i,1:2))/dt
              vr = sum(dx(1:2)*X(i,1:2))/r

              if(nr.lt.nrmax .and. nz.lt.nzmax) then
                  vel(dum,3) = vel(dum,3)+1
                  vel(dum,4) = vel(dum,4)+vr
                  vel(dum,5) = vel(dum,5)+ (X(i,3)-Xi(i,3))/dt
              end if

          end if
      end do
      close(29)
  end do

  cnt = 0
  do i=1,nrmax
      do j=1,nzmax
          cnt = cnt+1
          if(vel(cnt,3) .ne. 0) then
              vel(cnt,4) = vel(cnt,4)/vel(cnt,3);
              vel(cnt,5) = vel(cnt,5)/vel(cnt,3);
          end if
      end do
  end do

  write(fn,"('veld_',I3.3,'.dat')")ii
  open(29, file=fn)
  do i=1,cnt
      write(29,"(5E20.10)") vel(i,:)
  end do
  close(29)
  
end program upsidedown_disp
