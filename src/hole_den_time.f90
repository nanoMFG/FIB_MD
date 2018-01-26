program hole_den_time

  integer :: i, j, k, cnt, nimpact, natm, ntop, nbot, reason, listlt, res_unit, Nt0, ions, zd
  real :: lb3, lb1, lb2,z, top, bot, buff, restart_time,r , rmin,pi, rmin2, m, kb
  real, allocatable, dimension(:,:) :: X, topx, botx, crysx, V, hole_den, hole_temp, ntime
  real, allocatable, dimension(:)  :: mass
  integer, allocatable, dimension (:) :: atype, donelist
  integer, allocatable, dimension(:,:):: num
  character(100) :: cmd, fn

  natm = 5121000
  nimpact = 54
  lb3 = 5.431073E-10*100
  lb2 = 5.431073E-10*80
  lb1 = 5.431073E-10*80
  ntop = 0
  nbot = 0
  top = 40.0E-9
  bot = 30.0E-9
  buff = 0.7
  rmin = 0.5e-9
  rmin2 = 5.0e-9  !12.0E-9
  pi = 3.141592654
  m = 46.637063E-27
  kb = 1.381E-23


  allocate(X(natm,3), atype(natm), crysx(natm,3), V(natm,3), mass(natm))
  allocate(topx(natm,8), botx(natm,8), donelist(natm))
  allocate(hole_den(nimpact,11), num(nimpact,11), hole_temp(nimpact,3), ntime(nimpact,1))

  donelist = 0

  !open(31, file="D/atoms000000000.dat")
  !do i=1,natm
  !    read(31,*) crysx(i,:)
  !end do
  !close(31)

  hole_den = 0.0
  num = 0
  hole_temp = 0.0

  do k=1,nimpact
      !ntop = 0; nbot = 0;
      write(cmd,"('ls D/restart.in',I6.6, '* > restart_list.dat')")k
      !print*, cmd
      call system (cmd)
      open(33, file='restart_list.dat')
      read(33,"(A28)",IOSTAT=reason)fn
      close (33)
      !read(fn(20:28),"(I9)")listlt
      !print*, listlt

      !write(fn,"('D/atoms',I9.9,'.dat')")listlt
      print*, trim(fn)
      open(res_unit,file=trim(fn),form='UNFORMATTED',status='OLD')
      !open(29,file=fn)
      read(res_unit)Nt0,restart_time
      read(res_unit)X,V,mass,atype
      read(res_unit)ions
      X(:,1) = X(:,1)-lb1/2.0
      X(:,2) = X(:,2)-lb2/2.0
      ntime(k,1) = restart_time
      hole_den(k,1) = k
      num(k,1) = k
      hole_temp(k,1) = k
      
      do i=1,natm
          !print*,X(i,3)
          if(X(i,3) .lt. top .and. X(i,3) .gt. bot) then
              r = sqrt(X(i,1)**2+X(i,2)**2)
              !print*, r
              if(r.lt.rmin) then
                  zd = ceiling(X(i,3)/lb3*10.0)
                  num(k,zd+1) = num(k,zd+1)+1
                  
              end if
              if (r.lt.rmin2) then
                  hole_temp(k,2) = hole_temp(k,2) + 1
                  hole_temp(k,3) = hole_temp(k,3) + m*(sum(V(i,:)**2))/3.0/kb
              end if
          end if
      end do
      !print*,pi*rmin**2*lb3/(5.431073E-10**3)*8.0/10.0
      hole_den(k,2:11) = num(k,2:11)/(pi*rmin**2*lb3/(5.431073E-10**3)*8.0/10.0)
      hole_temp(k,3) = hole_temp(k,3)/hole_temp(k,2)
      !hole_den(k,3) = hole_den(k,2)/(pi*rmin**2*lb3/(5.431073E-10**3)*8.0)
      close(res_unit)
  end do

  open(21, file='impact_vs_density.dat')
  open(22, file='impact_vs_temp.dat')

  do i=1,nimpact
      write(21,"(I6, 11E20.10)")num(i,1),ntime(i,1),hole_den(i,2:11)
      write(22,"(I6, 3E20.10)") num(i,1),ntime(i,1),hole_temp(i,2:3)
  end do

  close(21)
  close(22)
  
end program hole_den_time
