program radial

  integer :: i, j, k, cnt, nimpact, natm, ntop, nbot, reason, listlt, res_unit, Nt0, ions, zd,nz,n, nimpact_old, rd
  real :: lb3, lb1, lb2,z, top, bot, buff, restart_time,r , rmin,pi, rmin2, m, kb, lb_min, lb_max, dr, vol, topbotclip, ri, zi
  real, allocatable, dimension(:,:) :: X, topx, botx, crysx, V, hole_den, hole_temp, ntime, radial_density, radial_temp, Xi, rho
  real, allocatable, dimension(:)  :: mass
  integer, allocatable, dimension (:) :: atype, donelist
  integer, allocatable, dimension(:,:):: num
  character(100) :: cmd, fn, cn, rn
  real :: rho1, rho2, rho11, rho12, rho21, rho22, rhomup, rhombot

  natm = 5121000
  nimpact_old = 2000
  lb3 = 5.431073E-10*100
  lb2 = 5.431073E-10*80
  lb1 = 5.431073E-10*80
  topbotclip = 0.0e-9
  lb_min = topbotclip
  lb_max = lb3 - topbotclip
  dr = 1.0e-9
  nr = floor(lb1/2.0/dr)
  ntop = 0
  nbot = 0
  top = 4.5
  bot = 6.5
  buff = 0.7
  rmin = 1.0e-9
  rmin2 = 5.0e-9
  pi = 3.141592654
  m = 46.637063E-27
  kb = 1.381E-23
  nz = 10


  allocate(X(natm,3), atype(natm), crysx(natm,3), V(natm,3), mass(natm), Xi(natm,3))
  allocate(topx(natm,8), botx(natm,8), donelist(natm))
  !allocate(hole_den(nimpact,nz+1), num(nimpact,nz+1), hole_temp(nimpact,3), ntime(nimpact,1))
  !allocate(radial_density(nimpact, nr+1), radial_temp(nimpact), num(nimpact,nr+1), ntime(nimpact,1))

  donelist = 0
  num = 0

  ! finding nimpact
  write(cmd,"('ls D/restart.in0* > restart_list.dat')")
  call system(cmd)
  open(33, file='restart_list.dat')
  do k=1,nimpact_old
      read(33,"(A28)",IOSTAT=reason)cn
      if(reason.lt.0) exit
      read (cn(13:18),*)nimpact
  end do
  close(33)
  nimpact = nimpact-1


  !open(31, file="G/atoms000000000.dat")
  !do i=1,natm
  !    read(31,*) Xi(i,:)
  !end do
  !close(31)


  allocate(radial_density(nimpact, nr+1), radial_temp(nimpact, nr+1), num(nimpact,nr+1), ntime(nimpact,1), rho(nimpact,4))
  radial_density = 0.0
  radial_temp = 0.0
  !print*, nr
  !call exit(0)

  !open(41, file='rhomix.dat')

  do k=1,nimpact
      write(cmd,"('ls D/restart.in',I6.6, '* > restart_list.dat')")k
      call system (cmd)
      !print*,cmd
      open(33, file='restart_list.dat')
      read(33,"(A28)",IOSTAT=reason)fn
      close (33)
      print*, trim(fn)

      open(res_unit,file=trim(fn),form='UNFORMATTED',status='OLD')
      read(res_unit)Nt0,restart_time
      read(res_unit)X,V,mass,atype
      read(res_unit)ions

      write(rn,"('XZfiles/XZ_',I4.4,'.dat')")k
      open (22,file=rn)
      do j=1,natm
          if(X(j,1).lt.(lb1/2.0+0.5e-9) .and. X(j,1).gt.(lb1/2.0-0.5e-9)) then
              if(X(j,3).gt.4.0e-7) X(j,3) = X(j,3) - lb3*10
              write(22,"(3E12.4, I2, E11.4)")X(j,:),atype(j),m*sum(V(j,:)**2)/3.0/kb
          end if
      end do
      ! needed this line below for the position of ions
      !write(22,"(3E12.4, I2, E11.4)")X(natm-1000+k,:),atype(natm-1000+k),m*sum(V(natm-1000+k,:)**2)/3.0/kb
      close(22)
  end do


end program radial
