program midsection

  integer :: i, j, k, cnt, nimpact, natm, ntop, nbot, reason, listlt, res_unit, Nt0, ions, zd
  real :: lb3, lb1, lb2,z, top, bot, buff, restart_time,r , rmin,pi,m,kb
  real, allocatable, dimension(:,:) :: X, topx, botx, crysx, V, hole_den, outx
  real, allocatable, dimension(:)  :: mass
  integer, allocatable, dimension (:) :: atype, donelist
  integer, allocatable, dimension(:,:):: num
  character(100) :: cmd, fn

  natm = 401000
  nimpact = 700
  lb3 = 5.431073E-10*20
  lb2 = 5.431073E-10*50
  lb1 = 5.431073E-10*50
  ntop = 0
  nbot = 0
  top = 4.5
  bot = 6.5
  buff = 0.5E-9
  rmin = 0.5e-9
  pi = 3.141592654
  m = 46.637063E-27
  kb = 1.381E-23


  allocate(X(natm,3), atype(natm), crysx(natm,3), V(natm,3), mass(natm),outx(natm,7))
  allocate(topx(natm,8), botx(natm,8), donelist(natm))
  allocate(hole_den(nimpact,11), num(nimpact,11))

  donelist = 0

  !open(31, file="D/atoms000000000.dat")
  !do i=1,natm
  !    read(31,*) crysx(i,:)
  !end do
  !close(31)

  hole_den = 0.0
  num = 0

  do k=240,nimpact
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
      close (res_unit)
      
      cnt = 0
      do i=1,natm
          if(X(i,2).lt.(lb2/2.0+buff) .and. X(i,2).gt.(lb2/2.0-buff)) then
              if(X(i,3).gt.8.0E-8 .or. X(i,3).lt.2.0E-8) then
                  if(X(i,3).gt.8.0E-8) X(i,3) = X(i,3)-10.0*lb3
                  cnt = cnt+1
                  outx(cnt,1:3) = X(i,:)
                  outx(cnt,4:6) = V(i,:)
                  outx(cnt,7) = m*(sum(V(i,:)**2))/3.0/kb
              end if
          end if
      end do

      write(fn,"('out2/mid_',I4.4,'.dat')")k
      open(21, file=fn)

      do i=1,cnt
          write(21,"(7E20.10)")outx(i,:)
      end do
      close(21)
  end do

end program midsection
