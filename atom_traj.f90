program atom_traj

  integer :: i, j, k, cnt, nimpact, natm, ntop, nbot, reason, listlt, res_unit, Nt0, ions, zd, na
  real :: lb3, lb1, lb2,z, top, bot, buff, restart_time,r , rmin,pi
  real, allocatable, dimension(:,:) :: X, topx, botx, crysx, V, hole_den, xlist, vlist
  real, allocatable, dimension(:)  :: mass
  integer, allocatable, dimension (:) :: atype, donelist, list
  integer, allocatable, dimension(:,:):: num
  character(100) :: cmd, fn
  
  natm = 401000
  nimpact = 755
  lb3 = 5.431073E-10*20
  lb2 = 5.431073E-10*50
  lb1 = 5.431073E-10*50
  ntop = 0
  nbot = 0
  top = 4.5
  bot = 6.5
  buff = 0.7
  rmin = 0.5e-9
  pi = 3.141592654
  
  
  allocate(X(natm,3), atype(natm), crysx(natm,3), V(natm,3), mass(natm))
  allocate(topx(natm,8), botx(natm,8), donelist(natm), list(natm))
  allocate(hole_den(nimpact,11), num(nimpact,11))
  
  donelist = 0
  
  open(31, file="D/atoms000000000.dat")
  do i=1,natm
      read(31,*) crysx(i,:)
  end do
  close(31)
  
  hole_den = 0.0
  num = 0
  
  open(31, file="atom_traj_list.dat")
  do i=1,natm
      read(31,*,IOSTAT=reason)list(i)
      if(reason .lt. 0) then
          na = i-1
          print*, 'Done list loading, number of atoms = ', na
          exit
      end if
  end do
  
  allocate(xlist(nimpact, na*3), vlist(nimpact, na*3))
  
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
      close(res_unit)
      
      do j=1,na
          i= list(j)
          xlist(k,(j-1)*3+1:(j-1)*3+3) = X(i,:)
          vlist(k,(j-1)*3+1:(j-1)*3+3) = V(i,:)
      end do
  end do
  
  do i=1,na
      write(fn, "('atom_traj/atom_traj_',I8.8,'.dat')")list(i)
      open(31,file=fn)
      do k = 1,nimpact
          write(31,"(6E20.10, 2I9)")xlist(k,(i-1)*3+1:(i-1)*3+3), vlist(k,(i-1)*3+1:(i-1)*3+3), k, list(i)
      end do
      close(31)
  end do
  
end program atom_traj


