program count_sputter

  integer :: i, j, k, cnt, nimpact, natm, ntop, nbot, reason, listlt, res_unit, Nt0, ions
  real :: lb,z, top, bot, buff, restart_time
  real, allocatable, dimension(:,:) :: X, topx, botx, crysx, V
  real, allocatable, dimension(:)  :: mass
  integer, allocatable, dimension (:) :: atype, donelist
  character(100) :: cmd, fn
  
  natm = 401000
  nimpact = 680
  lb = 5.431073E-10*20
  ntop = 0
  nbot = 0
  top = 4.5
  bot = 6.5
  buff = 0.7

  
  allocate(X(natm,3), atype(natm), crysx(natm,3), V(natm,3), mass(natm))
  allocate(topx(natm,8), botx(natm,8), donelist(natm))

  donelist = 0
  
  open(31, file="D/atoms000000000.dat")
  do i=1,natm
      read(31,*) crysx(i,:)
  end do
  close(31)
  
  
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
      
      
      do i=1,natm
          !read(29,*)X(i,:), atype(i)
          
          if(X(i,3) .lt. (top+buff)*lb .and. X(i,3) .gt. (lb+2.0E-9) .and. donelist(i).eq.0) then
              if(atype(i) .lt. 3) then
                  ntop = ntop + 1
                  topx(ntop,1:3) = X(i,1:3)
                  topx(ntop,4:6) = crysx(i,1:3)
                  topx(ntop,7) = real(ntop)
                  topx(ntop,8) = real(k)
                  donelist(i) = 1
              end if
          elseif (X(i,3) .gt. (bot-buff)*lb .and. X(i,3) .lt. (lb*10.0-2.0E-9) .and. donelist(i).eq.0) then
              if(atype(i) .lt. 3) then
                  nbot = nbot + 1
                  botx(nbot,1:3) = X(i,1:3)
                  botx(nbot,4:6) = crysx(i,1:3)
                  botx(nbot,7) = real(nbot)
                  botx(nbot,8) = real(k)
                  print*, real(k)
                  donelist(i) = 1
              end if
          end if
      end do
      !close (29)
      close(res_unit)
      print*, k, ntop, nbot
      open(91, file='topbotsputter.dat', POSITION='APPEND')
      write(91,"(3I9)") k, ntop, nbot
      close(91)
  end do

  open (91, file='topx.dat')
  do i=1,ntop
      write(91,"(8E20.10)") topx(i,:)
  end do
  close(91)

  open (91, file='botx.dat')
  do i=1,nbot
      write(91,"(8E20.10)") botx(i,:)
  end do
  close(91)
      

end program count_sputter
