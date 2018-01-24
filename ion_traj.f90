program ion_traj

  integer :: i,j,k, lt, maxlt, dlt,natm,cnt,cnt2
  real :: lb3, lb3p
  character(30) :: fn
  integer, allocatable,dimension(:) :: list
  real, allocatable, dimension(:,:) :: X, xl,yl,zl

  maxlt = 4000000
  dlt = 20000
  natm = 401000
  lb3 = 5.431073E-10*20*10
  lb3p = lb3/2.0
  nlist = 429
  allocate(X(natm,3), xl(nlist,maxlt/dlt+1), yl(nlist,maxlt/dlt+1),zl(nlist,maxlt/dlt+1), list(nlist))
  

  open (29, file = 'atom_traj_indx.dat')
  do i=1,nlist
      read(29,*) list(i)
  end do
  close(29)
  print*, list(:)

  cnt = 1
  j = list(cnt)
  cnt2 = 0
  do lt = 0,maxlt,dlt
      cnt2 = cnt2+1
      write(fn,"('D/atoms',I9.9,'.dat')")lt
      open(29, file=fn)
      cnt = 1
      j = list(cnt)
      do i=1,natm
          if(i.eq.j) then
              
              read(29,*),xl(cnt,cnt2),yl(cnt,cnt2),zl(cnt,cnt2)
              if(zl(cnt,cnt2) .gt. lb3p) zl(cnt,cnt2) = zl(cnt,cnt2)-lb3
              cnt = cnt+1
              if(cnt .le. nlist)then
                  j = list(cnt)
              else
                  exit
              end if
          else
              read(29,*)X(i,:)
          end if
      end do
      close(29)
      print*,lt
  end do

  do i=1,nlist
      write(fn,"('ion_traj/traj_', I6.6,'.dat')")i
      open(29, file=fn)
      do j=1,cnt2
          write(29,"(3E20.10)")xl(i,j),yl(i,j),zl(i,j)
      end do
      close(29)
      print*, i
  end do

end program ion_traj
      
      
  
