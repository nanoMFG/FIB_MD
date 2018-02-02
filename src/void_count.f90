program void_count

  integer :: i, j, k, m, n, reason, natm, atm, nx, nz, cx, cz, dum, listcnt, tempcnt
  real :: lb1, lb2, lb3, lb, temp
  real, allocatable, dimension(:,:) :: X, hole_ini, list, hole_list, hole, V
  real, allocatable, dimension(:) :: temperature
  character(100) :: fn

  natm = 400000
  n = 12958
  lb = 5.431073E-10
  nx = 50
  nz = 20
  listcnt = 0
  temp = 0.0
  V = 0.0
  temperature = 0.0
  tempcnt = 0
  
  allocate(X(natm,3),hole(nx*nz,4), hole_ini(nx*nz,4), list(nx*nz,4),hole_list(n,3), V(natm,3), temperature(natm))
  hole = 0.0
  do cx=1,nx
      do cz=1,nz
          dum = (cx-1)*nz + cz
          hole_ini(dum,1) = cx*lb-0.5*lb
          hole_ini(dum,2) = cz*lb-0.5*lb
          !print*, hole(dum,1:2)
      end do
  end do
  print*, dum

  open(22, file='hole_volume.dat')
  
  do j=1,n
      hole = 0.0
      write(fn,"('out/mid_',I9.9,'.dat')")j*1000
      open(21, file=fn)
      temp = 0.0
      temperature = 0.0
      tempcnt = 0
      do k=1,natm
          read(21,*,IOSTAT=reason) X(k,:), V(k,:), temperature(k)
          cx = ceiling(abs(X(k,1))/lb)
          cz = ceiling(abs(X(k,3))/lb)
          if(cx .gt. nx) cx = nx
          if(cz .gt. nz) cz = nz
          dum = (cx-1)*nz + cz
          hole(dum,3) = hole(dum,3)+1.0

          if(cx .gt. 20 .and. cx .lt. 30) then
              tempcnt = tempcnt+1
              temp = temp + temperature(k)
          end if
          
          if(dum .gt. nx*nz) print*, 'error!!'
          if(reason.lt.0) then
              atm = k-1
              !print*,'Found', atm, 'atoms'
              exit
          end if
      end do
      close(21)

      temp = temp/tempcnt
      
      listcnt = 0
      do cx=1,nx
          do cz=2,nz-1
              dum = (cx-1)*nz + cz
              hole(dum,4) = hole(dum,3)/14.0
              if(hole(dum,4).lt. 0.5) then
                  listcnt = listcnt+1
                  list(listcnt,1) = cx
                  list(listcnt,2) = cz
              end if
          end do
      end do

      hole_list(j,1) = j*1000
      hole_list(j,2) = listcnt
      hole_list(j,3) = temp
      
                  

      write(22, "(3E20.10)") hole_list(j,:)
              

      print*,j
      
  end do

  
  
  
  close(22)

  open(23, file='hole_test.dat')
  do i=1,nx*nz
      write(23,"(4E20.10)")hole_ini(i,1:2),hole(i,3:4)
  end do
  
end program void_count
