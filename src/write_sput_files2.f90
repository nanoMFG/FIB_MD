program write_sput_files

  integer :: i,j,k,m
  character(30) :: fn, fn2
  real,dimension(1993,2) :: indx
  integer,dimension(1993,2) :: index
  real,dimension(1993,503) :: x,y,z
  real :: lb3
  real, dimension(3) :: dum

  lb3 = 108.62146E-9

  open(19, file='sput_index.dat')
  do i=1,1993
      read(19,*)indx(i,:)
  end do

  index = NINT(indx)
  m=0
  k=1

  do i=0,2008000,4000
      m=m+1
      k=1
      write(fn,"('D/atoms00',I7.7,'.dat')"),i
      open(19,file=fn)
      do j=1,400000
          read(19,*)dum
          if(index(k,1).eq.j .and. k.le.1993) then
              if (dum(3).gt. 100.0E-9) dum(3) = dum(3)-lb3
              x(k,m) = dum(1)
              y(k,m) = dum(2)
              z(k,m) = dum(3)
              k=k+1
          end if
      end do
      close(19)
      print*,i
  end do

  do k=1,1993
      write(fn,"('sput_files/traj_',I7.7,'.dat')"),k
      open(19,file=fn)
      do m=1,503
          write(19,"(3E20.10,I7)"),x(k,m),y(k,m),z(k,m),index(k,2)
      end do
      close(19)
      print*,k
  end do
    
end program write_sput_files

      
