program write_sput_files

  integer :: i,j,k
  character(30) :: fn, fn2
  real,allocatable, dimension(:,:) :: indx
  integer, allocatable, dimension(:,:) :: index
  real, allocatable,dimension(:,:) :: X
  real :: lb3

  lb3 = 108.62146E-9

  allocate(X(401000,3))
  allocate(indx(1993,2),index(1993,2))
  open(19, file='sput_index.dat')
  do i=1,1993
      read(19,*)indx(i,:)
  end do

  index = NINT(indx)

  k = 1
  do i=0,2008000,4000
      k = 1
      write(fn,"('D/atoms00',I7.7,'.dat')"),i
      write(fn2,"('sput_files/sput_',I7.7,'.dat')"),i
      open(19,file=fn)
      open(20,file=fn2)
      do j=1,400000
          read(19,*)X(j,:)
          if(index(k,1).eq.j .and. k.le.1993) then
              if (X(j,3).gt. 100.0E-9) X(j,3) = X(j,3)-lb3
              write(20,"(3E20.10,2I7)")X(j,:),index(k,2), j
              k = k+1
          end if
      end do
      close(19)
      close(20)
      print*,i
  end do

  deallocate(X, indx, index)
  
end program write_sput_files

      
