program test

  integer :: i,j
  real :: x
  integer(2) :: p
  real, allocatable, dimension(:,:) :: pp
  !real, dimension(2000000000) :: pp
  integer, dimension(8) :: seed1,seed2

  x = 240.3123445
  write(*,"(F5.2)")x
  

  print*,'             |'
  print*,'           .\:/.'
  print*,'KNOCK    --=>*<=--'
  print*,"           //:\\"

  do i=1,0
      print*, i
  end do
  
  !allocate(pp(2000000000))
  !call DATE_AND_TIME(values=seed1)
  !print*, seed1(8)
  !do j=1,100
  !    do i=1,20000000
  !        pp(i)=123.0
  !    end do
  !end do
  !call DATE_AND_TIME(values=seed2)
  !print*, seed2(8)
  !print*, 'time in do loop', seed2(7)-seed1(7), 'seconds',seed2(8)-seed1(8)

  !call DATE_AND_TIME(values=seed1)
  !do j=1,100
  !    pp(:)=123.0
  !end do
  !call DATE_AND_TIME(values=seed2)
  !print*, 'time in direct indexing', seed2(7)-seed1(7), 'seconds',seed2(8)-seed1(8)

  !print*,i
  
  !x = 4.5E-10
  !do i = 1,100
  !    
  !    if(MOD(i,5).eq.0) cycle
  !    print*,i
  !    p = i
  !    pp(i) = i
      
  !end do

  !print*, MAXVAL(pp), MAXVAL(pp(:))

  !write(*,"(E10.2E3)")x

  !print*,MIN(5,2,3,1,4)

  
end program test
