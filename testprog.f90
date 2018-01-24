program test

  integer :: i, j, k
  real :: sum
  
  print*, 'Testing...'

  sum = 0.0
  do i=1,10
      sum = sqrt(2.0) + sum
  end do

  print*, 'Sum = ', sum
  
end program test
