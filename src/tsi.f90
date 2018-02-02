PROGRAM atmtopdb

  integer, parameter      :: Nm=90000
  real(8), dimension(Nm,3) :: X
  real(8)                 :: Xm1,Xm2,Xm3
  integer, dimension(Nm)  :: atype
  integer                 :: i,it
  character*3 ctyp(4)

  ctyp(1) = 'O'
  ctyp(2) = 'N'
  ctyp(3) = 'N'
  ctyp(4) = 'N'

  Xm1 = 0.; Xm2 = 0.; Xm3 = 0.
  do i = 1,Nm
     read(5,*,END=100)X(i,1),X(i,2),X(i,3),atype(i)
  end do
100 continue
  it = i -1

  X = X*1.E10/2.

  do i = 1,it
     write(6,101) i,ctyp(atype(i)),i,X(i,1), &
             X(i,2),X(i,3)
  end do
  
101 FORMAT('ATOM  ',I5,1X,' ',A,'     ',1X,I4,2X,3F8.3)

END PROGRAM atmtopdb



