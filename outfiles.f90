program outfiles

  integer :: i,j,k, Natm, startf, endf, df, atype, p
  character (30) :: fn
  real :: ke,kb,ec
  real, dimension(3) :: lb
  !real, allocatable,dimension(:,:) :: X, V
  real, dimension(401000,3) :: X,V
  real, dimension(401000) :: T
  !real, allocatable, dimension(:) :: T

  lb = (/ 2.7155365E-8, 2.7155365E-8, 1.0862146E-7  /)
  Natm = 401000
  startf = 4312000
  endf = 5868000
  df = 4000
  kb = 1.3806503E-23
  ec = 1.60217646E-19

  !allocate(X(Natm,3),V(Natm,3), T(Natm))
  
  do k=startf, endf, df
      write(fn,"('D/atoms00', I7.7,'.dat')")k
      open(23, file=fn)
      do i = 1,Natm
          read(23,*)X(i,:), atype, p, ke, V(i,:)
          T(i) = ke*ec/(1.5*kb)
      end do
      close(23)

      
      write(fn,"('outfiles/out_',I7.7,'.dat')")k
      open(23, file=fn)
      do i=1,Natm
          if(X(i,3).lt.4.0E-8 .or. X(i,3).gt.8.0E-8) then
              if(X(i,3).gt.4.0E-8) X(i,3)=X(i,3)-lb(3)
              write(23,"(4E15.6,I8)")X(i,:),T(i), i
          end if
      end do
      print*, "done ",real((k-startf))/real((endf-startf))*100.0,"%"

  end do

end program outfiles
      
          
  
