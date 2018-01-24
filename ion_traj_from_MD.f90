program iontraj

  integer :: i, j, k, nfile,reason, tcount
  real :: m, ec
  character(100) :: fn, cmd
  real, dimension(100000,12) :: b

  m = 1.157777417E-25
  ec = 1.602e-19
  nfile = 4
  
  do i=1,nfile
      write(cmd,"('ls iontraj/ion_traj_', I4.4, '_* > ionlist.dat')")i
      print*, cmd
      call system(cmd)
      open(33, file='ionlist.dat')
      tcount = 1
      
      do j=1,10
          read(33, "(A100)", IOSTAT=reason)fn
          if(reason .lt. 0) then
              exit
          end if
          print*,fn

          open(35, file=fn)
          do k=1,100000
              read(35,*,IOSTAT=reason)b(tcount,:)
              if(reason .lt. 0) then
                  
                  exit
              end if
              tcount = tcount+1
          end do
          print*, b(tcount-1,12), tcount-1
          
          
          close(35)
          
      end do

      
      
      close(33)

  end do
  
  
end program iontraj
