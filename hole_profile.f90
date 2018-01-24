program hole_profile

  integer :: i,j,k,nz, nt,natm, z, thetai,lt,cnt,slt,maxlt,dlt
  real :: uc, box, theta,dz, pi
  real, dimension(3) :: lb, ctr, xl
  integer, dimension(3) :: nc
  real, allocatable, dimension (:,:) :: hole, X
  real, allocatable, dimension (:) :: avg_hole
  integer, allocatable, dimension (:,:) :: holei
  character(30) :: fn

  pi = 4.0*atan(1.0)
  natm = 401000
  nc = (/ 50, 50, 20 /)
  uc = 5.431073E-10
  lb = nc*uc
  nz = 20
  nt = 24
  dtheta = 360.0/real(nt)
  box = 7.0E-9
  slt = 4000; maxlt = 4000000; dlt = 4000;

  dz = lb(3)/real(nz)
  ctr = lb/2.0

  allocate(hole(nz,nt),holei(nz,nt),X(natm,3),avg_hole(nz))
  
  do lt=slt,maxlt,dlt

      hole = 1.0
      holei = 0
      
      write(fn,"('D/atoms00',I7.7,'.dat')")lt
      print*,fn
      open(29, file=fn)
      
      do i=1,natm
          read(29,*)xl(:)
          !print*,i,xl
          !X(i,:) = xl
          xl(1:2) = xl(1:2)-ctr(1:2)
          if(xl(1).gt.(-1.0)*box .and. xl(1).lt.box .and. xl(2).gt.(-1.0)*box .and. xl(2).lt.box .and. xl(3).ge.0 .and. xl(3).le.lb(3)) then
              z = ceiling(xl(3)/dz)
              if(z.eq.0) z = 1
              theta = atan(xl(2)/xl(1))
              if(xl(1).lt.0) theta = theta+pi
              if(xl(2).lt.0 .and. xl(1).ge.0) theta = theta+2.0*pi
              theta = theta*180.0/pi
              thetai = ceiling(theta/dtheta)
              if(thetai .eq. 0) thetai = 1
              rxy = sqrt(xl(1)*xl(1)+xl(2)*xl(2))
              if(rxy .lt. hole(z,thetai)) then
                  !print*,i,z,theta
                  holei(z,thetai) = i
                  hole(z,thetai) = rxy
              end if
          end if
      end do
      close(29)
      
      avg_hole = 0.0
      write(fn,"('profile/hole_',I7.7,'.dat')")lt
      open(30,file=fn)
      write(fn,"('profile/avg_hole_',I7.7,'.dat')")lt
      open(31, file=fn)
      do i=1,nz
          cnt = 0
          do j=1,nt
              if(holei(i,j) .gt. 0) then
                  cnt = cnt+1
                  write(30,"(3E20.10)")hole(i,j)*cosd(j*dtheta), hole(i,j)*sind(j*dtheta), i*dz
                  avg_hole(i)=avg_hole(i)+hole(i,j)
              end if
          end do
          avg_hole(i) = avg_hole(i)/cnt
          write(31,"(2E20.10)")avg_hole(i),i*dz
      end do
      close(30)
      close(31)
      print*,maxlt/lt*100.0

  end do
      
end program hole_profile
  
  
               
