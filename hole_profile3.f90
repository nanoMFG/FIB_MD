program hole_profile

  integer :: i,j,k,nz, nt,natm, z, thetai,lt,cnt,slt,maxlt,dlt, ndr, nrxy
  real :: uc, box, theta,dz, pi, dr,r, density_ref
  real, dimension(3) :: lb, ctr, xl
  integer, dimension(3) :: nc
  real, allocatable, dimension (:,:) :: hole, X, density, vol
  real, allocatable, dimension (:) :: avg_hole
  integer, allocatable, dimension (:,:) :: holei,nrad
  character(30) :: fn

  pi = 4.0*atan(1.0)
  natm = 401000
  nc = (/ 50, 50, 20 /)
  uc = 5.431073E-10
  lb = nc*uc
  nz = 20
  nt = 24
  dr = 1.0E-10
  density_ref = uc*uc*uc/8.0
  ndr = int(min(lb(1)/2.0,lb(2)/2.0)/dr)
  dtheta = 360.0/real(nt)
  box = 7.0E-9
  slt = 0; maxlt = 1000000; dlt = 4000;

  dz = lb(3)/real(nz)
  ctr = lb/2.0

  allocate(hole(nz,nt),holei(nz,nt),X(natm,3),avg_hole(nz), nrad(nz,ndr),vol(nz,ndr),density(nz,ndr))
  do i=1,nz
      r = 0.0
      do j = 1,ndr
          r = j*dr
          vol(i,j) = pi*( r*r-(r-dr)*(r-dr) ) * dz
          print*,vol(i,j)
      end do
  end do
  write(33,*) density_ref
  
  do lt=slt,maxlt,dlt

      hole = 1.0
      holei = 0
      
      write(fn,"('D/atoms00',I7.7,'.dat')")lt
      print*,fn
      open(29, file=fn)
      nrad = 0
      do i=1,natm
          read(29,*)xl(:)
          xl(1:2) = xl(1:2)-ctr(1:2)
          rxy = sqrt(xl(1)*xl(1)+xl(2)*xl(2))
          nrxy = ceiling(rxy/dr)
          z = ceiling(xl(3)/dz); if(z.eq.0) z = 1
          nrad(z,nrxy) = nrad(z,nrxy)+1
      end do

      write(fn,"('profile/rad_hole_',I7.7,'.dat')")lt
      open(31, file=fn)

      write(32,*)nrad
      density = 0.0
      do i=1,nz
          do j=1,ndr
              if(nrad(i,j).ne.0) density(i,j) = vol(i,j)/nrad(i,j)
              !if(density(i,j) .gt. 1.2*density_ref) then
              !    write(31,"(3E20.10)")j*dr,i*dz,density(i,j)
              !    exit
              !end if
          end do
       
      end do
      write(34,*)density
      close(31)
      print*,lt/maxlt*100.0
      print*, ndr
      stop
  end do
      
end program hole_profile
  
  
               
