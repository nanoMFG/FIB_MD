program velocity_avg

  integer :: i,j,k,nx,ny,nz,ux,uy,uz,lt,natm,ncell,n,dum,maxlt,dlt,slt,reason,nimpact,fl,listlt,klt
  real :: dx,dy,dz,uc
  real, dimension(3) :: lb, lbs
  integer, dimension(3) :: ijk,npc,nc
  real, allocatable, dimension(:,:) :: X,V,F,xl,vl
  real, allocatable, dimension(:)   :: kin
  integer, allocatable, dimension(:):: atype,p,id,cnt
  character(100) :: fn, cmd

  natm = 400000; maxlt = 12000; dlt = 4000; slt = 0;
  uc = 5.4310729e-10
  nc(1) = 50; nc(2) = 50; nc(3) = 20;
  lb = uc*nc + 1.0E-10
  npc = nc/10
  lbs = lb/npc
  ncell = npc(1)*npc(2)*npc(3)
  nimpact = 400
  klt = 0
  
  allocate(atype(natm),id(natm),p(natm),cnt(ncell))
  allocate(X(natm,3),V(natm,3),F(natm,3),kin(natm),xl(ncell,3),vl(ncell,3))

  cnt = 0
  n=0
  do i=1,npc(1)
      do j=1,npc(2)
          do k=1,npc(3)
              n=n+1
              xl(n,1) = i*lbs(1)-lbs(1)/2.0
              xl(n,2) = j*lbs(2)-lbs(2)/2.0
              xl(n,3) = k*lbs(3)-lbs(3)/2.0
          end do
      end do
  end do

  do k=1,nimpact
      write(cmd,"('ls D/restart.in',I6.6, '* > restart_list.dat')")k
      print*, cmd
      call system (cmd)
      open(33, file='restart_list.dat')
      !do fl=1,1000
      read(33,"(A28)",IOSTAT=reason)fn
      close (33)
      !    if(reason<0) exit
      read(fn(20:28),"(I9)")listlt
      !end do
      print*, listlt
      !stop
  
      do lt=klt,listlt,dlt
          write(fn,"('D/atoms',I9.9,'.dat')")lt
          open(29,file=fn)
          do i=1,natm
              read(29,*)X(i,:), atype(i), p(i), kin(i), V(i,:), F(i,:), id(i)
              ijk = INT(abs(X(i,:))/lbs)
              if(ijk(3) .lt. npc(3)) then
                  dum = ijk(1) + ijk(2)*Npc(1) + ijk(3)*Npc(1)*Npc(2) + 1
                  cnt(dum) = cnt(dum)+1
                  vl(dum,:) = vl(dum,:)+V(i,:)
                  !write(31,"(6E20.10,4I4,I9)")X(i,1:3),vl(dum,:),cnt(dum),ijk(:),dum
              end if
          end do
          close(29)
      end do
      klt = lt
      print*, 'here'
      !if(mod(lt,4000).eq.0) then
      write(fn,"('meanvel/avg_vel_',I7.7,'.dat')")k
      open(30,file=fn)
      do i=1,ncell
          vl(i,:) = vl(i,:)/cnt(i)
          if(cnt(i).eq.0) vl(i,:) = 0.0
          write(30,"(6E20.10)")xl(i,:),vl(i,:)
      end do
      vl = 0.0
      cnt = 0
      close(30)
      !end if
      print*,k


      write(fn,"('D/atoms',I9.9,'.dat')")listlt
      open(29,file=fn)
      do i=1,natm
          read(29,*)X(i,:), atype(i), p(i), kin(i), V(i,:), F(i,:), id(i)
          ijk = INT(abs(X(i,:))/lbs)
          if(ijk(3) .lt. npc(3)) then
              dum = ijk(1) + ijk(2)*Npc(1) + ijk(3)*Npc(1)*Npc(2) + 1
              cnt(dum) = cnt(dum)+1
              vl(dum,:) = vl(dum,:)+V(i,:)
              !write(31,"(6E20.10,4I4,I9)")X(i,1:3),vl(dum,:),cnt(dum),ijk(:),dum
          end if
      end do
      close(29)
      write(fn,"('meanvel2/avg_vel_',I7.7,'.dat')")k
      open(30,file=fn)
      do i=1,ncell
          vl(i,:) = vl(i,:)/cnt(i)
          if(cnt(i).eq.0) vl(i,:) = 0.0
          write(30,"(6E20.10)")xl(i,:),vl(i,:)
      end do
      vl = 0.0
      cnt = 0
      close(30)
      
      
  end do
      
end program velocity_avg

  
