program buuble

  integer :: i,j,k,p,q,r, natm, nimpact,nmax,nxmax,nymax,nzmax,nx,ny,nz,n, ncryst, Nt0,ions,vcount,bigc,id,bi,bf,bfnew,binew,xx,yy,zz,reason
  real :: uc, ev,kb,lx,ly,lz, dc,pi,m,restart_time,thresh,thresh2
  real, allocatable, dimension(:,:) :: X,V,bbl_dtls,blist_vol
  real, allocatable, dimension(:) :: mass,blist
  integer, allocatable, dimension(:) :: atype
  character(100) :: fn,rn,cmd

  natm = 401000
  nimpact = 750
  ev = 1.60217646e-19;
  kb = 1.3806503e-23;
  uc = 5.431073E-10
  lx=uc*50.0; ly=uc*50.0; lz=uc*20.0;
  dc = uc
  pi = 3.141592654
  m = 46.637063E-27
  thresh = 0.3
  thresh2 = 1
  nxmax = ceiling(lx/dc);
  nymax = ceiling(ly/dc);
  nzmax = ceiling(lz/dc);
  nmax = nxmax*nymax*nzmax
  ncryst = dc*dc*dc/(uc*uc*uc)*8.0
  print*, nxmax, nymax, nzmax
  
  allocate(X(natm,3), V(natm,3),bbl_dtls(nmax,11),blist(nmax),blist_vol(nmax,4),atype(natm),mass(natm))

  do impact=1,nimpact
      
  
      bbl_dtls = 0.0;
      blist = 0.0
      blist_vol = 0.0

      write(cmd,"('ls D/restart.in00',I4.4,'_* > restart_list.dat')")impact
      call system(cmd)
      open(33,file='restart_list.dat')
      read(33,"(A50)",IOSTAT=reason)rn
      close(33)
      print*,trim(rn)
      open(30,file=rn,form='UNFORMATTED',status='OLD')
      read(30) Nt0, restart_time
      read(30)X,V,mass,atype
      read(30)ions
      close(30)
      
      do i=1,natm
          if(X(i,3).lt.lz .and. X(i,3).gt.0.0) then
              nx = ceiling(X(i,1)/dc)
              ny = ceiling(X(i,2)/dc)
              nz = ceiling(X(i,3)/dc)
              if(nx .le. 0) nx = 1
              if(ny .le. 0) ny = 1
              if(nz .le. 0) nz = 1
              if(nx .gt. nxmax) nx = nxmax
              if(ny .gt. nymax) ny = nymax
              if(nz .gt. nzmax) nz = nzmax
              if(nx*ny*nz .gt. nmax) print*, 'ERROR!'
              n = (nx-1)*nymax*nzmax + (ny-1)*nzmax + nz
              
              bbl_dtls(n,4) = bbl_dtls(n,4) + 1      ! count atoms
              bbl_dtls(n,7) = 0                      ! atom id
              bbl_dtls(n,8) = bbl_dtls(n,8) + m*(sum(V(i,:)**2))/3.0/kb ! temperature addition
              
          end if
      end do
      
      !open(31,file='details_bubble.dat')
      vcount = 0
      do nx=1,nxmax
          do ny=1,nymax
              do nz=1,nzmax
                  n=(nx-1)*nymax*nzmax + (ny-1)*nzmax + nz;
                  bbl_dtls(n,1:3) = (/ nx, ny, nz /);
                  bbl_dtls(n,5) = bbl_dtls(n,4)/ncryst;
                  if(bbl_dtls(n,4).gt.0.0) bbl_dtls(n,9) = bbl_dtls(n,8)/bbl_dtls(n,4)
                  if(bbl_dtls(n,5).lt.thresh) then
                      bbl_dtls(n,6) = 0.0
                      vcount = vcount+1
                  end if
                  if(bbl_dtls(n,5).gt.thresh) bbl_dtls(n,6) = 1.0
                  
                  !write(31,"(11E15.6)")bbl_dtls(n,:)
                  !write(*,"(11E15.6)")bbl_dtls(n,:)
              end do
          end do
      end do
      !close(31)
      
      bigc = 0
      id = 0
      bi=1
      bf=0
      bfnew=1
      binew=1
      
      do k=1,nmax
          if(bbl_dtls(k,6).eq.0.0 .and. bbl_dtls(k,10).eq.0.0) then
              bigc=bigc+1
              bf = bf+1
              id = id+1
              bbl_dtls(k,7) = id
              blist(bigc) = k
              
              !blist_vol(id,1:3) = blist_vol(id,1:3) + bbl_dtls(k,1:3)
              !blist_vol(id,4) = blist_vol(id,4)+1
              
              !print*, k
              
              do
                  
                  do i=bi,bf
                      !binew = i+1
                      !bfnew = bf
                      indx = blist(i)
                      if(bbl_dtls(indx,10).eq.0) then
                          bbl_dtls(indx,10)=1;
                          bbl_dtls(indx,11)=1;
                          blist_vol(id,1:3) = blist_vol(id,1:3) + bbl_dtls(indx,1:3)
                          blist_vol(id,4) = blist_vol(id,4)+1
                          
                          xx = bbl_dtls(indx,1);
                          yy = bbl_dtls(indx,2);
                          zz = bbl_dtls(indx,3);
                          do p=-1,1
                              do q=-1,1
                                  do r=-1,1
                                      nx=xx+p; ny=yy+q; nz=zz+r;
                                      if(nx.le.nxmax .and. nx.gt.0 .and. ny.le.nymax .and. ny.gt.0 .and. nz.le.nzmax .and. nz.gt.0) then
                                          if(p.eq.0 .and. q.eq.0 .and. r.eq.0) cycle
                                          n = (nx-1)*nymax*nzmax + (ny-1)*nzmax + nz;
                                          if(bbl_dtls(n,6).eq.0 .and. bbl_dtls(n,10).eq.0 .and. bbl_dtls(n,11).eq.0) then
                                              !bfnew=bfnew+1
                                              bigc = bigc+1;
                                              !print*, bigc, bfnew, id
                                              blist(bigc) = n;
                                              bbl_dtls(n,7) = id;
                                              bbl_dtls(n,11)=1;
                                          end if
                                      end if
                                  end do
                              end do
                          end do
                      end if
                  end do
                  
                  bi = bf+1
                  bf = bigc
                  !print*, bi,bf,bigc,id
                  if(bi.gt.bf) exit
                  
                  
              end do
              
          end if
      end do
      
      ! disabling for now
      !open(31,file='details_bubble.dat')
      !do nx=1,nxmax
      !    do ny=1,nymax
      !        do nz=1,nzmax
      !            n=(nx-1)*nymax*nzmax + (ny-1)*nzmax + nz;
      !            write(31,"(11E15.6)")bbl_dtls(n,:)
      !        end do
      !    end do
      !end do
      !close(31)

      write(fn,"('bubble_dim/bubble_vol_',I4.4,'.dat')")impact
      open(32,file=fn)
      do i=1,id
          if(blist_vol(i,4).ge.thresh2) then
              write(32,"(5E15.6)") blist_vol(i,1:3)/blist_vol(i,4), blist_vol(i,4), blist_vol(i,4)**(1./3.)
              write(*,"(5E15.6)") blist_vol(i,1:3)/blist_vol(i,4), blist_vol(i,4), blist_vol(i,4)**(1./3.)
          end if
      end do
      close(32)
      
      
  end do

  
  
end program buuble
