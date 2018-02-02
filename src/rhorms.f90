program act

  integer :: i,j,k,p,q,r,nrmax,nmax,nxmax,nymax,nzmax,nx,ny,nz,n, ncryst,bigc,id,bi,bf,bfnew,binew,xx,yy,zz,reason,indx,impact,lt, natm, nimpact, ions, l, Nt0, radid
  integer :: newnx, newny, newnz,newn, idnew, bigci, idi, bii, bfi, cnt, nsurf
  real :: uc, lx,ly,lz, dc,thresh,thresh2,m,time, restart_time,kb, lzdiv, lzmax, act_surf, thresh3, wsurf
  real, allocatable,dimension(:,:) :: X,V, multiplier
  real, allocatable, dimension (:,:) :: bbl_dtls,blist_vol,clist_vol
  real, allocatable, dimension(:) :: blist, mass, blisti
  integer,allocatable,dimension(:):: atype
  character(100) :: fn,rn,cmd
  integer, allocatable,dimension(:,:) :: surfid
  real :: w1, w2, w3, w4, w5, w6, w7, w8, w9
  real,allocatable,dimension(:)::sumden, sumtemp, sumsqden, sumsqtemp, stdden, stdtemp, sumden2, sumtemp2
  integer, allocatable,dimension(:)::sumnum, stdnum

  natm = 5121000
  nimpact = 1000
  uc = 5.431073E-10
  lx=uc*80.0; ly=uc*80.0; lz=uc*100.0;
  nrmax = ceiling(lx/2.0/1.0e-9);
  lzdiv = uc*800.0; lzmax = uc*1000.0
  dc = uc*1.0;
  m = 46.637063E-27
  thresh = 0.3 !0.02
  thresh2 = 4.0
  kb = 1.381E-23
  nxmax = ceiling(lx/dc);
  nymax = ceiling(ly/dc);
  nzmax = ceiling(lz/dc);
  nmax = nxmax*nymax*nzmax
  ncryst = dc*dc*dc/(uc*uc*uc)*8.0
  thresh3 = lx/dc*ly/dc*uc*20.0/dc/3.0
  print*, nmax,ncryst

  allocate(X(natm,3),V(natm,3),mass(natm),atype(natm))
  allocate(bbl_dtls(nmax,8),blist(nmax),blist_vol(nmax,4), blisti(nmax),clist_vol(nmax,4))
  allocate(sumden(nrmax), sumtemp(nrmax), sumsqden(nrmax), sumsqtemp(nrmax), stdden(nrmax), stdtemp(nrmax), sumden2(nrmax), sumtemp2(nrmax), sumnum(nrmax), stdnum(nrmax))

  open(41,file='rmsrho.dat')
  open(42,file='rmstemp.dat')

  do l =1,nimpact !00,100 !,nimpact !1,nimpact

      write(cmd,"('ls D/restart.in',I6.6, '* > restart_list.dat')")l
      call system (cmd)
      open(33, file='restart_list.dat')
      read(33,"(A28)",IOSTAT=reason)fn
      close (33)
      print*, trim(fn)
      open(res_unit,file=trim(fn),form='UNFORMATTED',status='OLD')
      read(res_unit)Nt0,restart_time
      read(res_unit)X,V,mass,atype
      read(res_unit)ions
      lt = Nt0
      time = restart_time
      X = X+0.2E-10

      bbl_dtls = 0.0;

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
              bbl_dtls(n,6) = bbl_dtls(n,6) + m*(sum(V(i,:)**2))/3.0/kb ! temperature addition

          end if
      end do

      sumden = 0.0; sumtemp = 0.0; sumnum = 0;
      sumsqden = 0.0; sumsqtemp = 0.0;
      stdden = 0.0; stdtemp = 0.0; stdnum = 0;
      sumden2 = 0.0; sumtemp2 = 0.0;

      !write(rn,"('for_rhorms/bubble_details_',I4.4,'_',I9.9,'.dat')")l,lt
      !open(43,file=rn)
      do nx=1,nxmax
          do ny=1,nymax
              do nz=1,nzmax
                  n=(nx-1)*nymax*nzmax + (ny-1)*nzmax + nz;
                  bbl_dtls(n,1:3) = (/ (nx-0.5)*dc, (ny-0.5)*dc, (nz-0.5)*dc /);
                  bbl_dtls(n,5) = bbl_dtls(n,4)/ncryst;
                  if(bbl_dtls(n,4).gt.0.0) bbl_dtls(n,7) = bbl_dtls(n,6)/bbl_dtls(n,4)
                  bbl_dtls(n,8) = sqrt((bbl_dtls(n,1)-lx/2.0)**2 + (bbl_dtls(n,2)-ly/2.0)**2)
                  !write(43,"(8E15.6)")bbl_dtls(n,:)

                  radid = ceiling(bbl_dtls(n,8)/1.0e-9)
                  if(radid .eq. 0) radid = 1

                  if(radid .le. nrmax) then
                      sumden(radid) = sumden(radid) + bbl_dtls(n,5)
                      sumsqden(radid) = sumsqden(radid) + (bbl_dtls(n,5))**2
                      sumtemp(radid) = sumtemp(radid) + bbl_dtls(n,7)
                      sumsqtemp(radid) = sumsqtemp(radid) + (bbl_dtls(n,7))**2
                      sumnum(radid) = sumnum(radid)+1
                  end if

              end do
          end do
      end do


      do i=1,nrmax
          do j = i,1,-1
              stdden (i) = stdden(i) + sumsqden(j)
              sumden2(i) = sumden2(i) + sumden(j)
              stdtemp(i) = stdtemp(i) + sumsqtemp(j)
              sumtemp2(i) = sumtemp2(i) + sumtemp(j)
              stdnum (i) = stdnum(i) + sumnum(j)
          end do

          stdden (i) = sqrt(stdden(i)/stdnum(i) - (sumden2(i)/stdnum(i))**2)
          stdtemp(i) = sqrt(stdtemp(i)/stdnum(i) - (sumtemp2(i)/stdnum(i))**2)
      end do

      write(41,"(I5, <nrmax>G15.5)")l,stdden(:)
      write(42,"(I5, <nrmax>E15.5)")l,stdtemp(:)

  end do

  close(41)
  close(42)



end program act
