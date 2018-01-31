MODULE temp

  USE prms
  USE data
  USE tools
  USE stats
  USE parallel
  USE par

  IMPLICIT none 

  integer   ::  Ntemp1,Ntemp2,NtempTot1,NtempTot2
  integer,allocatable,dimension(:)   ::  T1list,T2list

CONTAINS

  SUBROUTINE bub_dim (X,V,lt,impact,time)

    integer :: i,j,k,p,q,r,nmax,nxmax,nymax,nzmax,nx,ny,nz,n, ncryst,bigc,id,bi,bf,bfnew,binew,xx,yy,zz,reason,indx,impact,lt
    real :: uc, lx,ly,lz, dc,thresh,thresh2,m,time
    real, dimension(Natm,3) :: X,V
    real, allocatable, dimension (:,:) :: bbl_dtls,blist_vol
    real, allocatable, dimension(:) :: blist
    character(100) :: fn,rn,cmd

    call updateroot(X)
    call updateroot(V)

    if(myid .eq. 0) then

        uc = 5.431073E-10
        lx=uc*50.0; ly=uc*50.0; lz=uc*20.0;
        dc = uc;
        m = 46.637063E-27
        thresh = 0.3
        thresh2 = 1
        nxmax = ceiling(lx/dc);
        nymax = ceiling(ly/dc);
        nzmax = ceiling(lz/dc);
        nmax = nxmax*nymax*nzmax
        ncryst = dc*dc*dc/(uc*uc*uc)*8.0
        
        if(.not.allocated(bbl_dtls)) then
            allocate(bbl_dtls(nmax,11),blist(nmax),blist_vol(nmax,4))
        end if
        
        bbl_dtls = 0.0;
        blist = 0.0
        blist_vol = 0.0
        
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
        
        do nx=1,nxmax
            do ny=1,nymax
                do nz=1,nzmax
                    n=(nx-1)*nymax*nzmax + (ny-1)*nzmax + nz;
                    bbl_dtls(n,1:3) = (/ nx, ny, nz /);
                    bbl_dtls(n,5) = bbl_dtls(n,4)/ncryst;
                    if(bbl_dtls(n,4).gt.0.0) bbl_dtls(n,9) = bbl_dtls(n,8)/bbl_dtls(n,4)
                    if(bbl_dtls(n,5).lt.thresh) then
                        bbl_dtls(n,6) = 0.0
                    end if
                    if(bbl_dtls(n,5).gt.thresh) bbl_dtls(n,6) = 1.0
                end do
            end do
        end do
        
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
        write(fn,"(I4.4,'_det/bubble_details_',I4.4,'_',I9.9,'.dat')")impact,impact,lt
        open(31,file=fn)
        do nx=1,nxmax
            do ny=1,nymax
                do nz=1,nzmax
                    n=(nx-1)*nymax*nzmax + (ny-1)*nzmax + nz;
                    if(bbl_dtls(n,6).eq.0) then
                        write(31,"(11E15.6)")bbl_dtls(n,:)
                    end if
                end do
            end do
        end do
        close(31)
        
        write(fn,"(I4.4,'_dim/bubble_vol_',I4.4,'_',I9.9,'.dat')")impact,impact,lt
        open(32,file=fn)
        do i=1,id
            if(blist_vol(i,4).ge.thresh2) then
                write(32,"(6E16.8)") blist_vol(i,1:3)/blist_vol(i,4), blist_vol(i,4), blist_vol(i,4)**(1./3.), time
                !write(*,"(6E15.6)") blist_vol(i,1:3)/blist_vol(i,4), blist_vol(i,4), blist_vol(i,4)**(1./3.), time
            end if
        end do
        close(32)

    end if
        
  END SUBROUTINE bub_dim
  
  SUBROUTINE controlsidetemp (X,V,F, Ttar)
    real, dimension(Natm,3)     :: X,V,F
    real                        :: tempside, tempside_l, s, Ttar, nm
    integer                     :: kk_l, ierr, kk, i, ii
    integer, dimension(Nl)      :: tempside_list

    tempside_list = 0
    kk_l = 0
    tempside_l = 0.0
    tempside = 0.0
    nm = 1.0e-9
    do ii=1,Nl
        i = il(ii)
        if(   ((X(i,1).lt.sidewidth*nm) .or. (X(i,1).gt.(Lb(1)-sidewidth*nm)))  .or. ((X(i,2).lt.sidewidth*nm) .or. (X(i,2).gt.(Lb(2)-sidewidth*nm)))  .and.  ((X(i,3).lt.outz(1)*Lb(3)/10.0) .or. (X(i,3).gt.outz(4)*Lb(3)/10.0))) then
            tempside_l = tempside_l + mass(i)*(V(i,1)**2   + V(i,2)**2   + V(i,3)**2  )
            kk_l = kk_l+1
            tempside_list(kk_l) = i
            
        end if
    end do

    !call MPI_ALLREDUCE(kk_l, kk, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    !tempside_l = tempside_l/3./REAL(kk)/kB
    !call MPI_ALLREDUCE(tempside_l, tempside, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Disabling previous lines, because scaling in local processors, no need for global scaling.
    if(kk_l .gt. 0) then
        tempside_l = tempside_l/3./REAL(kk_l)/kB
    
        !s = SQRT(1. + Ts/Tau*(Ttar/tempside-1.))
        s = SQRT(1. + Ts/Tau*(Ttar/tempside_l-1.))
    
        do ii = 1,kk_l
            i = tempside_list(ii)
            !if(myid.eq.0) write(51,*)V(i,:)
            V(i,:) = s*V(i,:)
            !if(myid.eq.0) write(52,*)V(i,:)
        end do
    end if
    
  END SUBROUTINE controlsidetemp

  SUBROUTINE cntsputter(X,V,F,lt,impact)

    real, dimension(Natm,3)     :: X,V,F
    integer                     :: i,ii,ierr,lt,impact, sputtercnt, sputtercnt_l
    real                        :: buffer

    sputtercnt_l = 0
    sputtercnt = 0
    buffer = 0.7

    do ii=1,Nl
        i = il(ii)
        if(X(i,3) .gt. ((atomz(1)-buffer)*Lb(3)/10.0) .and. X(i,3) .lt. ((atomz(1)+buffer)*Lb(3)/10.0) .and. atype(i).lt.3) then
            sputtercnt_l = sputtercnt_l + 1
        end if
    end do
    call MPI_REDUCE(sputtercnt_l,sputtercnt,1,MPI_INTEGER,MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if(myid .eq. 0) then
        open(91, file='sputtercnt.dat', POSITION='APPEND')
        write(91,"(2I9)") impact, sputtercnt
        close(91)
    end if



  END SUBROUTINE CNTSPUTTER
  
  
  SUBROUTINE frzsptrdatms (X,V,F,lt)
    real, dimension(Natm,3)     :: X,V,F
    integer                     :: i,ii,ierr,lt
    real                        :: buffer

    sputter_index_l = 0
    buffer = 0.3
    !sputter_index = 0
    do ii=1,Nl
        i = il(ii)
        if(   ((X(i,3).gt.outz(1)*Lb(3)/10.0) .and. (X(i,3).lt.outz(2)*Lb(3)/10.0))   ) then
            !sputter_index_l(i) = 1
            X(i,3) = atomz(1)*Lb(3)/10.0
            V(i,:) = 0.0
            F(i,:) = 0.0
            !    call si_nlist(X,lt)
        elseif (  ((X(i,3).gt.outz(3)*Lb(3)/10.0) .and. (X(i,3).lt.outz(4)*Lb(3)/10.0))    ) then
            !sputter_index_l(i) = 1
            X(i,3) = atomz(2)*Lb(3)/10.0
            V(i,:) = 0.0
            F(i,:) = 0.0
            !    call si_nlist(X,lt)
        elseif( ((X(i,3).gt.(atomz(1)-buffer)*Lb(3)/10.0) .and. (X(i,3).lt.(atomz(1)+buffer)*Lb(3)/10.0))  ) then
        !    sputter_index_l(i) = 1
            V(i,:) = 0.0
            F(i,:) = 0.0
        elseif( ((X(i,3).gt.(atomz(2)-buffer)*Lb(3)/10.0) .and. (X(i,3).lt.(atomz(2)+buffer)*Lb(3)/10.0))  ) then
            !    sputter_index_l(i) = 1
            V(i,:) = 0.0
            F(i,:) = 0.0
        elseif( ((X(i,3).gt.(ionz-buffer)*Lb(3)/10.0) .and. (X(i,3).lt.(ionz+buffer)*Lb(3)/10.0))  ) then
        !    sputter_index_l(i) = 1
        !    X(i,3) = ionz/10.0*Lb(3)
            V(i,:) = 0.0
            F(i,:) = 0.0
        !    call si_nlist(X,lt)
        end if
    end do
    !call MPI_ALLREDUCE(sputter_index_l, sputter_index, Natm, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  END SUBROUTINE frzsptrdatms
  
  
  
  SUBROUTINE compute_TProf(V,Tbar)
    real, dimension(Natm,3)  :: V
    real, dimension(Nlat(1))    :: Tbar
    integer, dimension(Nlat(1)) :: nTbar
    real, dimension(Nlat(1))    :: Tbarl
    integer, dimension(Nlat(1)) :: nTbarl

    integer                  :: i,ii
    real                     :: Tatm
    integer                  :: ierr

    nTbarl = 0
    Tbarl = 0.

    do ii = 1,Nl
       i = il(ii)
       nTbarl(Lijk(i,1)) = nTbarl(Lijk(i,1)) + 1
       Tatm = mass(i)*(V(i,1)**2 + V(i,2)**2 + V(i,3)**2)/(3.*kB)
       Tbarl(Lijk(i,1)) = Tbarl(Lijk(i,1))+Tatm
    end do

    call MPI_REDUCE(Tbarl,  Tbar,  Nlat(1), MPI_REAL8,   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(nTbarl, nTbar, Nlat(1), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    Tbar = Tbar/REAL(nTbar)

  END SUBROUTINE compute_TProf

  SUBROUTINE TempList(X)
    real, dimension(Natm,3)  :: X

    integer     :: ierr
    integer     :: i

    real,dimension(Natm)    :: T1tmp,T2tmp

    if (myid.eq.0) write(out_unit,*) "INIT TEMP CONTROL"

    Ntemp1 = 0
    Ntemp2 = 0

    if (whichtemp.eq.1) then
       do i = 1,Natm
          if (Lijk(i,1).le.Nltemp1+1.and. &
                                  Lijk(i,1).gt.1.and.P(i).eq.myid) then 
             Ntemp1 = Ntemp1 + 1
             T1tmp(Ntemp1) = i
          else if (Lijk(i,1).gt.Nlat(1)-Nltemp2-1.and. &
                                  Lijk(i,1).lt.Nlat(1).and.P(i).eq.myid) then
             Ntemp2 = Ntemp2 + 1
             T2tmp(Ntemp2) = i
          end if
       end do
    else if (whichtemp.eq.2) then
       do i = 1,Natm
          if (Lijk(i,1).le.Nltemp1.and.P(i).eq.myid) then  
             Ntemp1 = Ntemp1 + 1
             T1tmp(Ntemp1) = i
          else if (Lijk(i,1).gt.Nlat(1)/2.and.Lijk(i,1).le.Nlat(1)/2+Nltemp2.and.P(i).eq.myid) then
             Ntemp2 = Ntemp2 + 1
             T2tmp(Ntemp2) = i
          end if
       end do
    end if

    allocate (T1list(Ntemp1),T2list(Ntemp2))

    T1list = T1tmp(1:Ntemp1)
    T2list = T2tmp(1:Ntemp2)

    call MPI_ALLREDUCE(Ntemp1,NtempTot1,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(Ntemp2,NtempTot2,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    !call writeTsten commented by Kallol

  END SUBROUTINE TempList

  SUBROUTINE writeTsten

    character(30)  :: fn

    write(fn,"('Tstencil.1.',I3.3)")myid
    open(1,file=fn)
    write(1,"(I10)")Ntemp1,T1list
    close(1)

    write(fn,"('Tstencil.2.',I3.3)")myid
    open(1,file=fn)
    write(1,"(I10)")Ntemp2,T2list
    close(1)

    write(fn,"('Tstencil.0.',I3.3)")myid
    open(1,file=fn)
    write(1,"(I10)")NtempTot1,NtempTot2
    close(1)

  END SUBROUTINE writeTsten

  SUBROUTINE TempControl(time,X,V,lt)
    real                     :: time
    integer                  :: lt
    real, dimension(Natm,3)  :: X,V

    integer     :: i
    integer     :: ierr
    real        :: Temp1,Temp2  ! temps of zone 1/2
    real, dimension(4)    :: KE12ab,KE12abl   ! KE of zone 1/2 before and after

    if (lt.lt.Tinit) then
       if (lt.lt.Tinit-1) then
          call temperature(V,Temp1)
          call adjusttemp(V,Temp1,(Ttar1+Ttar2)/2.)
       end if
       call constrainmom(V,mass)

    else if (whichtemp.gt.0) then

       call zonetemperature(V,Ntemp1,NtempTot1,T1list,Temp1,KE12ab(1))
       call zonetemperature(V,Ntemp2,NtempTot2,T2list,Temp2,KE12ab(2))

       if (whichthrm.eq.2) then
          call adjusttempB(V,Ntemp1,T1list,Temp1,Ttar1)
          call adjusttempB(V,Ntemp2,T2list,Temp2,Ttar2)
       else if (whichthrm.eq.1) then
          call adjusttempAnd(V,Ntemp1,T1list,Temp1,Ttar1)
          call adjusttempAnd(V,Ntemp2,T2list,Temp2,Ttar2)
       end if

       call zonetemperature(V,Ntemp1,NtempTot1,T1list,Temp1,KE12ab(3))
       call zonetemperature(V,Ntemp2,NtempTot2,T2list,Temp2,KE12ab(4))

       call MPI_REDUCE(KE12ab,KE12abl,4,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       if (zen_out.eq.1.and.myid.eq.0) then
          write(zen_unit,"(3E20.10)") time,KE12abl(3)-KE12abl(1),KE12abl(4)-KE12abl(2)
       end if

    end if

  END SUBROUTINE TempControl


  SUBROUTINE adjusttempAnd(V,Ntemp,Tlist,Temp,Ttar)
    real, dimension(Natm,3)  :: V

    integer                  :: Ntemp
    integer,dimension(Ntemp) :: Tlist    
    real                     :: Temp
    real                     :: Ttar

    real                     :: s
    integer                  :: i,ii

    integer                  :: l
    real, dimension(3)       :: Vb
    
    real                     :: ran1

    real   :: PE,Bltz


    do i = 1,Ntemp                       
       if (ran1(ranseed).lt.nu) then          
 
          ii = Tlist(i)
 
          do l = 1,3
             do 
                Vb(l) = -Vmax + 2.*Vmax*ran1(ranseed)
                Bltz =EXP(-mass(ii)*Vb(l)*Vb(l)/(2.*kB*Ttar))
                if (ran1(ranseed).lt.Bltz) exit
             end do
          end do          
 
          V(ii,:) = Vb

       end if

    end do

  END SUBROUTINE adjusttempAnd

 SUBROUTINE adjusttempB(V,Ntemp,Tlist,Temp,Ttar)
    real, dimension(Natm,3)  :: V

    integer                  :: Ntemp
    integer,dimension(Natm)  :: Tlist    
    real                     :: Temp
    real                     :: Ttar

    real                     :: s
    integer                  :: i,ii

    s = SQRT(1. + Ts/Tau*(Ttar/Temp-1.))

!    print *,myid,": ",Ttar,Temp,s
    do i = 1,Ntemp
       ii = Tlist(i)
       V(ii,:) = s*V(ii,:)
    end do

  END SUBROUTINE adjusttempB

  SUBROUTINE adjusttemp(V,Temp,Ttar)
    real, dimension(Natm,3)  :: V
   
    real                     :: Temp
    real                     :: Ttar

    real                     :: s
    
    integer                  :: i,ii

    s = SQRT(1. + Ts/Tau*(Ttar/Temp-1.))
    
    !print *,myid,": ",Ttar,Temp,s
    !if (myid .eq. 0) print*,s, Temp
    do ii = 1,Nl
        i = il(ii)
        if(atype(i).lt.3)  V(i,:) = s*V(i,:)
    end do

  END SUBROUTINE adjusttemp

  SUBROUTINE zonetemperature(V,Ntemp,NtempTot,Tlist,Temp,KE)
    real, dimension(Natm,3)  :: V

    integer                  :: Ntemp,NtempTot
    integer,dimension(Natm)  :: Tlist    
    real                     :: Temp,KE,KE_l

    integer                  :: ierr
    integer                  :: i,ii

!                       3
!       3/2 kB T   =   SUM ( 1/2 m v_i^2 )
!                      i=1

    KE_l = 0.
    do i = 1,Ntemp
       ii = Tlist(i)
       KE_l = KE_l + mass(ii)*(V(ii,1)**2 + V(ii,2)**2 + V(ii,3)**2  )
    end do

    call MPI_ALLREDUCE(KE_l,  KE,  1, MPI_REAL8,   MPI_SUM, MPI_COMM_WORLD, ierr)

    Temp = KE/3./REAL(NtempTot)/kB
    KE = KE/2.

  END SUBROUTINE zonetemperature
  
!  Compute the mean temperature 
  SUBROUTINE temperature(V,Temp)
    real, dimension(Natm,3)  :: V
    real                     :: Temp,Temp_l
    integer                  :: i,ii

    integer                  :: ierr

!                       3
!       3/2 kB T   =   SUM ( 1/2 m v_i^2 )
!                      i=1

    Temp_l = 0.
    do ii = 1,Nl
        i = il(ii)
        if(atype(i).lt.3) then
            Temp_l = Temp_l + mass(i)*(V(i,1)**2   + V(i,2)**2   + V(i,3)**2  )
        end if
    end do
    Temp_l = Temp_l/3./REAL(Nsg)/kB
    !Temp_l = Temp_l/3./REAL(Nl)/kB  !Kallol

    call MPI_ALLREDUCE(Temp_l, Temp, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  END SUBROUTINE temperature
  
END MODULE temp

