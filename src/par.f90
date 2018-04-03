!------------------------------------------------------------------------------
!
! `PAR` Source File
!
! par.f90 source function file. This file contains functions for assigning
! atoms to processors, creating maps for data assigned to each processor,
! and distributing neighborlist information to relevant processors
!
!------------------------------------------------------------------------------

MODULE par

  USE prms
  USE parallel

  IMPLICIT none

  integer :: Nl                     ! number of atoms on local processor
  integer :: Nsr,Nsra               ! number of processor to send force

  integer, allocatable, dimension(:) :: il     ! global index of local atoms
  integer, allocatable, dimension(:) :: Psr    ! processors to communicate with
  integer, allocatable, dimension(:) :: Psra   ! processors to communicate with
  integer, allocatable, dimension(:)   :: Nisf ! number of forces to send to each proc
  integer, allocatable, dimension(:)   :: Nisa ! number of atoms to send to each proc
  integer, allocatable, dimension(:,:) :: Isf  ! nlist2 of Fi's to send
  integer, allocatable, dimension(:,:) :: Isa  ! nlist2 of atoms to send
  integer, allocatable, dimension(:)   :: P, P_l, P_big      ! local processor
  integer, allocatable, dimension(:,:) :: P_send

  integer :: Ms   ! number of s/r
  integer, allocatable, dimension(:,:) :: procSR
  real, dimension(3)      :: Lbs
  real                                 :: buff

CONTAINS

  SUBROUTINE initparaops(X)

    real, dimension(Natm,3) :: X
    integer :: ierr

    Lbs(1) = (Lb(1)        + 1.E-10)/REAL(Npc(1))
    Lbs(2) = (Lb(2)        + 1.E-10)/REAL(Npc(2))
    Lbs(3) = (Lb(3)/10.0   + 1.E-10)/REAL(Npc(3))
    buff = listfac*al*MAX(sigmaGe,sigmaSi)

    allocate (P(Natm), P_l(Natm), P_big(Natm), P_send(Nl,Np))
    call assignprocstaub(X)
    call findlocals

  END SUBROUTINE initparaops

  SUBROUTINE initparaops_l(X,V,F)

    real, dimension(Natm,3) :: X,V,F
    integer :: ierr

    P_l = 0
    if (allocated(P_send)) then
        deallocate(P_send)
        allocate(P_send(Nl,Np))
    end if
    call assignprocs_l(X,V,F)
    call findlocals

  END SUBROUTINE initparaops_l

  SUBROUTINE closeparaops

    deallocate(il,Psr,Nisf,Isf,Isa)

  END SUBROUTINE closeparaops

  SUBROUTINE findlocals

    integer  :: i,ierr
    integer, dimension(Natm) :: il_l

    Nl = 0
    do i = 1,Natm
       if (P(i).eq.myid) then
          Nl = Nl + 1
          il_l(Nl) = i
       end if
   end do
   if (allocated(il)) deallocate(il)
   allocate(il(Nl))
   il = il_l(1:Nl)

   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  END SUBROUTINE findlocals

  SUBROUTINE updateglobal(G)
    real, dimension(Natm,3) :: G

    integer      :: i
    integer      :: ierr

    do i = 1,Natm
        call MPI_BCAST(G(i,1), 1, MPI_REAL8, P(i), MPI_COMM_WORLD, ierr)
        call MPI_BCAST(G(i,2), 1, MPI_REAL8, P(i), MPI_COMM_WORLD, ierr)
        call MPI_BCAST(G(i,3), 1, MPI_REAL8, P(i), MPI_COMM_WORLD, ierr)
    end do

  END SUBROUTINE updateglobal

  SUBROUTINE updateroot(G)
    real, dimension(Natm,3) :: G
    real, dimension(3,Natm) :: Gl
    integer, dimension(Natm) :: indx
    integer :: i,ii,Ni,j, ierr
    integer, dimension(MPI_STATUS_SIZE) :: status

    if(myid .ne. 0) then
        do ii=1,Nl
            i=il(ii)
            Gl(:,ii) = G(i,:)
        end do
        call MPI_SEND(Nl,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
        call MPI_SEND(il,Nl,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(Gl,3*Nl,MPI_REAL8,0,2,MPI_COMM_WORLD,ierr)
    else
        do j=1,Np-1
            call MPI_RECV(Ni,1,MPI_INTEGER,j,0,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(indx,Ni,MPI_INTEGER,j,1,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(Gl,3*Ni,MPI_REAL8,j,2,MPI_COMM_WORLD,status,ierr)
            do ii=1,Ni
                i=indx(ii)
                G(i,:) = Gl(:,ii)
            end do
        end do
    end if

  END SUBROUTINE updateroot

  SUBROUTINE find_p_big (dummy_x)
    real, dimension(Natm,3) :: dummy_x
    integer                 :: i,j,l,ierr,dum
    integer, dimension(3)   :: ijk

    do l = 1,3
        dummy_x(:,l) = MOD(dummy_x(:,l) + 2.*Lb(l),Lb(l))
    end do
    do i = 1,Natm
        if(P_big(i).lt.0) then
            ijk = INT(abs(dummy_x(i,:))/Lbs) + 1
            do l = 1,3
                if(ijk(l) .gt. Npc(l)) ijk(l) = Npc(l)
            end do
            dum = ijk(1)-1 + (ijk(2)-1)*Npc(1) + (ijk(3)-1)*Npc(1)*Npc(2)
            if(dum.eq.myid) P_big(i)=dum
        end if
    end do
  END SUBROUTINE find_p_big

  SUBROUTINE assignprocstaub(X)

    real, dimension(Natm,3)   :: X, dummy_x
    integer                   :: i,j,l,ierr,maxdim,dum
    integer, dimension(3)     :: ijk
!    real, dimension(3)        :: q,r,s	!modified for dimension by Josh
    integer                   :: q, r, s
!    integer :: ri,qi,si	!Added by Josh
!    q = [-1, 1, 2]
!    r = q
!    s = q
    P_big = -1

    if (myid .eq. 0) then
        do i = 1,Natm
            ijk = INT(abs(X(i,:))/Lbs) + 1
            do l = 1,3
                if(ijk(l) .gt. Npc(l)) ijk(l) = Npc(l)
            end do
            P(i) = ijk(1)-1 + (ijk(2)-1)*Npc(1) + (ijk(3)-1)*Npc(1)*Npc(2)
        end do
    end if
    call MPI_BCAST(P, Natm, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    do i=1,Natm
        if(P(i).eq.myid) P_big(i)=P(i)
    end do
    dummy_x = X
    if(Npc(1).gt.1) then
        do s=-1,1,2
            dummy_x(:,1)=X(:,1)+buff*s
            if(Npc(2).gt.1) then
                do q=-1,1,2
                    dummy_x(:,2)=X(:,2)+buff*q
                    if(Npc(3).gt.1) then
                        do r=-1,1,2
                            dummy_x(:,3)=X(:,3)+buff*r
                            call find_p_big(dummy_x)
                        end do
                    else
                        call find_p_big(dummy_x)
                    end if
                end do
            else
                call find_p_big(dummy_x)
            end if
        end do
    end if

  END SUBROUTINE assignprocstaub


  SUBROUTINE find_p_send(dummy_x, n_send_l, i_send)
    real, dimension(Natm,3)           :: dummy_x
    integer                           :: i,j,l,ierr,maxdim,ii,dum
    integer, dimension(3)             :: ijk
    integer, dimension(Np,Np)         :: n_send_l
    integer, dimension(2*Natm/Np,Np)  :: i_send


    do ii = 1,Nl
        i=il(ii)
        ijk = INT(abs(dummy_x(i,:))/Lbs) + 1
        do l = 1,3
            if(ijk(l) .gt. Npc(l)) ijk(l) = Npc(l)
        end do
        dum = ijk(1)-1 + (ijk(2)-1)*Npc(1) + (ijk(3)-1)*Npc(1)*Npc(2)
        if(P(i).ne.dum .and. P_l(i).ne.dum .and. P_send(ii,dum+1).lt.0) then
            P_send(ii,dum+1) = dum
            n_send_l(myid+1,dum+1) = n_send_l(myid+1,dum+1)+1
            i_send (n_send_l(myid+1,dum+1),dum+1) = i
        end if
    end do
  END SUBROUTINE find_p_send

  SUBROUTINE assignprocs_l(X,V,F)

    real, dimension(Natm,3)         :: X,V,F, dummy_x
    real, dimension(3,2*Natm/Np)      :: Xs, Fs, Vs
    integer                         :: i,j,l,ierr,maxdim,ii,dum
    integer, dimension(3)           :: ijk
    integer, dimension(Np,Np)       :: n_send, n_send_l
    integer, dimension(2*Natm/Np,Np)  :: i_send
    integer, dimension(2*Natm/Np)     :: i_recv
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: q,r,s

    P_big = -1
    P_send = -1

    n_send = 0
    n_send_l = 0
    i_send = 0

    do ii = 1,Nl
        i=il(ii)
        ijk = INT(abs(X(i,:))/Lbs) + 1
        do l = 1,3
            if(ijk(l) .gt. Npc(l)) ijk(l) = Npc(l)
        end do
        P_l(i) = ijk(1)-1 + (ijk(2)-1)*Npc(1) + (ijk(3)-1)*Npc(1)*Npc(2)
        if(P(i).ne.P_l(i)) then
            n_send_l(myid+1,P_l(i)+1) = n_send_l(myid+1,P_l(i)+1)+1
            i_send (n_send_l(myid+1,P_l(i)+1),P_l(i)+1) = i
        end if
    end do

    do i=1,Natm
        if(P(i).eq.myid) P_big(i)=P(i)
    end do
    dummy_x = X
    if(Npc(1).gt.1) then
        do s=-1,1,2
            dummy_x(:,1)=X(:,1)+buff*REAL(s)
            if(Npc(2).gt.1) then
                do q=-1,1,2
                    dummy_x(:,2)=X(:,2)+buff*REAL(q)
                    if(Npc(3).gt.1) then
                        do r=-1,1,2
                            dummy_x(:,3)=X(:,3)+buff*REAL(r)
                            call find_p_big(dummy_x)
                            call find_p_send(dummy_x, n_send_l, i_send)
                        end do
                    else
                        call find_p_big(dummy_x)
                        call find_p_send(dummy_x, n_send_l, i_send)
                    end if
                end do
            else
                call find_p_big(dummy_x)
                call find_p_send(dummy_x, n_send_l, i_send)
            end if
        end do
    end if

    call MPI_ALLREDUCE(P_l,P,Natm, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(n_send_l,n_send,Np*Np, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    !if(myid .eq. 1)then
    !    do i=1,Np
    !        write(44,"(33I6)")n_send(i,:), i
    !    end do
    !end if
    !write(44,"(8I7)")i_send(1,:)

    do i=0,Np-1
        if(i .eq. myid) then
            do j=0,Np-1
                if(j .ne. myid .and. n_send(i+1,j+1) .gt. 0) then
                    do ii=1,n_send(i+1,j+1)
                        Xs(:,ii) = X(i_send(ii,j+1),:)
                        Vs(:,ii) = V(i_send(ii,j+1),:)
                        Fs(:,ii) = F(i_send(ii,j+1),:)
                    end do
                    !print*,myid,'sends',j,i_send(1,j+1),n_send(i+1,j+1)
                    call MPI_SEND(i_send(1,j+1),n_send(i+1,j+1),MPI_INTEGER, j , 0, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(Xs,3*n_send(i+1,j+1),MPI_REAL8, j , 1, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(Vs,3*n_send(i+1,j+1),MPI_REAL8, j , 2, MPI_COMM_WORLD, ierr)
                    call MPI_SEND(Fs,3*n_send(i+1,j+1),MPI_REAL8, j , 3, MPI_COMM_WORLD, ierr)
                end if
            end do
        elseif(n_send(i+1,myid+1) .gt. 0) then
            call MPI_RECV(i_recv,n_send(i+1,myid+1),MPI_INTEGER, i , 0, MPI_COMM_WORLD, status,ierr)
            call MPI_RECV(Xs,3*n_send(i+1,myid+1),MPI_REAL8, i , 1, MPI_COMM_WORLD, status,ierr)
            call MPI_RECV(Vs,3*n_send(i+1,myid+1),MPI_REAL8, i , 2, MPI_COMM_WORLD, status,ierr)
            call MPI_RECV(Fs,3*n_send(i+1,myid+1),MPI_REAL8, i , 3, MPI_COMM_WORLD, status,ierr)
            !print*,myid,'recvs',i,i_recv(1:n_send(i+1,myid+1)),n_send(i+1,myid+1)
            do ii=1,n_send(i+1,myid+1)
                X(i_recv(ii),:) = Xs(:,ii)
                V(i_recv(ii),:) = Vs(:,ii)
                F(i_recv(ii),:) = Fs(:,ii)
            end do
        end if
    end do

    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  END SUBROUTINE assignprocs_l

  SUBROUTINE collectsends

    integer                            :: iw,ll
    integer                            :: Nsrl
    integer, dimension(Np)             :: PsL
    integer, dimension(0:Np-1,0:Np-1)  :: Smat
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer                            :: ierr

    if (myid.eq.0) then

       Smat = 0.
       do ll = 1,Nsr
          Smat(0,Psr(ll)) = ll
       end do
       !print*, 'inside collectsends 1'
       do iw = 1,Np-1
          call MPI_RECV(Nsrl,1,MPI_INTEGER,iw,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(PsL,Nsrl,MPI_INTEGER,iw,1,MPI_COMM_WORLD,status,ierr)
          do ll = 1,Nsrl
             Smat(iw,PsL(ll)) = ll
          end do
       end do
       !write(88,"('     |',19I5)") (ll,ll=0,Np-1)
       !write(88,"('--------------------------------------------------')")
       !do iw = 0,Np-1
       !   write(88,"(I5,'|',19I5)") iw,(Smat(iw,ll),ll=0,Np-1)
       !end do
    else
       call MPI_SEND(Nsr,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
       call MPI_SEND(Psr,Nsr,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
    end if

    call MPI_BCAST(Smat, Np*Np, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call makesendmaps(Smat)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr); call MPI_FINALIZE(ierr); stop

  END SUBROUTINE collectsends

  SUBROUTINE makesendmaps(Smat)
    integer, dimension(0:Np-1,0:Np-1) :: Smat

    integer, dimension(2*Np,2)        :: procSR_l
    integer, dimension(Np)            :: inNs
    integer                           :: ll,lll
    integer                           :: iv,iw,is
    integer                           :: Ms_l
    integer                           :: ierr

    Ms_l = 0
    procSR_l = 0
    do is = 1,2*Np
       ll = 0
       inNs = 0
       do iw = 0,Np-1
          do iv = 0,Np-1
             if (Smat(iv,iw).gt.0) then
                if (notinlist(iw,inNs,ll).and.notinlist(iv,inNs,ll)) then
                   inNs(ll+1) = iw; inNs(ll+2) = iv; ll = ll + 2
                   Ms_l = is
                   if (iv.eq.myid) procSR_l(is,1) = Smat(iv,iw)
                   if (iw.eq.myid) then
                      do lll = 1,Nsra
                         if (Psra(lll).eq.iv) then
                            procSR_l(is,2) = lll
                         end if
                      end do
                   end if
                   Smat(iv,iw) = -2
                end if
             end if
          end do
       end do
    end do

    call MPI_ALLREDUCE(Ms_l,Ms,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

    if (Ms.gt.0) then
        if (allocated(procSR)) deallocate(procSR)
        allocate (procSR(Ms,2))
        procSR = procSR_l(1:Ms,:)
    else
        if (allocated(procSR)) deallocate(procSR)
        allocate (procSR(1,2))
        procSR = 0
        Ms = 0
    end if

    !do is = 1,Ms
    !    write(90+myid,"(I5,' -- ',2I5)")is,procSR(is,:)
    !end do
!!$    do is = 1,Ms
!!$       write(100+myid,"(I5,' -- ',2I5)")is,Psr(procSR(is,1)),Psra(procSR(is,2))
!!$    end do

  END SUBROUTINE makesendmaps

  SUBROUTINE sendposN(X)

    real, dimension(Natm,3) :: X
    real, dimension(3,Nl)    :: Xs

    integer                 :: ll,ii,is
    integer                 :: ierr
    integer                 :: req
    integer, dimension(MPI_STATUS_SIZE) :: status

    do is = 1,Ms

       if (procSR(is,2).ne.0) then
          ll = procSR(is,2)
          do ii = 1,NIsa(ll)
             Xs(:,ii) = X(Isa(ii,ll),:)
          end do
          call MPI_SEND(Xs,3*NIsa(ll),MPI_REAL8,Psra(ll),0,MPI_COMM_WORLD,ierr)
       else if (procSR(is,1).ne.0) then
          ll = procSR(is,1)
          call MPI_RECV(Xs,3*NIsf(ll),MPI_REAL8,Psr(ll),0,MPI_COMM_WORLD,status,ierr)
          do ii = 1,NIsf(ll)
             X(Isf(ii,ll),:) = Xs(:,ii)
          end do
       end if
    end do

  END SUBROUTINE sendposN

  SUBROUTINE combfN(F, lt)

    real, dimension(Natm,3)   :: F
    real, dimension(3,Nl)     :: Fs, Fs_recv, Fs_send
    integer                 :: ll,ii,is,lt
    integer                 :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: status

    Fs = 0.0 ! kallol
    do is = 1,Ms

        if (procSR(is,1).ne.0) then
            ll = procSR(is,1)
            do ii = 1,NIsf(ll)
                Fs(:,ii) = F(Isf(ii,ll),:)
            end do
            call MPI_SEND(Fs,3*NIsf(ll),MPI_REAL8,Psr(ll),0,MPI_COMM_WORLD,ierr)

        else if (procSR(is,2).ne.0) then
            ll = procSR(is,2)
            call MPI_RECV(Fs,3*NIsa(ll),MPI_REAL8,Psra(ll),0,MPI_COMM_WORLD,status,ierr)
            do ii = 1,NIsa(ll)
                F(Isa(ii,ll),:) = F(Isa(ii,ll),:) + Fs(:,ii)
            end do
        end if

    end do

  END SUBROUTINE combfN

  FUNCTION notinlist(l,P,N) RESULT(notinit)
    integer  :: l,N
    integer, dimension(N) :: P

    logical  :: notinit

    integer  :: i

    notinit = .TRUE.
    if (N.eq.0) then
        notinit = .TRUE.
    else
        do i = 1,N
            if (P(i).eq.l) notinit = .FALSE.
        end do
    end if

  END FUNCTION notinlist

!findIsa does something with passing information for pairs or triplets that contain
!atoms from more than one processor
  SUBROUTINE findIsa

    integer :: zero=0
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: request

    integer, dimension(Np)      :: Psral,NIsal
    integer, dimension(Nl,Np)   :: Isal

    integer             :: iw,l,ll
    integer             :: ierr
    integer             :: Nin
    integer             :: lr,i

    NIsal = 0
    Nsra = 0
    ll = 0

    do iw = 0,Np-1
        if (iw.eq.myid) then
            do l = 0,Np-1
                if (.not.notinlist(l,Psr,Nsr)) then
                    ll = ll + 1
                    call MPI_SEND(NIsf(ll),1,MPI_INTEGER,l,0,MPI_COMM_WORLD,ierr)
                    call MPI_SEND(Isf(1,ll),NIsf(ll),MPI_INTEGER,l,1,MPI_COMM_WORLD,ierr)
                else
                    if(myid .ne. l) then ! prevent sending info to itself
                        call MPI_ISEND(zero,1,MPI_INTEGER,l,0,MPI_COMM_WORLD, request, ierr)
                    end if
                end if
            end do
        else

            call MPI_RECV(Nin,1,MPI_INTEGER,iw,0,MPI_COMM_WORLD,status,ierr)
            if (Nin.ne.0) then
                Nsra = Nsra + 1
                Psral(Nsra) = iw
                NIsal(Nsra) = Nin
                call MPI_RECV(Isal(1,Nsra),Nin,MPI_INTEGER,iw,1,MPI_COMM_WORLD,status,ierr)
            end if
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end do

    if(allocated(Psra)) then
        deallocate (Psra, NIsa, Isa)
    end if
    allocate (Psra(Nsra),NIsa(Nsra))
    allocate(Isa(MAXVAL(NIsal),Nsra))
    Psra = Psral(1:Nsra)
    NIsa = NIsal(1:Nsra)
    Isa = Isal(1:MAXVAL(NIsal),1:Nsra)
    !write(70+myid,"('--------')")
    !write(70+myid,"(4I4)")Psra(1:Nsra)
    !write(70+myid,"(4I4)")NIsa(1:Nsra)
    !do i=1,MAXVAL(NIsal)
    !    write(70+myid,"(4I4)")Isa(i,:)
    !end do
    !write(70+myid,"('--------')")

  END SUBROUTINE findIsa

!findIsf finds interactions pairs/triplets that have extend beyond single processors
  SUBROUTINE findIsf(mnlist3, mnlistcnt, Mbrs, mss, mnlist2, ML2, msslj, mnlistlj, MLlj, mlistcnt)

    integer                           :: Mbrs,ML2,ML3,MLlj, mlistcnt, mnlistcnt
    integer, dimension(mnlistcnt,3)   :: mnlist3
    integer, dimension(Nl,2)          :: mss,msslj
    integer, dimension(mlistcnt)      :: mnlist2
    integer, dimension(MLlj)          :: mnlistlj

    integer, dimension(Np)            :: Psrl
    integer, dimension(0:Np-1)        :: NIsfl
    integer, dimension(Nl,0:Np-1)     :: Isfl
    integer                           :: ii,l,ll,j,i

!initialize
    Nsr = 0
    NIsfl = 0
    Psrl = -1
    Isfl = 0

    !loop through all atoms assigned to processor
    do ii = 1,Nl
        !pair interactions, using refined list mnlist2
        do l = mss(ii,1),mss(ii,2)
            j = mnlist2(l)
            !check if second atom is assigned to processor
            if (P(j).ne.myid) then
                if (notinlist(P(j),Psrl,Np)) then
                    Nsr = Nsr + 1
                    Psrl(Nsr) = P(j)
                end if
                if (notinlist(j,Isfl(1:NIsfl(P(j)),P(j)),NIsfl(P(j)))) then !I have doubt Isfl(1,P(j)) why 1?
                    NIsfl(P(j)) = NIsfl(P(j)) + 1
                    Isfl(NIsfl(P(j)),P(j)) = j
                end if
            end if
        end do
    end do

    !finding forces for the L-J members with inert ions

    do ii = 1,Nl
        do l = msslj(ii,1),msslj(ii,2)
            j = mnlistlj(l)
            if (P(j).ne.myid) then
                if (notinlist(P(j),Psrl,Np)) then
                    Nsr = Nsr + 1
                    Psrl(Nsr) = P(j)
                end if
                if (notinlist(j,Isfl(1:NIsfl(P(j)),P(j)),NIsfl(P(j)))) then !I have doubt Isfl(1,P(j)) why 1?
                    NIsfl(P(j)) = NIsfl(P(j)) + 1
                    Isfl(NIsfl(P(j)),P(j)) = j
                end if
            end if
        end do
    end do

    do l = 1,Mbrs
        do ll = 2,3
            j = mnlist3(l,ll)
            if (P(j).ne.myid) then
                if (notinlist(P(j),Psrl,Np)) then
                    Nsr = Nsr + 1
                    Psrl(Nsr) = P(j)
                end if
                if (notinlist(j,Isfl(1:NIsfl(P(j)),P(j)),NIsfl(P(j)))) then
                    NIsfl(P(j)) = NIsfl(P(j)) + 1
                    Isfl(NIsfl(P(j)),P(j)) = j
                end if
            end if
        end do
    end do

    if (allocated(Psr)) then
        deallocate (Psr, NIsf, Isf)
    end if
    allocate (Psr(Nsr),NIsf(Nsr),Isf(MAXVAL(NIsfl),Nsr))
    ll = 0

    do l = 1,Np
       if (NIsfl(l-1).ne.0) then
          ll = ll + 1
          Psr(ll) = P(Isfl(1,l-1))
          NIsf(ll) = NIsfl(l-1)
          Isf(:,ll) = Isfl(1:MAXVAL(NIsfl),l-1)
       end if
   end do

   !write(60+myid,"('--------')")
   !write(60+myid,"(4I4)")Psr(1:Nsr)
   !write(60+myid,"(4I4)")NIsf(1:Nsr)
   !do i=1,MAXVAL(NIsfl)
   !    write(60+myid,"(4I4)")Isf(i,:)
   !end do
   !write(60+myid,"('--------')")

!    call MPI_BARRIER(MPI_COMM_WORLD,ll); call MPI_FINALIZE(ll); stop
  END SUBROUTINE findIsf

  SUBROUTINE si_nlistL(nlist3, ML3, Nbrs, mnlist3, Mbrs, nss, nlist2, ML2, mss, &
      mnlist2, nlistlj, nsslj, mnlistlj, msslj, MLlj, mnlistcnt, mlistcnt)

    integer                         :: Nbrs,ML3,MLlj, mnlistcnt, mlistcnt
    integer, dimension(ML3,3)       :: nlist3

    integer                         :: Mbrs
    integer, dimension(mnlistcnt,3) :: mnlist3

    integer, dimension(Natm,2)      :: nss, nsslj
    integer                         :: ML2
    integer, dimension(ML2)         :: nlist2
    integer, dimension(MLlj)        :: nlistlj, mnlistlj

    integer, dimension(Nl,2)        :: mss, msslj
    integer, dimension(mlistcnt)    :: mnlist2

    integer :: i,l,ii

    call ntom_list3(nlist3, ML3, Nbrs, mnlist3, Mbrs, mnlistcnt)
    call ntom_list2(nss, nlist2, ML2, mss, mnlist2, nsslj, nlistlj, msslj, mnlistlj, MLlj, mlistcnt)
    call findIsf(mnlist3, mnlistcnt, Mbrs, mss, mnlist2, ML2, msslj, mnlistlj, MLlj, mlistcnt)
    call findIsa
    call collectsends
  END SUBROUTINE si_nlistL

!ntom_list3() extracts all triplets where the first atom is on the assigned processor
  SUBROUTINE ntom_list3(nlist3,ML3,Nbrs,mnlist3,Mbrs,mnlistcnt)

    integer                         :: Nbrs, ML3, mnlistcnt
    integer, dimension(ML3,3)       :: nlist3

    integer                         :: Mbrs
    integer, dimension(mnlistcnt,3) :: mnlist3

    integer    :: i, l

    !intialize counters
    mnlist3 = 0
    Mbrs = 0

!loop through the local neighborlists and check if first atom is assigned to current
!processor
    do l = 1,Nbrs
        if (P(nlist3(l,1)).eq.myid) then
            Mbrs = Mbrs + 1
            if(Mbrs .ge. mnlistcnt) print*, "********** Error, writing to unallocated space inside ntom_list3 ********", myid, Mbrs, mnlistcnt
            mnlist3(Mbrs,:) = nlist3(l,:)
       end if
   end do

  END SUBROUTINE ntom_list3

!ntom_list2() extracts all pair interactions where at least the first atom is on the assigned
!processor
  SUBROUTINE ntom_list2(nss, nlist2, ML2, mss, mnlist2, nsslj, nlistlj, msslj, mnlistlj, MLlj, mlistcnt)

    integer, dimension(Natm,2)        :: nss, nsslj
    integer                           :: ML2, MLlj, mlistcnt
    integer, dimension(ML2)           :: nlist2
    integer, dimension(MLlj)          :: nlistlj, mnlistlj

    integer, dimension(Nl,2)          :: mss, msslj
    integer, dimension(mlistcnt)      :: mnlist2

    integer                           :: i, ii

!initialize counters, pointers
    ii = 0
    mss(:,1) = 0
    mss(:,2) = -1
    mnlist2 = 0

    do i = 1,Natm
        if (P(i).eq.myid) then !added by kallol, definition mismatch of ii
            ii = ii + 1
            if (ii.eq.1) then
                mss(ii,1) = 1
                msslj(ii,1) = 1
            else
                mss(ii,1) = mss(ii-1,2) + 1
                msslj(ii,1) = msslj(ii-1,2) + 1
            end if

            mss(ii,2) = mss(ii,1) + nss(i,2) - nss(i,1)
            mnlist2(mss(ii,1):mss(ii,2)) = nlist2(nss(i,1):nss(i,2))

            msslj(ii,2) = msslj(ii,1) + nsslj(i,2) - nsslj(i,1)
            mnlistlj(msslj(ii,1):msslj(ii,2)) = nlistlj(nsslj(i,1):nsslj(i,2))

            if(mss(ii,2).ge. mlistcnt) print*, "***************** An error has occured inside ntom_list2*****************", myid, mss(ii,1), mss(ii,2), mlistcnt
        end if
    end do
    !print*,myid, mss(ii,1), mss(ii,2), mlistcnt

  END SUBROUTINE ntom_list2

END MODULE par
