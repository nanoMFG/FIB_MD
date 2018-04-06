!------------------------------------------------------------------------------
!
! `Si` Source File
!
! si.f90 source function file. This file contains subroutines for the actual
! generation of neighborlists. It also contains the subroutines that use these
! neighborlists to calculate interatomic forces and potential energies using
! stillinger-weber and molliere functions.
!
!------------------------------------------------------------------------------

MODULE stillweb

  USE prms
  USE par
  USE nlistmod

  IMPLICIT none

  integer                        :: ML3, ML2, MLlj, mnlistcnt, mlistcnt
  integer,allocatable, dimension(:,:)  :: nlist3 !3 body neighborlist
  integer,allocatable, dimension(:,:)  :: mnlist3 !3 body neighborlist that each processor uses

  integer                       :: Nbrs, Mbrs !Nbrs number of 3body neighbors

  integer,allocatable,dimension(:)   :: nlist2, nlistlj !2 body neighborlists, Si-Si and Si-Ga
  integer,allocatable,dimension(:,:) :: nss, nsslj !index pointers for 2 body lists

  integer,allocatable,dimension(:)   :: mnlist2, mnlistlj
  integer,allocatable,dimension(:,:) :: mss, msslj

!  Lookup tables
  real                              :: umax                ! max table parameter
  real,allocatable,dimension(:)     :: table2bd            ! table for 2-body potential
  real,allocatable,dimension(:)     :: u                   ! array of table absissas


  integer,dimension(3) :: Nc                ! number of cells for neighborlist (**computed**)
  real                 :: rc                ! cuttoff for neighborlist
  real                 :: rlj2              ! cuttoff for L-J neighborlist
  real                 :: rclj              ! nlist2 cuttoff for L-J neighborlist

CONTAINS

!init_nlist initializes the neighborlist for each processor,
!needs to be called every time atoms are reallocated
!Nl is # atoms on each processor, Nlj is # ions fired
  SUBROUTINE init_nlist

    ML3 = Nl*200   !Natm*150
    MLlj = Nlj*1500  !Nlj*10000 ! changed value kallol
    ML2 = Nl*40        !Natm*50 ! kallol prev 20
    mnlistcnt = 150*Nl
    mlistcnt = 15*Nl

    rc = listfac*al*MAX(sigmaGe,sigmaSi) !SW cutoff radius plus 10% buffer
    Nc = MAX(1,INT(Lb/rc))               !cell size to split atoms into

    rlj2 = 2.5*MAX(sigmaGe,sigmaSi)      !cutoff radius for moliere
    rclj = listfac*rlj2                  !moliere cutoff with 10% buffer
    rlj2 = rlj2**2                       !squared

    if (allocated(mss)) then
        deallocate(nlist3, nlist2, nlistlj, mss, msslj, mnlist3, mnlist2, mnlistlj)
        allocate (nlist3(ML3,3),nlist2(ML2))
        allocate (nlistlj(MLlj)) !check the number, also check with par.f90
        allocate (mss(Nl,2), msslj(Nl,2))
        allocate (mnlist3(mnlistcnt,3),mnlist2(mlistcnt))
        allocate (mnlistlj(MLlj)) !check the hard coded number
    else
        allocate (nss(Natm,2))
        allocate (nsslj(Natm,2))
        allocate (nlist3(ML3,3),nlist2(ML2))
        allocate (nlistlj(MLlj)) !check the number, also check with par.f90
        allocate (mss(Nl,2), msslj(Nl,2))
        allocate (mnlist3(mnlistcnt,3),mnlist2(mlistcnt))
        allocate (mnlistlj(MLlj)) !check the hard coded number
    end if

  END SUBROUTINE init_nlist


!si_nlist(X) is the top level call for generating the neighborlists for silicon
!interactions. First calls the global neighborlist generator, then the per-processor
!nlist2 generator

  SUBROUTINE si_nlist(X)
    real   ,dimension(Natm,3)                               :: X   ! particle positions
    integer :: ierr !,lt

    call si_nlistG(X)
!    print*, "si_nlistG done"
    call si_nlistL(nlist3, ML3, Nbrs, mnlist3, Mbrs, nss, nlist2, ML2, mss, &
          mnlist2, nlistlj, nsslj, mnlistlj, msslj, MLlj, mnlistcnt, mlistcnt)

    !call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! debugging, take it off in production runs

  END SUBROUTINE si_nlist

!si_nlistG(X) generates the neighborlist for each processor, both Si-Si and Si-Ga
!uses current positions X
  SUBROUTINE si_nlistG(X)
    real   ,dimension(Natm,3)                               :: X   ! particle positions
    integer,dimension(Natm)                                 :: LL  ! next in nlist2 pointer
    integer,dimension(0:Nc(1)+1,0:Nc(2)+1,0:Nc(3)+1)        :: HOC ! Head-of-Chain pointer
                                                                ! with perioidic continuation
    real,dimension(3)     :: xx
    integer               :: ic1,ic2,ic3 !,lt
    real                  :: rij2,rv2, rv2lj
    integer               :: i,j,k,l,m,n,n1,n2,n3
    integer               :: nbr_cnt, nbr_cntlj
    logical               :: flagi, flagj

    !parameters used
    rv2 = rc*rc !cutoff radius squared, for silicon SW
    rv2lj = rclj*rclj !cutoff squared for ion interactions

    !generate LL and HOC, which describe how atoms are split apart into NC cells
    call chainlist(LL,HOC,Nc,X)

    nbr_cnt = 1
    nss(:,2) = -1
    nss(:,1) = 0
    nbr_cntlj = 1
    nsslj(:,2) = -1
    nsslj(:,1) = 0
    flagi = .true.
    flagj = .true.

    do ic3 = 1,Nc(3)
        do ic2 = 1,Nc(2)
            do ic1 = 1,Nc(1)
                i = HOC(ic1,ic2,ic3)         ! start with atom at the head of the linked nlist2
                do while(i.ne.0)             ! loop over all atoms in the linked nlist2 if any

                    flagi = .true.
                      !P_big(i) is the processor assigned to atom i; check whether the processor, myid,
                      !locally assembling this nlist2, is associated with one of the atoms linked together
                    if(P_big(i).ne. myid) flagi = .false.

                    nss(i,1) = nbr_cnt
                    nsslj(i,1) = nbr_cntlj

                    !loop through atoms in the neighboring cells to find neighbors
                    do n1 = -1,1
                        do n2 = -1,1
                            do n3 = -1,1
                                !take atom j as the head in cell atom for cell ic1+n1,ic2+nc2,ic3+nc3
                                j = HOC(ic1+n1,ic2+n2,ic3+n3)
                                do while(j.ne.0)

                                    flagj = .true.
                                    if (P_big(j).ne. myid) flagj = .false. !again, check which processor is assigned

                                    if ( (flagi .or. flagj) .and. (i.ne.j) ) then
                                        xx = X(i,:) - X(j,:) - NINT( (X(i,:) - X(j,:))*iLb)*Lb
                                        rij2 = xx(1)*xx(1)+xx(2)*xx(2)+xx(3)*xx(3)
                                        if (MAX(atype(i),atype(j)).lt.3) then !if not ions, record for `nlist2`
                                            if ( (rij2 .lt. rv2) ) then !check against nlist3 cutoff radius
                                                nlist2(nbr_cnt) = j       !record atom 2's index
                                                nbr_cnt = nbr_cnt + 1              ! increment number of neighbors found
                                            end if
                                        else !if an ion is involved with pair, use separate `nlistlj`
                                            if ( (rij2 .lt. rv2lj) ) then
                                                nlistlj(nbr_cntlj) = j
                                                nbr_cntlj = nbr_cntlj+1
                                                if (nbr_cntlj .eq. MLlj) then  !check the previously hardcoded value
                                                    print *, "******nbr_Cntlj exceeded****** (too many ion neighbors)"
                                                end if
                                            end if
                                        end if
                                    end if

                                    j = LL(j) !iterate

                                end do
                            end do
                        end do
                    end do

                    nss(i,2) = nbr_cnt-1 !last nlist2 index for i's neighbors
                    nsslj(i,2) = nbr_cntlj-1 !last nlist2 index for ion i's neighbors
                    i = LL(i)

                end do

            end do
        end do
    end do

!  This bit finds all threesomes.
!
!   RULES:  -2nd and third atoms have index > first
!           -first get 3-somes from present i-row
!           -look at lower (higher i) rows only for determining second neighbor
!           -check present i-row to make sure no repeat
    Nbrs = 0 !number of triplets found

    do i = 1,Natm
        do l = nss(i,1),nss(i,2)
            j = nlist2(l)
            if (j.gt.i) then
                ! get threesomes from present i-row
                do m = l+1,nss(i,2)
                    k = nlist2(m)
                    if (k.gt.i) then
                        ! add threesome
                        Nbrs = Nbrs + 1
                        nlist3(Nbrs,1) = MIN(i,j,k)
                        nlist3(Nbrs,2) = j+k+i - MIN(i,j,k) - MAX(i,j,k)
                        nlist3(Nbrs,3) = MAX(i,j,k)
                    end if
                end do
                ! get threesomes from later (only larger i) rows
                do m = nss(j,1),nss(j,2)
                    k = nlist2(m)
                    if (k.gt.i) then
                        ! make sure not repeating 3-some from present i-row
                        do n = nss(i,1),nss(i,2)
                            if (k.eq.nlist2(n)) goto 11
                        end do
                        ! add threesome
                        Nbrs = Nbrs + 1
                        nlist3(Nbrs,1) = MIN(i,j,k)
                        nlist3(Nbrs,2) = j+k+i - MIN(i,j,k) - MAX(i,j,k)
                        nlist3(Nbrs,3) = MAX(i,j,k)
                    end if
11                  continue
                end do
            end if
        end do
    end do

! remove j < i neighbors
    do i = 1,Natm
       do l = nss(i,1),nss(i,2)
          j = nlist2(l)
          do while (j.lt.i.and.nss(i,2).gt.l)
             nlist2(l) = nlist2(nss(i,2))
             j = nlist2(l)
             nss(i,2) = nss(i,2) - 1
          end do
          if (nss(i,2).eq.l.and.j.lt.i) nss(i,2) = nss(i,2) - 1

       end do
   end do

   ! remove j < i neighbors
   do i = 1,Natm
       do l = nsslj(i,1),nsslj(i,2)
           j = nlistlj(l)
           do while (j.lt.i.and.nsslj(i,2).gt.l)
               nlistlj(l) = nlistlj(nsslj(i,2))
               j = nlistlj(l)
               nsslj(i,2) = nsslj(i,2) - 1
           end do
           if (nsslj(i,2).eq.l.and.j.lt.i) nsslj(i,2) = nsslj(i,2) - 1

       end do
   end do

   !write(*,"(7I10)")myid, Nbrs, ML3, nbr_cnt, ML2, nbr_cntlj, MLlj

   if(Nbrs .ge. ML3) then
       write(*,"('Nbrs (triplets) exceeded ML3 at ', I5, ' th processor, Nbrs = ', I11, 'ML3 = ', I11)")myid, Nbrs, ML3
   end if
   if(nbr_cnt .ge. ML2) then
       write(*,"('nbr_cnt (pairs) exceeded ML2 at ', I5, ' th processor, nbr_cnt = ', I11, 'ML2 = ', I11)")myid, nbr_cnt, ML2
   end if
   if(nbr_cntlj .ge. MLlj) then
       write(*,"('nbr_cntlj (ion involved pairs) exceeded MLlj at ', I5, ' th processor, nbr_cntlj = ', I11, 'MLlj = ', I11)")myid, nbr_cntlj, MLlj
   end if

   !print*,myid, Nbrs, nbr_cnt, nbr_cntlj, Nl

  END SUBROUTINE si_nlistG

  SUBROUTINE si_potential(X,Utot,u)

    real, dimension(Natm,3) :: X
    real, dimension(Natm)   :: u
    real                    :: Utot

    real                    :: rij,rik,rjk, rij2
    real                    :: rijs,riks,rjks
    real                    :: r

    real                    :: ctjik,ctijk,ctikj
    real, dimension(3)      :: xxij,xxik,xxjk

    integer                 :: i,j,k,l,ii

    real                    :: ul

    Utot = 0.
    u = 0.
    do ii = 1,Nl
        !Stillinger-Weber 2 body potentital
        i = il(ii)
        do l = mss(ii,1),mss(ii,2)
            j = mnlist2(l)
            xxij = X(i,:) - X(j,:) - NINT( (X(i,:) - X(j,:))*iLb)*Lb
            rij = SQRT( xxij(1)**2 + xxij(2)**2 + xxij(3)**2)
            r = rij*isigma(atype(i),atype(j))
            if (r .lt. 0.999*al) then
                ul =  0.5*eps(atype(i),atype(j))*A*(B*r**(-psi) - r**(-qsi))*EXP(1./(r-al))
                u(i) = u(i) + ul
                u(j) = u(j) + ul
                Utot = Utot + 2.*ul
            end if
        end do

        !Moliere potential
        do l = msslj(ii,1),msslj(ii,2)
            j = mnlistlj(l)
            xxij = X(i,:) - X(j,:) - NINT( (X(i,:) - X(j,:))*iLb)*Lb
            rij2 = SUM(xxij*xxij)
            if (rij2.lt.rlj2) then
                rij = SQRT(rij2)
                ul = ie0*Zmol(atype(i))*Zmol(atype(j))*ec*ec/rij*(  &
                     Amol*EXP(-alpha_mol*rij/aamol(atype(i),atype(j))) &
                     + Bmol*EXP(-beta_mol*rij/aamol(atype(i),atype(j))) &
                     + Cmol*EXP(-gamma_mol*rij/aamol(atype(i),atype(j))) )
                u(i) = u(i) + ul/2.
                u(j) = u(j) + ul/2.
                Utot = Utot + ul
            end if
        end do


    end do

    !Stillinger-Weber 3 body potential
    do l = 1,Mbrs

       i = mnlist3(l,1)
       j = mnlist3(l,2)
       k = mnlist3(l,3)

       xxij = X(i,:) - X(j,:) - NINT( (X(i,:) - X(j,:))*iLb)*Lb
       rij = SQRT( xxij(1)**2 + xxij(2)**2 + xxij(3)**2)
       rijs = rij*isigma(atype(i),atype(j))

       xxik = X(i,:) - X(k,:) - NINT( (X(i,:) - X(k,:))*iLb)*Lb
       rik = SQRT( xxik(1)**2 + xxik(2)**2 + xxik(3)**2)
       riks = rik*isigma(atype(i),atype(k))

       xxjk = X(j,:) - X(k,:) - NINT( (X(j,:) - X(k,:))*iLb)*Lb
       rjk = SQRT( xxjk(1)*xxjk(1) + xxjk(2)*xxjk(2) + xxjk(3)*xxjk(3))
       rjks = rjk*isigma(atype(j),atype(k))

       ctjik = (xxij(1)*xxik(1) + xxij(2)*xxik(2) + xxij(3)*xxik(3))/rij/rik
       ctijk = -(xxij(1)*xxjk(1) + xxij(2)*xxjk(2) + xxij(3)*xxjk(3))/rij/rjk
       ctikj = (xxik(1)*xxjk(1) + xxik(2)*xxjk(2) + xxik(3)*xxjk(3))/rik/rjk

       ul = 1./3.*leps(atype(i),atype(j),atype(k))  &
            *(h(rijs,riks,ctjik) + h(rijs,rjks,ctijk) + h(riks,rjks,ctikj))
       u(i) = u(i) + ul
       u(j) = u(j) + ul
       u(k) = u(k) + ul
       Utot = Utot + 3.*ul

    end do

  !  Utot = SUM(u)!
    !print *,myid,": Utot = ",Utot

  END SUBROUTINE si_potential

!function definition for part of S-W potential
  FUNCTION h(u,v,w)
    real        :: h,u,v,w

    if (u.lt.0.999*al.and.v.lt.0.999*al) then
       h = EXP(gamma/(u-al) + gamma/(v-al))*(w + 1./3.)**2
    else
       h = 0.
    end if

  END FUNCTION h


! the Stillinger-Weber force
  SUBROUTINE si_force(Xp,F,lt)
    real, dimension(Natm,3) :: Xp
    real, dimension(Natm,3) :: F

    real                    :: rij,rik,rjk, rij2
    real                    :: rijs,riks,rjks
    real                    :: df, fc
    real, dimension(3)      :: fi,fj,fk
    real, dimension(3)      :: xxij,xxik,xxjk
    real, dimension(3)      :: dxxij,dxxik,dxxjk

    real                    :: ctjik,ctijk,ctikj
    real, dimension(3)      :: dctjik,dctijk,dctikj

    real                    :: aa
    real                    :: fjik,fijk,fikj
    real                    :: Erijs,Eriks,Erjks
    real                    :: dErijs,dEriks,dErjks
    real                    :: Wjik,Wijk,Wikj

    real                    :: irijsal,iriksal,irjksal

    real                    :: isigij,isigik,isigjk,lepsijk
    integer                 :: i,j,k,l,m,ii,lt
    integer                 :: ierr

    F = 0.
    do ii = 1,Nl
        i = il(ii)
        do l = mss(ii,1),mss(ii,2)
            j = mnlist2(l)
            !if(sputter_index_l(i).gt.0 .or. sputter_index_l(j).gt.0 .or. i.gt.Nsg+impact .or. j.gt.Nsg+impact) cycle !kallol
            if (i.gt.Nsg+impact .or. j.gt.Nsg+impact) cycle
            xxij = Xp(i,:) - Xp(j,:) - NINT( (Xp(i,:) - Xp(j,:))*iLb)*Lb
            rij = SQRT( xxij(1)**2 + xxij(2)**2 + xxij(3)**2)
            rijs = rij*isigma(atype(i),atype(j))
            df = itable2bd(rijs,i,j,lt)

            fi = -df*xxij/rij*epssig(atype(i),atype(j))
            F(i,:) = F(i,:) + fi
            F(j,:) = F(j,:) - fi
        end do

        do l = msslj(ii,1),msslj(ii,2)
            j = mnlistlj(l)
            !if(sputter_index_l(i).gt.0 .or. sputter_index_l(j).gt.0 .or. i.gt.Nsg+impact .or. j.gt.Nsg+impact) cycle !kallol
            if (i.gt.Nsg+impact .or. j.gt.Nsg+impact) cycle
            xxij = Xp(i,:) - Xp(j,:) - NINT( (Xp(i,:) - Xp(j,:))*iLb)*Lb
            rij2 = SUM(xxij*xxij)
            if (rij2.lt.rlj2) then
                rij = SQRT(rij2)
                aa = aamol(atype(i),atype(j))
                fc = ie0*Zmol(atype(i))*Zmol(atype(j))*ec*ec/rij*(  ( &
                     Amol*EXP(-alpha_mol*rij/aa) &
                     + Bmol*EXP(-beta_mol*rij/aa) &
                     + Cmol*EXP(-gamma_mol*rij/aa) )/rij   &
                     + ( &
                     Amol*alpha_mol/aa*EXP(-alpha_mol*rij/aa) &
                     + Bmol*beta_mol/aa*EXP(-beta_mol*rij/aa) &
                     + Cmol*gamma_mol/aa*EXP(-gamma_mol*rij/aa) ))
                fi = fc*xxij/rij
                F(i,:) = F(i,:) + fi
                F(j,:) = F(j,:) - fi
            end if
        end do


    end do

    do l = 1,Mbrs

       i = mnlist3(l,1)
       j = mnlist3(l,2)
       k = mnlist3(l,3)

       !if(sputter_index_l(i).gt.0 .or. sputter_index_l(j).gt.0 .or. sputter_index_l(k).gt.0 .or. i.gt.Nsg+impact .or. j.gt.Nsg+impact .or. k.gt.Nsg+impact) cycle !kallol
       if(i.gt.Nsg+impact .or. j.gt.Nsg+impact .or. k.gt.Nsg+impact) cycle
       isigij = isigma(atype(i),atype(j))
       isigik = isigma(atype(i),atype(k))
       isigjk = isigma(atype(j),atype(k))
       lepsijk = leps(atype(i),atype(j),atype(k))

       xxij = Xp(i,:) - Xp(j,:) - NINT( (Xp(i,:) - Xp(j,:))*iLb)*Lb
       rij = 1./SQRT( xxij(1)*xxij(1) + xxij(2)*xxij(2) + xxij(3)*xxij(3))
       dxxij = xxij*rij
       rijs = isigij/rij

       xxik = Xp(i,:) - Xp(k,:) - NINT( (Xp(i,:) - Xp(k,:))*iLb)*Lb
       rik = 1./SQRT( xxik(1)*xxik(1) + xxik(2)*xxik(2) + xxik(3)*xxik(3))
       dxxik = xxik*rik
       riks = isigik/rik

       xxjk = Xp(j,:) - Xp(k,:) - NINT( (Xp(j,:) - Xp(k,:))*iLb)*Lb
       rjk = 1./SQRT( xxjk(1)*xxjk(1) + xxjk(2)*xxjk(2) + xxjk(3)*xxjk(3))
       dxxjk = xxjk*rjk
       rjks = isigjk/rjk

       ctjik = (dxxij(1)*dxxik(1) + dxxij(2)*dxxik(2) + dxxij(3)*dxxik(3))
       dctjik = dxxij*rik + dxxik*rij - ctjik*(dxxij*rij + dxxik*rik)
       ctijk = -(dxxij(1)*dxxjk(1) + dxxij(2)*dxxjk(2) + dxxij(3)*dxxjk(3))
       dctijk = - dxxjk*rij - ctijk*dxxij*rij
       ctikj = (dxxik(1)*dxxjk(1) + dxxik(2)*dxxjk(2) + dxxik(3)*dxxjk(3))
       dctikj = + dxxjk*rik - ctikj*dxxik*rik

       Erijs = 0.; Eriks = 0.; Erjks = 0.

       irijsal = 1./(rijs-al)
       iriksal = 1./(riks-al)
       irjksal = 1./(rjks-al)

       if (rijs.lt.0.999*al) Erijs = EXP(gamma*irijsal)
       if (riks.lt.0.999*al) Eriks = EXP(gamma*iriksal)
       if (rjks.lt.0.999*al) Erjks = EXP(gamma*irjksal)

       if (rijs.lt.0.999*al.and.rjks.lt.0.999*al) then
          fijk = 1.
       else
          fijk = 0.
       end if
       if (riks.lt.0.999*al.and.rjks.lt.0.999*al) then
          fikj = 1.
       else
          fikj = 0.
       end if
       if (rijs.lt.0.999*al.and.riks.lt.0.999*al) then
          fjik = 1.
       else
          fjik = 0.
       end if

       dErijs =  -gamma*irijsal*irijsal
       dEriks =  -gamma*iriksal*iriksal
       dErjks =  -gamma*irjksal*irjksal
       Wjik = (ctjik + thrd)*(ctjik + thrd)
       Wijk = (ctijk + thrd)*(ctijk + thrd)
       Wikj = (ctikj + thrd)*(ctikj + thrd)


       fi = -lepsijk*(    &
              Erijs*(  fjik*Eriks*( Wjik*(dErijs*dxxij*isigij  &
                                   + dEriks*dxxik*isigik )  &
                             + (ctjik+thrd)*dctjik*2. )   &
                     + fijk*Erjks*(dErijs*Wijk*dxxij*isigij  &
                             + (ctijk+thrd)*dctijk*2. ) )   &
            + fikj*Eriks*Erjks*( dEriks*Wikj*dxxik*isigik  &
                     + (ctikj+thrd)*dctikj*2.  ) )


       dctjik = - dxxik*rij + ctjik*dxxij*rij
       dctijk = - dxxij*rjk + dxxjk*rij + ctijk*(dxxij*rij-dxxjk*rjk)
       dctikj = + dxxik*rjk - ctikj*dxxjk*rjk

       fj = -lepsijk*(    &
            + Erijs*(fjik*Eriks*( -dErijs*Wjik*dxxij*isigij  &
                          + (ctjik+thrd)*dctjik*2. )   &
                   + fijk*Erjks*(Wijk*(-dErijs*dxxij*isigij &
                                 + dErjks*dxxjk*isigjk ) &
                           + (ctijk+thrd)*dctijk*2.  ) )  &
            + fikj*Eriks*Erjks*( dErjks*Wikj*dxxjk*isigjk  &
                          + (ctikj+thrd)*dctikj*2. ) )

       fk = -fi - fj

       F(i,:) = F(i,:) + fi
       F(j,:) = F(j,:) + fj
       F(k,:) = F(k,:) + fk


   end do
   call combfN(F, lt)

  END SUBROUTINE si_force

  FUNCTION itable2bd(ue,ii,jj,lt)
    real      ::  itable2bd
    real      ::  ue,rij
    integer   ::  i,ii,jj,lt,j
    real      ::  d1,d2,d3
    real, dimension(3)::xxij

    i = INT(ue/umax*REAL(Ntable-1)) + 1
    if (i .gt. Ntable) then
        write(*,"('ERRRRRRROR: proc',I3,' atom',I8,' &',I8, ' tstep = ',I9,' z =',2E11.3)")myid, ii,jj,lt,X(ii,3),X(jj,3)
        itable2bd = 0.0
    else
        d1 = (ue-u(i))
        d2 = (ue-u(i+1))
        d3 = (ue-u(i+2))

        itable2bd = (0.5*d2*d3*table2bd(i) - d1*d3*table2bd(i+1) + 0.5*d1*d2*table2bd(i+2))
    end if

  END FUNCTION itable2bd

  SUBROUTINE mktables

    integer   ::  i
    real      ::  du

    allocate(table2bd(Ntable))
    allocate(u(Ntable))

    umax = 2.5*listfac*al
    du = umax/REAL(Ntable-1)
    do i = 2,Ntable
       u(i) = umax*REAL(i-1)/REAL(Ntable-1)
       du = umax/REAL(Ntable-1)
       if (u(i).lt..999*al) then
          table2bd(i) = -A*EXP(1./(u(i)-al))*( &
                  (B*u(i)**(-psi) - u(i)**(-qsi))/(u(i)-al)**2 &
                  + (psi*B*u(i)**(-psi-1) - qsi*u(i)**(-qsi-1))   ) / du**2
       else
          table2bd(i) = 0.
      end if
    end do

  END SUBROUTINE mktables

END MODULE stillweb
