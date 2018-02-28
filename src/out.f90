!------------------------------------------------------------------------------
!
! `OUT` Source File
!
! out.f90 source function file. This contains all of the subroutines used for
! initializing data outputs and then actually writing output data when called
! during the simulation.
!
!------------------------------------------------------------------------------

MODULE out

  USE prms
  USE stats
  USE temp
  USE data
  USE parallel
  USE par


  IMPLICIT none

  integer                       :: NTbar
  real,allocatable,dimension(:) :: TbarS

  integer                         :: NYen
  real,allocatable,dimension(:,:) :: YenS

  integer                         :: NmaxS
  real,allocatable,dimension(:,:) :: MaxwS

CONTAINS

  SUBROUTINE writerestart(X,V,lt)
    real, dimension(Natm,3) :: X,V
    integer                 :: lt
    real                    :: time
    character(30)           :: fn

    call updateroot(X)
    call updateroot(V)

    if (myid.eq.0) then
        write(fn,"('D/restart.in_old')")
        open(res_unit,file=fn,form='UNFORMATTED')
        write(res_unit)lt,restart_time
        write(res_unit)X,V,mass,atype
        write(res_unit)impact
        close(res_unit)

        write(fn,"('D/restart.in',I6.6,'_',I9.9)")impact,lt
        open(res_unit,file=fn,form='UNFORMATTED')
        write(res_unit)lt,restart_time
        write(res_unit)X,V,mass,atype
        write(res_unit)impact
        close(res_unit)
    end if

  END SUBROUTINE writerestart

  SUBROUTINE writeatoms(X,V,F,lt,time,temp)
    real, dimension(Natm,3) :: X,V,F
    real, dimension(Natm) :: kin_eng
    integer  :: lt, fon
    integer  :: i,k,l,j,ii
    integer :: ierr
    real    :: keprint, temp, TE, time
    character(30)           :: fn
    character(2)            :: tempatype
    !call updateroot(X)
    !call updateroot(F)
    !call updateroot(V)

    if (myid.eq.0) then		!added by Josh, converted output to xyz format
!        write(fn,"('data/mdrun2.xyz')") !moved to general initio / closeio
!        open(atm_unit,file=fn,action='write',position='append')
        write(atm_unit,"(I9)"),Natm !heads the section with # of atoms
        write(atm_unit,"(1A12,F12.5,3A12,2A8,3A16,A8,F14.5)")'Time (ps) = ',time*1e12,'x (A)','y (A)','z (A)','type','proc#','Vx','Vy','Vz','atom#',temp
        do i = 1,Natm
            kin_eng(i) = 0.5*mass(i)*SUM(V(i,:)**2)/ec
      	    if(atype(i).eq.1) then
			          tempatype = "Si"
            elseif(atype(i).eq.3) then
		            tempatype = "Ga"
            end if
    	      !write(atm_unit,"(3E20.10E3,A5,I3,7E20.10E3, I8)")X(i,:)*1e10,tempatype,P(i),kin_eng(i),V(i,:),F(i,:),i
        	  write(atm_unit,"(3F12.5,A8,I8,3E14.5, I8)")X(i,:)*1e10,tempatype,P(i),V(i,:),i
        end do
!        close(atm_unit)
        !write a line tracking the current ion's trajectory and system temperature
    		write(*,"(I10,6E20.10)")lt,temp,V(Nsg+impact,3)
            !write(*,"(I10,4E20.10,I10)")lt,X(Nsg+impact,3),V(Nsg+impact,3), MAXVAL(kin_eng), X(MAXLOC(kin_eng),3), MAXLOC(kin_eng)
        end if

    if (myid.eq.-1) then
        write(fn,"('G/atoms',I9.9,'.dat')")lt
        open(atm_unit,file=fn)
        do i = 1,Natm
            kin_eng(i) = 0.5*mass(i)*SUM(V(i,:)**2)/ec
            write(atm_unit,"(3E20.10E3,2I3,7E20.10E3, I8)")X(i,:),atype(i),P(i),kin_eng(i),V(i,:),F(i,:),i
        end do
        close(atm_unit)
        write(*,"(I10,4E20.10,I10)")lt,X(Nsg+impact,3),V(Nsg+impact,3), MAXVAL(kin_eng), X(MAXLOC(kin_eng),3), MAXLOC(kin_eng)
    end if

    if (myid .eq. -1) then
        !added by Kallol
        if(lt .gt. 0) then
            write(fn, "('D/proc_',I2.2,'_',I9.9,'.dat')"),myid,lt
            open(atm_unit, file=fn)
            do ii = 1,Nl
                i = il(ii)
                write(atm_unit,"(3E20.10,3I3,7E20.10)")X(i,:),atype(i),P(i),i,0.5*mass(i)*SUM(V(i,:)**2),V(i,:),F(i,:)
                do l = mss(ii,1),mss(ii,2)
                    j = mlist(l)
                    if(P(j).eq.myid) then
                        !write(atm_unit,"(3E20.10,2I3)")X(j,:),atype(j),P(j)
                    else
                        write(atm_unit,"(3E20.10,3I3,7E20.10)")X(j,:),atype(j),20,j,0.5*mass(j)*SUM(V(j,:)**2),V(j,:),F(j,:)
                    end if
                end do
            end do
            do l = 1,Mbrs
                i = mnlist(l,1)
                j = mnlist(l,2)
                k = mnlist(l,3)
                if(P(i).ne.myid) then
                    write(atm_unit,"(3E20.10,3I3,7E20.10)")X(i,:),atype(i),30,i,0.5*mass(i)*SUM(V(i,:)**2),V(i,:),F(i,:)
                end if
                if(P(j).ne.myid) then
                    write(atm_unit,"(3E20.10,3I3,7E20.10)")X(j,:),atype(j),30,j,0.5*mass(j)*SUM(V(j,:)**2),V(j,:),F(j,:)
                end if
                if(P(k).ne.myid) then
                    write(atm_unit,"(3E20.10,3I3,7E20.10)")X(k,:),atype(k),30,k,0.5*mass(k)*SUM(V(k,:)**2),V(k,:),F(k,:)
                end if
            end do
        end if
        close(atm_unit)

    end if

  END SUBROUTINE writeatoms

  SUBROUTINE writeabn(X,V,lt)
    real, dimension(Natm,3) :: X,V
    integer  :: lt
    integer  :: i,k

    character(30)           :: fn

    if (myid.eq.0) then
       write(fn,"('D/abn.out',I9.9)")lt
       open(atm_unit,file=fn,form='UNFORMATTED')
       write(atm_unit)X,V
       close(atm_unit)
    end if

  END SUBROUTINE writeabn

  SUBROUTINE initio

    character(10)  fn

    !if (myid.eq.0) write(out_unit,*) "INITIALIZING I/O"

    write(fn,"(I10.10)")Nt0

    if (myid.eq.0) then
       if (zen_out.gt.0) open(zen_unit,file='G/zen.out'//fn)
       if (tmp_out.gt.0) open(tmp_unit,file='G/tmp.out'//fn)
       if (eng_out.gt.0) open(eng_unit,file='G/eng.out'//fn)
       if (atm_out.gt.0) open(atm_unit,file='data/mdrun2.xyz')
    end if

  END SUBROUTINE initio

  SUBROUTINE closeio

    if (myid.eq.0) then
       if (tmp_out.gt.0) close(tmp_unit)
       if (eng_out.gt.0) close(eng_unit)
       if (zen_out.gt.0) close(zen_unit)
       if (atm_out.gt.0) close(atm_unit)
    end if

  END SUBROUTINE closeio

  SUBROUTINE iostuff(X,V,F,lt,time,temp)
    real, dimension(Natm,3)  :: X,V,F
    real, dimension(Natm)    :: q
    real                     :: time,temp
    integer                  :: lt

    if (res_out.gt.0.and.MOD(lt,res_out).eq.0.or.   &
         abn_out.gt.0.and.MOD(lt,abn_out).eq.0.or.   &
         atm_out.gt.0.and.MOD(lt,atm_out).eq.0)  then
          call updateroot(X)
          call updateroot(V)
          call temperature(V,temp)
          if (res_out.gt.0.and.MOD(lt,res_out).eq.0) then
              call writerestart(X,V,lt)
          end if
          if (abn_out.gt.0.and.MOD(lt,abn_out).eq.0) then
!              call updateglobal(V)
              call writeabn(X,V,lt)
          end if
          if (atm_out.gt.0.and.MOD(lt,atm_out).eq.0) then
              call writeatoms(X,V,F,lt,time,temp)
          end if
    end if
    call diagnostics(X,V,lt,time)

  END SUBROUTINE iostuff

  SUBROUTINE initmeans

    !if (myid.eq.0) write(out_unit,*) "INIT MEANS"

    allocate (TbarS(Nlat(1)))
    NTbar = 0
    TbarS = 0.

  END SUBROUTINE initmeans

  SUBROUTINE diagnostics(X,V,lt,time)
    real, dimension(Natm,3)  :: X,V
    real, dimension(Natm)    :: q
    real                     :: time
    integer                     :: lt

    real,dimension(Nlat(1))      :: Tbar
    real,dimension(Nlat(1),5)    :: Yen

    real                      :: T,KE,PE           ! Temperatures + Energies
    real,dimension(3)         :: Pm                ! Momentum
    real,dimension(0:Nmaxw,3) :: Maxw              ! Maxwellian distributions
    integer,dimension(Natm)   :: Lall              ! list of all atoms
    real                      :: e1,e2
    real,dimension(0:Nmaxg)   :: g, rg             ! Radial distribution fct
    character(30)             :: fn
    integer                   :: ib,i,j

    if (tmp_out.gt.0.and.MOD(lt,tmp_out).eq.0) then
       call temperature(V,T)
       if (myid.eq.0) write(tmp_unit,"(I11, E20.10,F30.10)") lt, time,T
    end if

    if (eng_out.gt.0.and.MOD(lt,eng_out).eq.0) then
       call kineticenergy(V,KE)
       call potentialenergy(X,PE)
       call momentum(V,Pm)
       KE = KE/epsSi/REAL(Natm)
       PE = PE/epsSi/REAL(Natm)
       write(eng_unit,"(E20.10,3F25.15,3E10.2)") time,KE,PE,KE+PE,Pm(:)
    end if

    if (Tbr_out.gt.0) then
       call compute_Tprof(V,Tbar)
       TbarS = (Tbar + REAL(NTbar)*TbarS)/REAL(NTbar+1)
       NTbar = NTbar + 1
       if (MOD(lt,Tbr_out).eq.0.and.myid.eq.0) then
          write(fn,"('D/Tbar.out',I9.9)")lt
          open(Tbr_unit,file=fn)
          do j = 1,Nlat(1)
             write(Tbr_unit,"(3E20.10)")Xlat(j),TbarS(j)
          end do
          close(Tbr_unit)
          NTbar = 0
          TbarS = 0.
       end if
    end if

  END SUBROUTINE diagnostics

END MODULE out
