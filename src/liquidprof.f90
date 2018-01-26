program solidprof

  integer :: i, j, k, ti,tf,df, Natm, fi, s_nr, s_nz, s_nt, s_cnt1, s_cnt2, s_cnt, sym_cnt, s_ri, s_ti, s_zi, s_dum
  real, allocatable, dimension (:,:) :: X, Xf, V, F, s_rr, s_2, s_3, s_xp1, s_xp2, s_vol, s_symm
  real, allocatable, dimension(:)  :: kin_eng
  real, dimension(3) :: Lb, s_ctr, diffx, s_xi
  integer, allocatable, dimension (:) :: atype, P, id, tst
  real :: s_ro, s_dr, s_dz, s_dt, pi, halfuc, s_rxi, s_txi, diffr, ev, kb
  character(30) :: fn

  Lb = (/ 2.7155365E-8, 2.7155365E-8, 1.0862146E-8 /)
  Natm = 401000
  ev = 1.60217646e-19
  kb = 1.3806503e-23
  halfuc = 5.431073E-10/2.0
  allocate(Xf(Natm,3), X(Natm,3), V(Natm,3), F(Natm,3), kin_eng(Natm), atype(Natm), P(Natm), id(Natm), tst(Natm))
  tst = 0
  
  tf = 4900000 !5868000
  ti = 4312000
  df = 4000
  
  write(fn,"('D/atoms',I9.9,'.dat')")tf
  open(29,file=fn)

  do i=1,Natm
      read(29,*) Xf(i,:) !,atype(i),P(i),kin_eng(i),V(i,:),F(i,:),id(i)
  end do
  close(29)


  s_ro = 1.0E-09
  s_dr = 1.0E-09
  s_dz = 1.0E-09
  s_dt = 30.0
  s_ctr = Lb/2.0
  s_nr = INT((s_ctr(1)-s_ro)/s_dr)
  s_nz = INT(Lb(3)/s_dz)+1
  s_nt = INT(360.0/s_dt)
  pi = 4.0*atan(1.0)
  s_cnt1 = 0
  s_cnt2 = 0

  print*, s_nr, s_nt, s_nz
  allocate(s_rr(s_nr*s_nt*s_nz,12), s_symm(s_nr*s_nz,12))
  s_rr = 0.0
  s_symm = 0.0
  
  s_cnt = 0
  do i=1,s_nz
      do j=1,s_nt
          do k=1,s_nr
              s_cnt = s_cnt+1
              s_rr(s_cnt,:)=(/ k*s_dr, j*s_dt, i*s_dz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
             
          end do
      end do
  end do
  
  fi = 1000000
  !do fi=ti,tf,df
      write(fn,"('D/atoms',I9.9,'.dat')")fi
      open(29,file=fn)
      do i=1,Natm
          read(29,*) X(i,:),atype(i),P(i),kin_eng(i) !,V(i,:),F(i,:),id(i)

          dt = (fi-ti+df)*0.2E-15
          s_xi(1:2) = X(i,1:2)-s_ctr(1:2); s_xi(3)=X(i,3)
          s_rxi = sqrt(sum(s_xi(1:2)**2))
          
          if(s_xi(1).eq.0) then
              if(s_xi(2).gt.0) then
                  s_txi = 90.0
              else
                  s_txi = 270.0
              end if
          else
              s_txi = atan(s_xi(2)/s_xi(1))
              if (s_xi(1).lt.0.0) s_txi = s_txi+pi
              if (s_xi(2).lt.0.0 .and. s_xi(1).ge.0.0) s_txi = s_txi+ 2*pi
              s_txi = s_txi*180.0/pi
          end if

          !print*, i,s_symm(1,1),'!!!!!!!!!!!!'
          !if(s_symm(1,1).gt. 1.0e-8) stop
          
          s_ri = int(s_rxi/s_dr)
          s_ti = int(s_txi/s_dt)
          s_zi = int(s_xi(3)/s_dz)
          if(s_ri .le. s_nr .and. s_ti .le. s_nt .and. s_zi .le. s_nz) then
              s_dum = s_ri + s_ti*s_nr + s_zi*s_nt*s_nr + 1
              s_rr(s_dum,12) = s_rr(s_dum,12) + 1
              s_rr(s_dum,4) = s_rr(s_dum,4) + kin_eng(i)*ev*2.0/3.0/kb
              !print*, s_rr(s_dum,4),s_rr(s_dum,12), kin_eng(i)*ev*2.0/3.0/kb
              !s_rr(s_dum,7) = s_rr(s_dum,7) + dt
              !if(dt > 1.0) then
              !    print *, fi, i, diffr, dt
              !end if
          end if
                  !print*, i, diffx/dt, diffr
      end do
      close(29)
      !print*, 'done ', tf - fi
      !end do


      sym_cnt = 0
      do i=1,s_nz
          do k=1,s_nr
              sym_cnt = sym_cnt+1
              s_symm(sym_cnt,:)=(/ k*s_dr, i*s_dz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
              print*, s_symm(sym_cnt,1:2)
          end do
      end do
      print*, sym_cnt, s_cnt
      
  sym_cnt = 0
  s_cnt = 0
  do i=1,s_nz
      do j=1,s_nt
          do k=1,s_nr
              sym_cnt = (i-1)*s_nr + k
              s_cnt = s_cnt+1

              if(s_rr(s_cnt,12) .gt. 0) then
                  s_rr(s_cnt,4)=s_rr(s_cnt,4)/s_rr(s_cnt,12)
                  !s_rr(s_cnt,7)=s_rr(s_cnt,7)/s_rr(s_cnt,12)
              end if

              print*, s_symm(sym_cnt,1:2)
              s_symm(sym_cnt,4) = s_symm(sym_cnt,4) + s_rr(s_cnt,4)

          end do
      end do
  end do

  print *, 'writing...'
  open(52, file='liquidprof.dat')
  do i=1,s_cnt
      write(52,"(12E20.10)")s_rr(i,:)
  end do
  close (52)

  open(53, file='symmliquid.dat')
  do i=1,sym_cnt
      write(53,"(12E20.10)")s_symm(i,1:3), s_symm(i,4:12)/(360.0/s_dt)
      print*, i, s_symm(i,1:2)
  end do
  close(53)

  print*, 'done exiting'
  call exit(0)
  !deallocate(s_rr, s_symm)
              
end program solidprof
  
