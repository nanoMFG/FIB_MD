program meanvel_relax

  integer :: i,j,k,Natm,fl,startf, endf, df, atype, p
  real :: dx,dy,dz,uc, ke
  real, dimension(3) :: lb, lbs
  integer, dimension(3) :: ijk,npc,nc
  real, allocatable, dimension(:,:) :: X,V,F,xl,vl
  real, allocatable, dimension(:)   :: kin
  integer, allocatable, dimension(:):: id,cnt
  character(100) :: fn, cmd
  real, allocatable, dimension(:,:)  :: s_rr, s_xp1, s_xp2, s_2, s_3, s_vol
  real                               :: s_dr, s_dt, s_dz, s_ro
  Integer                            :: s_nr, s_nt, s_nz, s_cnt1, s_cnt2
  real, dimension(3)                 :: s_ctr
  real                               :: s_rxi, s_rxj, s_rxk, s_txi, s_txj, s_txk
  integer                            :: s_ri, s_ti, s_zi, s_rj, s_tj, s_zj,s_rk, s_tk, s_zk, s_cnt
  real, dimension(3)                 :: s_xi, s_xj, s_xk, fi, fj, fk
    
  
  Natm = 401000
  uc = 5.4310729e-10
  nc(1) = 50; nc(2) = 50; nc(3) = 20;
  lb = uc*nc + 1.0E-10
  lb(3) = lb(3)*10.0
  startf = 4312000
  endf = 5868000
  df = 4000
  
  s_ro = 0.0
  s_dr = 1.0E-09
  s_dz = 1.0E-09
  s_dt = 30.0
  s_ctr = Lb/2.0
  s_nr = INT((s_ctr(1)-s_ro)/s_dr)
  s_nz = INT(Lb(3)/s_dz/10.0)+1
  s_nt = INT(360.0/s_dt)
  pi = 4.0*atan(1.0)
  s_cnt1 = 0
  s_cnt2 = 0

  print*, s_nr, s_nt, s_nz
  allocate(s_rr(s_nr*s_nt*s_nz,11), X(Natm,3),V(Natm,3))
  s_rr = 0.0
  
  s_cnt = 0
  do i=0,s_nz-1
      do j=0,s_nt-1
          do k=0,s_nr-1
              s_cnt = s_cnt+1
              s_rr(s_cnt,:)=(/ k*s_dr, j*s_dt, i*s_dz, 0.0, 0.0 /)
          end do
      end do
  end do
  
  do fl=startf, endf, df
      write(fn,"('D/atoms00', I7.7,'.dat')")fl
      open(23, file=fn)
      do i = 1,Natm

          read(23,*)X(i,:), atype, p, ke, V(i,:)
          
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
          s_ri = int(s_rxi/s_dr)
          s_ti = int(s_txi/s_dt)
          s_zi = int(s_xi(3)/s_dz)
          if(s_ri.le.s_nr .and. s_ti.le.s_nt .and. s_zi.le.s_nz) then
              s_dum = s_ri + s_ti*s_nr + s_zi*s_nt*s_nr + 1
              
              s_rr(s_dum,4) = s_rr(s_dum,4) + sum(V(i,1:2)*s_xi(1:2))/sqrt(sum(s_xi(1:2)**2))
              s_rr(s_dum,5) = s_rr(s_dum,5) + (-V(i,1)*sind((s_ti+1)*s_dt) + V(i,2)*cosd((s_ti+1)*s_dt))
              s_rr(s_dum,6) = s_rr(s_dum,6) + V(i,3)
             
              s_rr(s_dum,7) = s_rr(s_dum,7) + V(i,1)
              s_rr(s_dum,8) = s_rr(s_dum,8) + V(i,2)
              s_rr(s_dum,9) = s_rr(s_dum,9) + V(i,3)
             
              s_rr(s_dum,10) = s_rr(s_dum,10)+1.0
             
          end if
      end do
      close(23)

      s_rr(:,4:9) = s_rr(:,4:9)/s_rr(s_dum,10)
      s_cnt = 0
      do i=1,s_nz
          do j=1,s_nt
              do k=1,s_nr
                  s_cnt = s_cnt+1
                  s_rr(s_cnt,11) = s_rr(s_cnt,10)/( pi*(k*s_dr*k*s_dr-(k-1)*s_dr*(k-1)*s_dr)*s_dt/360.0*s_dz  )
              end do
          end do
      end do

      write(fn,"('meanvel_relax/vel_',I7.7,'.dat')")fl
      open(50,file=fn)
      do i=1,s_cnt
          write(50,"(11E15.6)")s_rr(i,:)
      end do
      close(50)

      print*, "done ",real((fl-startf))/real((endf-startf))*100.0,"%"
      
  end do
  
end program meanvel_relax
