MODULE stats 

  USE prms
  USE data
  USE stillweb
  USE parallel
  USE par

  IMPLICIT none 

CONTAINS  

  SUBROUTINE totalenergy(X,V,TE)
    real, dimension(Natm,3)  :: X,V
    real                     :: TE,PE,KE

    call potentialenergy(X,PE)
    call kineticenergy(V,KE)

    TE = PE + KE

  END SUBROUTINE totalenergy

  SUBROUTINE potentialenergy(X,PE)
    real, dimension(Natm,3)  :: X
    real, dimension(Natm)    :: u
    real                     :: PE,PE_l

    integer                  :: ierr

    call si_potential(X,PE_l,u)
    !print *,myid,":  ",PE_l

    call MPI_REDUCE(PE_l,PE,1,MPI_REAL8,MPI_SUM, 0, MPI_COMM_WORLD, ierr)    

  END SUBROUTINE potentialenergy

  SUBROUTINE kineticenergyL(V,k)
    real, dimension(Natm,3)  :: V
    real, dimension(Nlat(1)) :: k

    integer                  :: i

    k = 0.
    do i = 1,Natm
       k(Lijk(i,1)) = k(Lijk(i,1)) + mass(i)*SUM(V(i,:)**2)
    end do
    k = k /2.

  END SUBROUTINE kineticenergyL
    
  SUBROUTINE kineticenergy(V,KE)
    real, dimension(Natm,3)  :: V
    real                     :: KE,KE_l

    integer                  :: ierr
    integer                  :: i,ii

    KE_l = 0.
    do ii = 1,Nl
       i = il(ii)
       KE_l = KE_l + mass(i)*(V(i,1)**2 + V(i,2)**2 + V(i,3)**2)
     end do
    KE_l = KE_l/2.

    call MPI_REDUCE(KE_l,KE,1,MPI_REAL8,MPI_SUM, 0, MPI_COMM_WORLD, ierr)    

  END SUBROUTINE kineticenergy

!  Compute the momentum
  SUBROUTINE momentum(V,Pm)
    real, dimension(Natm,3)  :: V
    real, dimension(3)       :: Pm,Pm_l
    integer                  :: i,ii

    integer                  :: ierr

    Pm_l = 0.
    do ii = 1,Nl
       i = il(ii)
       Pm_l = Pm_l + mass(i)*V(i,:)
    end do
    call MPI_REDUCE(Pm_l,Pm,3,MPI_REAL8,MPI_SUM, 0, MPI_COMM_WORLD, ierr)    

  END SUBROUTINE momentum

END MODULE stats

