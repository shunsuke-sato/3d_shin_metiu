module global_variables
  implicit none
! mathematical constants
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! model parameters
  real(8), :: Ldist_m
  real(8), :: Rf_m, Rr_m, Rl_m
  real(8), :: mass_p

! grid coordinates
  integer :: nx_p, nx_e
  real(8) :: Lx_p, Lx_e  ! f(-Lx_p:Lx_p)
  real(8) :: hx_p, hx_e

  real(8),allocatable :: wfn(:,:,:)
  real(8),allocatable :: tpsi(:,:,:),htpsi(:,:,:)

  real(8),allocatable :: vpot(:,:,:)

  contains
!    
    function erf_x(x) result(y)
      real(8),intent(in) :: x
      real(8),parameter :: epsilon_s = 1d-3
      real(8) :: y

      if(abs(x) > epsilon_s)then
        y = erf(x)/x
      else
        y = 2d0/sqrt(pi)*( 1d0 - x**2/3d0 + x**4/10d0 - x**6/42d0 + x**8/216d0)
      end if

    end function erf_x


end module global_variables
!===============================================================================
program main
  use global_variables
  implicit none

  call input
  call initialization


  call ground_state

end program main
!===============================================================================
subroutine input
  use global_variables
  implicit none

! model parameter
  Ldist_m = 19d0
  Rf_m = 3.5d0
  Rr_m = 3.5d0
  Rl_m = 3.5d0
  mass_p = 1836d0


! grid parameter
  Lx_p = 9d0
  Lx_e = 20d0
  nx_p = 60
  nx_e = 200d0
  
  hx_p = Lx_p/nx_p
  hx_e = Lx_e/nx_e

end subroutine input
!===============================================================================
subroutine initialization
  use global_variables
  implicit none
  integer :: ixe1, ixe2, ixp
  real(8) :: xx_e1, xx_e2, xx_p

  allocate(wfn(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p))
  allocate(vpot(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p))

  allocate(tpsi(-nx_e-2:nx_e+2,-nx_e-2:nx_e+2,-nx_p-2:nx_p+2))
  tpsi = 0d0
  allocate(htpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p))


  do ixp = -nx_p,nx_p
    xx_p = hx_p*ixp
    do ixe2 = -nx_e,nx_e
      xx_e2 = hx_e*ixe2
      do ixe1 = -nx_e,nx_e
        xx_e1 = hx_e*ixe1


        vpot(ixe1,ixe2,ixp) = 1d0/abs(0.5d0*Ldist_m-xx_p) + 1d0/abs(0.5d0*Ldist_m+xx_p) &
          -erf_x(abs(xx_e1-xx_p)/Rf_m)/Rf_m &
          -erf_x(abs(xx_e1-0.5d0*Ldist_m)/Rr_m)/Rr_m &
          -erf_x(abs(xx_e1+0.5d0*Ldist_m)/Rl_m)/Rl_m &
          -erf_x(abs(xx_e2-xx_p)/Rf_m)/Rf_m &
          -erf_x(abs(xx_e2-0.5d0*Ldist_m)/Rr_m)/Rr_m &
          -erf_x(abs(xx_e2+0.5d0*Ldist_m)/Rl_m)/Rl_m 


      end do
    end do
  end do


end subroutine initialization
!===============================================================================
subroutine ground_state
  use global_variables
  implicit none
  integer :: ixp,ixe1,ixe2
  real(8) :: ss

! initialize
  do ixp = -nx_p,nx_p
    do ixe2 = -nx_e, nx_e
      do ixe1 = -nx_e, nx_e
        call random_number(wfn(ixe1,ixe2,ixp))
      end do
    end do
  end do

! symmetrize
  do ixp = -nx_p,nx_p
    do ixe2 = -nx_e, nx_e
      do ixe1 = -nx_e, nx_e
        tpsi(ixe1,ixe2,ixp) = wfn(ixe1,ixe2,ixp) + wfn(ixe2,ixe1,ixp)
      end do
    end do
  end do
  wfn(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p) = tpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)

! normalize
  ss = sum(wfn**2)*hx_p*hx_e**2
  wfn = wfn/sqrt(ss)

  call conjugate_gradient


end subroutine ground_state
!===============================================================================
subroutine conjugate_gradient
  use global_variables
  implicit none



end subroutine conjugate_gradient
!===============================================================================
subroutine hpsi
  use global_variables
  implicit none
  real(8),parameter :: c0 = -5d0/2d0, c1 = 4d0/3d0, c2 = -1d0/12d0
  real(8) :: l0p,l1p,l2p,l0e,l1e,l2e,l0pe
  integer :: ixp, ixe1, ixe2

  l0p = -0.5d0*c0/(hx_p**2*mass_p)
  l1p = -0.5d0*c1/(hx_p**2*mass_p)
  l2p = -0.5d0*c2/(hx_p**2*mass_p)

  l0e = -0.5d0*c0/(hx_e**2)
  l1e = -0.5d0*c1/(hx_e**2)
  l2e = -0.5d0*c2/(hx_e**2)

  l0pe = l0p + 2*l0e

  do ixp = -nx_p,nx_p
    do ixe2 = -nx_e, nx_e
      do ixe1 = -nx_e, nx_e


        htpsi(ixe1,ixe2,ixp) = l0pe*tpsi(ixe1,ixe2,ixp) &
          +c1p*(tpsi(ixe1,ixe2,ixp+1)+tpsi(ixe1,ixe2,ixp-1)) &
          +c2p*(tpsi(ixe1,ixe2,ixp+2)+tpsi(ixe1,ixe2,ixp-2)) &
          +c1e*(tpsi(ixe1+1,ixe2,ixp)+tpsi(ixe1-1,ixe2,ixp) &
               +tpsi(ixe1,ixe2+1,ixp)+tpsi(ixe1,ixe2-1,ixp) ) &
          +c2e*(tpsi(ixe1+2,ixe2,ixp)+tpsi(ixe1-2,ixe2,ixp) &
               +tpsi(ixe1,ixe2+2,ixp)+tpsi(ixe1,ixe2-2,ixp) )

      end do
    end do
  end do

  htpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p) = htpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p) &
    + vpot(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)*tpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)


end subroutine hpsi
!===============================================================================
!===============================================================================
!===============================================================================
