module global_variables
  implicit none
! mathematical constants
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! model parameters
  real(8) :: Ldist_m
  real(8) :: Rf_m, Rr_m, Rl_m, Re_m
  real(8) :: mass_p

! grid coordinates
  integer :: nx_p, nx_e
  real(8) :: Lx_p, Lx_e  ! f(-Lx_p:Lx_p)
  real(8) :: hx_p, hx_e

  real(8),allocatable :: wfn_e(:,:), wfn_p(:)
  real(8),allocatable :: tpsi_e(:,:),htpsi_e(:,:)
  real(8),allocatable :: tpsi_p(:),htpsi_p(:)

  real(8),allocatable :: vpot_e(:,:), vh_e(:)
  real(8),allocatable :: vpot_p(:), vh_p(:)
  real(8),allocatable :: vint_ep(:,:), vint_pe(:,:)

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
!-----------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialization


!  call ground_state

!  call reduced_density_matrix

end program main
!-----------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  write(*,*)"Start: input"
! model parameter
!  Ldist_m = 19d0
!  Rf_m = 5d0
!  Rr_m = 4d0
!  Rl_m = 3.1d0
!  Re_m = Rf_m

  Ldist_m = 19d0 
  Rf_m = 4d0
  Rr_m = 4d0
  Rl_m = 4d0
  Re_m = 5.65d0  ! 5.6d0  5.7d0
  mass_p = 1836d0



! grid parameter
  Lx_p = 0.5d0*(Ldist_m-1d0)
  Lx_e = 25d0
  nx_p = 9*8
  nx_e = 25*10
  
  hx_p = Lx_p/nx_p
  hx_e = Lx_e/nx_e

  write(*,*)"End: input"

end subroutine input
!-----------------------------------------------------------------------------
subroutine initialization
  use global_variables
  implicit none
  integer :: ixe1, ixe2, ixp
  real(8) :: xx_e1, xx_e2, xx_p

  write(*,*)"Start: initialization"

  allocate(wfn_e(-nx_e:nx_e,-nx_e:nx_e))
  allocate(wfn_p(-nx_p:nx_p))
  allocate(vpot_e(-nx_e:nx_e,-nx_e:nx_e))
  allocate(vpot_p(-nx_p:nx_p))
  allocate(vint_ep(-nx_e:nx_e,-nx_p:nx_p))
  allocate(vint_pe(-nx_p:nx_p,-nx_e:nx_e))

  allocate(tpsi_e(-nx_e-2:nx_e+2,-nx_e-2:nx_e+2))
  allocate(tpsi_p(-nx_p-2:nx_p+2))
  tpsi_e = 0d0; tpsi_p = 0d0
  allocate(htpsi_e(-nx_e:nx_e,-nx_e:nx_e))
  allocate(htpsi_p(-nx_p:nx_p))


  do ixe2 = -nx_e,nx_e
    xx_e2 = hx_e*ixe2
    do ixe1 = -nx_e,nx_e
      xx_e1 = hx_e*ixe1

      vpot_e(ixe1,ixe2) = &
        -erf_x(abs(xx_e1-0.5d0*Ldist_m)/Rr_m)/Rr_m &
        -erf_x(abs(xx_e1+0.5d0*Ldist_m)/Rl_m)/Rl_m &
        -erf_x(abs(xx_e2-0.5d0*Ldist_m)/Rr_m)/Rr_m &
        -erf_x(abs(xx_e2+0.5d0*Ldist_m)/Rl_m)/Rl_m &
        +erf_x(abs(xx_e2-xx_e1)/Re_m)/Re_m 

    end do
  end do


  do ixp = -nx_p,nx_p
    xx_p = hx_p*ixp

    vpot_p(ixp) = 1d0/abs(0.5d0*Ldist_m-xx_p) + 1d0/abs(0.5d0*Ldist_m+xx_p) 

  end do



  do ixp = -nx_p,nx_p
    xx_p = hx_p*ixp
    do ixe1 = -nx_e,nx_e
      xx_e1 = hx_e*ixe1

      vint_ep(ixe1,ixp) = -erf_x(abs(xx_e1-xx_p)/Rf_m)/Rf_m
      vint_pe(ixp,ixe1) = -erf_x(abs(xx_e1-xx_p)/Rf_m)/Rf_m


    end do
  end do

  write(*,*)"End: initialization"


end subroutine initialization
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
