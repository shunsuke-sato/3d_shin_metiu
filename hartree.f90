module global_variables
  implicit none
! mathematical constants
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! physical constants
  real(8),parameter :: fs = 0.024189d0

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

  real(8),allocatable :: rho_e(:), rho_p(:)
  real(8),allocatable :: xe_l(:), xp_l(:)

  complex(8),allocatable :: zwfn_e(:,:), zwfn_p(:)
  complex(8),allocatable :: ztpsi_e(:,:),zhtpsi_e(:,:)
  complex(8),allocatable :: ztpsi_p(:),zhtpsi_p(:)

  integer :: nt
  real(8) :: Tprop, dt
  real(8) :: E0, Tpulse, omega0


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


  call ground_state

  call reduced_density


  call time_propagation

end program main
!-----------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  write(*,*)"Start: input"
! model parameter
  Ldist_m = 19d0
  Rf_m = 5d0
  Rr_m = 4d0
  Rl_m = 3.1d0
  Re_m = Rf_m

!  Ldist_m = 19d0 
!  Rf_m = 4d0
!  Rr_m = 4d0
!  Rl_m = 4d0
!  Re_m = 4d0  ! 5.6d0  5.7d0
  mass_p = 1836d0



! grid parameter
  Lx_p = 0.5d0*(Ldist_m-1d0)
  Lx_e = 25d0
  nx_p = 9*4
  nx_e = 25*4
  
  hx_p = Lx_p/nx_p
  hx_e = Lx_e/nx_e


! time propagation
  dt = 0.02d0
  Tprop = 40d0*fs
  nt = aint(Tprop/dt) + 1

  E0 = 0.002d0
  Tpulse = 20d0*fs
  omega0 = 0.055d0




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

  allocate(rho_e(-nx_e:nx_e))
  allocate(rho_p(-nx_p:nx_p))
  allocate(vh_e(-nx_e:nx_e))
  allocate(vh_p(-nx_p:nx_p))

  allocate(tpsi_e(-nx_e-2:nx_e+2,-nx_e-2:nx_e+2))
  allocate(tpsi_p(-nx_p-2:nx_p+2))
  tpsi_e = 0d0; tpsi_p = 0d0
  allocate(htpsi_e(-nx_e:nx_e,-nx_e:nx_e))
  allocate(htpsi_p(-nx_p:nx_p))


  allocate(zwfn_e(-nx_e:nx_e,-nx_e:nx_e))
  allocate(zwfn_p(-nx_p:nx_p))

  allocate(ztpsi_e(-nx_e-2:nx_e+2,-nx_e-2:nx_e+2))
  allocate(ztpsi_p(-nx_p-2:nx_p+2))
  ztpsi_e = 0d0; ztpsi_p = 0d0
  allocate(zhtpsi_e(-nx_e:nx_e,-nx_e:nx_e))
  allocate(zhtpsi_p(-nx_p:nx_p))

  allocate(xe_l(-nx_e:nx_e),xp_l(-nx_p:nx_p))


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

      vint_ep(ixe1,ixp) = -1d0*erf_x(abs(xx_e1-xx_p)/Rf_m)/Rf_m
      vint_pe(ixp,ixe1) = -1d0*erf_x(abs(xx_e1-xx_p)/Rf_m)/Rf_m

    end do
  end do

  do ixp = -nx_p,nx_p
    xx_p = hx_p*ixp
    xp_l(ixp) = xx_p
  end do

  do ixe1 = -nx_e,nx_e
    xx_e1 = hx_e*ixe1
    xe_l(ixe1) = xx_e1
  end do



  write(*,*)"End: initialization"


end subroutine initialization
!-----------------------------------------------------------------------------
subroutine ground_state
  use global_variables
  implicit none
  integer :: iscf, nscf

  nscf = 200
  call initialize_wfn

  call calc_Hartree_pot


  do iscf = 1, nscf

    call cg_e(5)
    call calc_Hartree_pot
    call cg_p(5)
    call calc_Hartree_pot

  end do
  

end subroutine ground_state
!-----------------------------------------------------------------------------
subroutine initialize_wfn
  use global_variables
  implicit none
  integer :: ixp, ixe1, ixe2
  real(8) :: ss

  do ixp = -nx_p,nx_p
    call random_number(wfn_p(ixp))
  end do
  ss = sum(wfn_p**2)*hx_p
  wfn_p = wfn_p/sqrt(ss)


  do ixe2 = -nx_e, nx_e
    do ixe1 = -nx_e, nx_e
      call random_number(wfn_e(ixe1,ixe2))
    end do
  end do


  do ixe2 = -nx_e, nx_e
    do ixe1 = -nx_e, nx_e
      tpsi_e(ixe1,ixe2) = wfn_e(ixe1,ixe2) + wfn_e(ixe2,ixe1)
    end do
  end do

  wfn_e(-nx_e:nx_e,-nx_e:nx_e) = tpsi_e(-nx_e:nx_e,-nx_e:nx_e)
  ss = sum(wfn_e**2)*hx_e**2
  wfn_e = wfn_e/sqrt(ss)



end subroutine initialize_wfn
!-----------------------------------------------------------------------------
subroutine calc_Hartree_pot
  use global_variables
  implicit none
  integer :: ixp, ixe

  do ixp = -nx_p,nx_p
    rho_p(ixp) = wfn_p(ixp)**2
  end do

  do ixe = -nx_e,nx_e
    rho_e(ixe) = 2d0*sum(wfn_e(:,ixe)**2)*hx_e
  end do

  vh_e(:) = matmul(vint_ep,rho_p)*hx_p
  vh_p(:) = matmul(vint_pe,rho_e)*hx_e
  


end subroutine calc_Hartree_pot
!-----------------------------------------------------------------------------
subroutine cg_e(ncg)
  use global_variables
  implicit none
  integer,intent(in) :: ncg
  integer :: icg
  real(8) :: psi(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: xi(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: phi(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: phi_old(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: lambda, xixi, xixi_old, gamma
  real(8) :: ss, aa, bb, theta 

  psi = wfn_e

  tpsi_e(-nx_e:nx_e,-nx_e:nx_e) = psi(-nx_e:nx_e,-nx_e:nx_e)
  call hpsi_e
  lambda = sum(psi*htpsi_e)*hx_e**2
  xi = lambda*psi-htpsi_e

  xixi = sum(xi**2)*hx_e**2
  xixi_old = xixi

  do icg = 0, ncg

    if(icg == 0)then
      gamma = 0d0
      phi_old = 0d0
    else
      gamma = xixi/xixi_old
      xixi_old = xixi
    end if

    phi = xi + gamma* phi_old
    phi_old = phi
    ss = sum(phi*psi)*hx_e**2
    phi = phi -ss*psi
    ss = sum(phi**2)*hx_e**2
    phi = phi/sqrt(ss)

    bb = 2d0*sum(phi*htpsi_e)*hx_e**2
    tpsi_e(-nx_e:nx_e,-nx_e:nx_e) = phi(-nx_e:nx_e,-nx_e:nx_e)
    call hpsi_e
    aa = sum(phi*htpsi_e)*hx_e**2 -lambda
    aa = - aa
    if(aa /= 0d0)then
      theta = 0.5d0*atan(bb/aa)
    else
      theta = 0.25d0*pi
    end if
    psi = cos(theta)*psi + sin(theta)*phi

    ss = sum(psi**2)*hx_e**2
    psi = psi/sqrt(ss)
    if(icg == ncg)exit

    tpsi_e(-nx_e:nx_e,-nx_e:nx_e) = psi(-nx_e:nx_e,-nx_e:nx_e)
    call hpsi_e
    lambda = sum(psi*htpsi_e)*hx_e**2
    xi = lambda*psi - htpsi_e

    xixi = sum(xi**2)*hx_e**2

    write(*,*)'e',icg,lambda,xixi

  end do

   wfn_e = psi

end subroutine cg_e
!-----------------------------------------------------------------------------
subroutine cg_p(ncg)
  use global_variables
  implicit none
  integer,intent(in) :: ncg
  integer :: icg
  real(8) :: psi(-nx_p:nx_p)
  real(8) :: xi(-nx_p:nx_p)
  real(8) :: phi(-nx_p:nx_p)
  real(8) :: phi_old(-nx_p:nx_p)
  real(8) :: lambda, xixi, xixi_old, gamma
  real(8) :: ss, aa, bb, theta 

  psi = wfn_p

  tpsi_p(-nx_p:nx_p) = psi(-nx_p:nx_p)
  call hpsi_p
  lambda = sum(psi*htpsi_p)*hx_p
  xi = lambda*psi-htpsi_p

  xixi = sum(xi**2)*hx_p
  xixi_old = xixi

  do icg = 0, ncg

    if(icg == 0)then
      gamma = 0d0
      phi_old = 0d0
    else
      gamma = xixi/xixi_old
      xixi_old = xixi
    end if

    phi = xi + gamma* phi_old
    phi_old = phi
    ss = sum(phi*psi)*hx_p
    phi = phi -ss*psi
    ss = sum(phi**2)*hx_p
    phi = phi/sqrt(ss)

    bb = 2d0*sum(phi*htpsi_p)*hx_p
    tpsi_p(-nx_p:nx_p) = phi(-nx_p:nx_p)
    call hpsi_p
    aa = sum(phi*htpsi_p)*hx_p -lambda
    aa = - aa
    if(aa /= 0d0)then
      theta = 0.5d0*atan(bb/aa)
    else
      theta = 0.25d0*pi
    end if
    psi = cos(theta)*psi + sin(theta)*phi

    ss = sum(psi**2)*hx_p
    psi = psi/sqrt(ss)
    if(icg == ncg)exit

    tpsi_p(-nx_p:nx_p) = psi(-nx_p:nx_p)
    call hpsi_p
    lambda = sum(psi*htpsi_p)*hx_p
    xi = lambda*psi - htpsi_p

    xixi = sum(xi**2)*hx_p

    write(*,*)'p',icg,lambda,xixi

  end do

   wfn_p = psi

end subroutine cg_p
!-----------------------------------------------------------------------------
subroutine hpsi_e
  use global_variables
  implicit none
  real(8),parameter :: c0 = -5d0/2d0, c1 = 4d0/3d0, c2 = -1d0/12d0
  real(8) :: l0,l1,l2
  integer :: ixp, ixe1, ixe2

  l0 = -2d0*0.5d0*c0/(hx_e**2)
  l1 = -0.5d0*c1/(hx_e**2)
  l2 = -0.5d0*c2/(hx_e**2)

  do ixe2 = -nx_e, nx_e
    do ixe1 = -nx_e, nx_e

      htpsi_e(ixe1,ixe2) = l0*tpsi_e(ixe1,ixe2) &
          +l1*(tpsi_e(ixe1+1,ixe2)+tpsi_e(ixe1-1,ixe2) &
              +tpsi_e(ixe1,ixe2+1)+tpsi_e(ixe1,ixe2-1) ) &
          +l2*(tpsi_e(ixe1+2,ixe2)+tpsi_e(ixe1-2,ixe2) &
              +tpsi_e(ixe1,ixe2+2)+tpsi_e(ixe1,ixe2-2) )

    end do
  end do

  do ixe2 = -nx_e, nx_e
    do ixe1 = -nx_e, nx_e
      htpsi_e(ixe1,ixe2) = htpsi_e(ixe1,ixe2) &
        + (vpot_e(ixe1,ixe2)+vh_e(ixe1)+vh_e(ixe2))*tpsi_e(ixe1,ixe2)
    end do
  end do


end subroutine hpsi_e
!-----------------------------------------------------------------------------
subroutine hpsi_p
  use global_variables
  implicit none
  real(8),parameter :: c0 = -5d0/2d0, c1 = 4d0/3d0, c2 = -1d0/12d0
  real(8) :: l0,l1,l2
  integer :: ixp

  l0 = -0.5d0*c0/(hx_p**2*mass_p)
  l1 = -0.5d0*c1/(hx_p**2*mass_p)
  l2 = -0.5d0*c2/(hx_p**2*mass_p)

  do ixp = -nx_p, nx_p

    htpsi_p(ixp) = l0*tpsi_p(ixp) &
          +l1*(tpsi_p(ixp+1)+tpsi_p(ixp-1)) &
          +l2*(tpsi_p(ixp+2)+tpsi_p(ixp-2)) 

  end do

  do ixp = -nx_p, nx_p
    htpsi_p(ixp) = htpsi_p(ixp) &
      + (vpot_p(ixp)+vh_p(ixp))*tpsi_p(ixp)
  end do


end subroutine hpsi_p
!-----------------------------------------------------------------------------
subroutine reduced_density
  use global_variables
  implicit none
  integer :: ixp, ixe1,ixe2

  do ixp = -nx_p,nx_p
    rho_p(ixp) = wfn_p(ixp)**2
  end do

  do ixe1 = -nx_e,nx_e
    rho_e(ixe1) = 2d0*sum(wfn_e(:,ixe1)**2)*hx_e
  end do

  vh_e(:) = matmul(vint_ep,rho_p)*hx_p
  vh_p(:) = matmul(vint_pe,rho_e)*hx_e


  open(20,file="density_1e_h.out")
  do ixe1 = -nx_e, nx_e
    write(20,"(999e26.16e3)")ixe1*hx_e,rho_e(ixe1)*0.5d0,vh_e(ixe1)
  end do
  close(20)

  open(20,file="density_1p_h.out")
  do ixp = -nx_p, nx_p
    write(20,"(999e26.16e3)")ixp*hx_p,rho_p(ixp),vh_p(ixp)
  end do
  close(20)

  open(20,file="density_2pe_h.out")
  do ixp = -nx_p, nx_p
    do ixe1 = -nx_e, nx_e
      write(20,"(999e26.16e3)")ixp*hx_p,ixe1*hx_e,rho_p(ixp)*rho_e(ixe1)*0.5d0
    end do
    write(20,*)
  end do
  close(20)

  open(20,file="density_2ee_h.out")
  do ixe1 = -nx_e, nx_e
    do ixe2 = -nx_e, nx_e
      write(20,"(999e26.16e3)")ixe1*hx_e,ixe2*hx_e,wfn_e(ixe1,ixe2)**2
    end do
    write(20,*)
  end do
  close(20)

end subroutine reduced_density
!-----------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it

  zwfn_e = wfn_e
  zwfn_p = wfn_p


  do it = 0, nt
    
    call dt_evolve(it)

  end do


end subroutine time_propagation
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
