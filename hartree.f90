module global_variables
  implicit none
! mathematical constants
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! physical constants
  real(8),parameter :: fs = 1d0/0.024189d0

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
  real(8),allocatable :: Et(:)

! bo
  integer,parameter :: n_bo = 4
  real(8),allocatable :: wfn_bo(:,:,:,:)


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

  call calc_BO_surface
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
  real(8) :: pop_bo(n_bo), cohe_bo(n_bo,n_bo)
  real(8) :: Etot_t
  integer :: it

  zwfn_e = wfn_e
  zwfn_p = wfn_p

  call init_laser
  
  open(51,file='BO_pop_coherence.out')
  open(52,file='Etot_xe_xp_t.out')



  do it = 0, nt
    if(mod(it,200)==0)then
      call output_td_density(it)
    end if

    if(mod(it,10)== 0)then
      call calc_Etot_td(Etot_t,it)
      write(52,"(999e26.16e3)")it*dt,Etot_t,sum(rho_e*xe_l)/sum(rho_e) &
                                           ,sum(rho_p*xp_l)/sum(rho_p)
    end if


    if(mod(it,10) == 0)then
      call calc_bo_quantity(pop_bo,cohe_bo)
      write(51,"(999e26.16e3)")dt*it,pop_bo,cohe_bo(1,:)
    end if

    call dt_evolve(it)

  end do

  close(51)
  close(52)

end subroutine time_propagation
!-----------------------------------------------------------------------------
subroutine output_td_density(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ixp, ixe
  character(256) :: cit, cfilename

  do ixp = -nx_p,nx_p
    rho_p(ixp) = abs(zwfn_p(ixp))**2
  end do

  do ixe = -nx_e,nx_e
    rho_e(ixe) = 2d0*sum(abs(zwfn_e(:,ixe))**2)*hx_e
  end do

  write(cit,"(I9.9)")it
  cfilename = trim(cit)//'_rho_e.out'

  open(31,file=cfilename)
  write(31,"(A,2x,e26.16e3)")'# time=',it*dt/fs
  do ixe = -nx_e,nx_e
    write(31,"(999e26.16e3)")xe_l(ixe),rho_e(ixe)
  end do
  close(31)

  cfilename = trim(cit)//'_rho_p.out'

  open(31,file=cfilename)
  write(31,"(A,2x,e26.16e3)")'# time=',it*dt/fs
  do ixp = -nx_p,nx_p
    write(31,"(999e26.16e3)")xp_l(ixp),rho_p(ixp)
  end do
  close(31)


end subroutine output_td_density
!-----------------------------------------------------------------------------
! Propagator with the predictor-corrector method
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8) :: zwfn_e_pc(-nx_e:nx_e,-nx_e:nx_e)
  complex(8) :: zwfn_p_pc(-nx_p:nx_p)
  real(8) :: Et_t

! t -> t +dt/2
  Et_t = Et(it)
  call calc_Hartree_pot_td
  call dt_evolve_e(Et_t,dt*0.5d0)
  call dt_evolve_p(Et_t,dt*0.5d0)

  zwfn_e_pc = zwfn_e
  zwfn_p_pc = zwfn_p

! t+ dt/2 -> t + dt (predictor step)
  Et_t = Et(it+1)
  call calc_Hartree_pot_td
  call dt_evolve_e(Et_t,dt*0.5d0)
  call dt_evolve_p(Et_t,dt*0.5d0)

! t+ dt/2 -> t + dt (corrector step)
  zwfn_e = zwfn_e_pc
  zwfn_p = zwfn_p_pc
  Et_t = Et(it+1)
  call calc_Hartree_pot_td
  call dt_evolve_e(Et_t,dt*0.5d0)
  call dt_evolve_p(Et_t,dt*0.5d0)


end subroutine dt_evolve
!-----------------------------------------------------------------------------
subroutine dt_evolve_e(Et_in,dt_in)
  use global_variables
  implicit none
  real(8),intent(in) :: Et_in, dt_in
  integer,parameter :: nexp = 4
  integer :: iexp
  complex(8) :: zfact

  ztpsi_e(-nx_e:nx_e,-nx_e:nx_e) = zwfn_e(-nx_e:nx_e,-nx_e:nx_e)
  zfact = 1d0
  do iexp = 1, nexp
    zfact = zfact*(-zI*dt_in)/iexp
    call zhpsi_e(Et_in)
    zwfn_e = zwfn_e + zfact*zhtpsi_e
    ztpsi_e(-nx_e:nx_e,-nx_e:nx_e) = zhtpsi_e(-nx_e:nx_e,-nx_e:nx_e)
  end do


end subroutine dt_evolve_e
!-----------------------------------------------------------------------------
subroutine dt_evolve_p(Et_in,dt_in)
  use global_variables
  implicit none
  real(8),intent(in) :: Et_in, dt_in
  integer,parameter :: nexp = 4
  integer :: iexp
  complex(8) :: zfact

  ztpsi_p(-nx_p:nx_p) = zwfn_p(-nx_p:nx_p)
  zfact = 1d0
  do iexp = 1, nexp
    zfact = zfact*(-zI*dt_in)/iexp
    call zhpsi_p(Et_in)
    zwfn_p = zwfn_p + zfact*zhtpsi_p
    ztpsi_p(-nx_p:nx_p) = zhtpsi_p(-nx_p:nx_p)
  end do


end subroutine dt_evolve_p
!-----------------------------------------------------------------------------
subroutine zhpsi_e(Et_in)
  use global_variables
  implicit none
  real(8),intent(in) :: Et_in
  real(8),parameter :: c0 = -5d0/2d0, c1 = 4d0/3d0, c2 = -1d0/12d0
  real(8) :: l0,l1,l2
  integer :: ixp, ixe1, ixe2

  l0 = -2d0*0.5d0*c0/(hx_e**2)
  l1 = -0.5d0*c1/(hx_e**2)
  l2 = -0.5d0*c2/(hx_e**2)

!$omp parallel do private(ixe2,ixe1)
  do ixe2 = -nx_e, nx_e
    do ixe1 = -nx_e, nx_e

      zhtpsi_e(ixe1,ixe2) = l0*ztpsi_e(ixe1,ixe2) &
          +l1*(ztpsi_e(ixe1+1,ixe2)+ztpsi_e(ixe1-1,ixe2) &
              +ztpsi_e(ixe1,ixe2+1)+ztpsi_e(ixe1,ixe2-1) ) &
          +l2*(ztpsi_e(ixe1+2,ixe2)+ztpsi_e(ixe1-2,ixe2) &
              +ztpsi_e(ixe1,ixe2+2)+ztpsi_e(ixe1,ixe2-2) )

    end do
  end do

!$omp parallel do private(ixe2,ixe1)
  do ixe2 = -nx_e, nx_e
    do ixe1 = -nx_e, nx_e
      zhtpsi_e(ixe1,ixe2) = zhtpsi_e(ixe1,ixe2) &
        + (vpot_e(ixe1,ixe2)+vh_e(ixe1)+vh_e(ixe2) + Et_in*(xe_l(ixe1)+xe_l(ixe2)))&
        *ztpsi_e(ixe1,ixe2)
    end do
  end do


end subroutine zhpsi_e
!-----------------------------------------------------------------------------
subroutine zhpsi_p(Et_in)
  use global_variables
  implicit none
  real(8),intent(in) :: Et_in
  real(8),parameter :: c0 = -5d0/2d0, c1 = 4d0/3d0, c2 = -1d0/12d0
  real(8) :: l0,l1,l2
  integer :: ixp

  l0 = -0.5d0*c0/(hx_p**2*mass_p)
  l1 = -0.5d0*c1/(hx_p**2*mass_p)
  l2 = -0.5d0*c2/(hx_p**2*mass_p)

  do ixp = -nx_p, nx_p

    zhtpsi_p(ixp) = l0*ztpsi_p(ixp) &
          +l1*(ztpsi_p(ixp+1)+ztpsi_p(ixp-1)) &
          +l2*(ztpsi_p(ixp+2)+ztpsi_p(ixp-2)) 

  end do

  do ixp = -nx_p, nx_p
    zhtpsi_p(ixp) = zhtpsi_p(ixp) &
      + (vpot_p(ixp)+vh_p(ixp) -Et_in*xp_l(ixp))*ztpsi_p(ixp)
  end do


end subroutine zhpsi_p
!-----------------------------------------------------------------------------
subroutine init_laser
  use global_variables
  implicit none
  real(8) :: tt, xx
  integer :: it

  allocate(Et(0:nt+1))
  Et = 0d0

  do it = 0, nt+1
    tt = dt*it
    xx = tt - 0.5d0*Tpulse
    if(abs(xx) < 0.5d0*Tpulse)then
      Et(it) = E0*cos(pi*xx/Tpulse)**2*sin(omega0*xx)
    end if
    
  end do


  open(31,file='laser.out')
  do it = 0, nt+1
    write(31,"(999e26.16e3)")dt*it,Et(it)
  end do
  close(31)

end subroutine init_laser
!-----------------------------------------------------------------------------
subroutine calc_Hartree_pot_td
  use global_variables
  implicit none
  integer :: ixp, ixe

  do ixp = -nx_p,nx_p
    rho_p(ixp) = abs(zwfn_p(ixp))**2
  end do

  do ixe = -nx_e,nx_e
    rho_e(ixe) = 2d0*sum(abs(zwfn_e(:,ixe))**2)*hx_e
  end do

  vh_e(:) = matmul(vint_ep,rho_p)*hx_p
  vh_p(:) = matmul(vint_pe,rho_e)*hx_e
  


end subroutine calc_Hartree_pot_td
!-----------------------------------------------------------------------------
subroutine calc_BO_surface
  use global_variables
  implicit none
  integer,parameter :: ncg_bo = 100, nscf = 5
  integer :: ixp, iscf
  real(8) :: psi_t(-nx_e:nx_e,-nx_e:nx_e,n_bo)
  real(8) :: eps_t(n_bo),eps_bo(-nx_p:nx_p,n_bo)

  allocate(wfn_bo(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p,n_bo))

  do ixp = -nx_p, nx_p
    write(*,*)'ixp/nx_p',ixp,nx_p
    vh_e(:) = vint_ep(:,ixp)
    call random_number(psi_t)    
    do iscf = 1, nscf
      call cg_bo(psi_t, eps_t, ncg_bo)
!      call subspace_diag(psi_t, eps_t)
      call reordering(psi_t, eps_t)
    end do
    wfn_bo(:,:,ixp,:)  = psi_t
    eps_bo(ixp,:) = eps_t(:) + vpot_p(ixp)

  end do

  open(31,file='bo_surface.out')
  do ixp = -nx_p, nx_p
    write(31,"(999e26.16e3)")xp_l(ixp),eps_bo(ixp,:)
  end do
  close(31)

!  stop


end subroutine calc_BO_surface
!-----------------------------------------------------------------------------
subroutine cg_bo(psi_out, eps_out, ncg)
  use global_variables
  implicit none
  real(8),intent(inout) :: psi_out(-nx_e:nx_e,-nx_e:nx_e,n_bo)
  real(8),intent(out) :: eps_out(n_bo)
  integer,intent(in) :: ncg

  real(8),parameter :: res_c = 1d-10
  integer :: icg
  integer :: ib,ib2
  real(8) :: psi(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: xi(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: phi(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: phi_old(-nx_e:nx_e,-nx_e:nx_e)
  real(8) :: lambda, xixi, xixi_old, gamma
  real(8) :: ss, aa, bb, theta 
  real(8) :: psi_t(-nx_e:nx_e,-nx_e:nx_e)
  integer :: ix1, ix2

  do ib = 1, n_bo
    psi_t = psi_out(:,:,ib)
    do ix1 = -nx_e,nx_e
      do ix2 = -nx_e,nx_e
        psi(ix1,ix2) = 0.5d0*(psi_t(ix1,ix2)+psi_t(ix2,ix1))
      end do
    end do
    ss = sum(psi**2)*hx_e**2
    psi = psi/sqrt(ss)
    do ib2 = 1, ib-1
      ss = sum(psi*psi_out(:,:,ib2))*hx_e**2
      psi = psi -ss*psi_out(:,:,ib2)
    end do
    ss = sum(psi**2)*hx_e**2
    psi = psi/sqrt(ss)


    tpsi_e(-nx_e:nx_e,-nx_e:nx_e) = psi(-nx_e:nx_e,-nx_e:nx_e)
    call hpsi_e
    lambda = sum(psi*htpsi_e)*hx_e**2
    xi = lambda*psi-htpsi_e
    do ib2 = 1, ib-1
      ss = sum(xi*psi_out(:,:,ib2))*hx_e**2
      xi = xi -ss*psi_out(:,:,ib2)
    end do

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
      psi_t = psi
      do ix1 = -nx_e,nx_e
        do ix2 = -nx_e,nx_e
          psi(ix1,ix2) = 0.5d0*(psi_t(ix1,ix2)+psi_t(ix2,ix1))
        end do
      end do
      
      ss = sum(psi**2)*hx_e**2
      psi = psi/sqrt(ss)
      do ib2 = 1, ib-1
        ss = sum(psi*psi_out(:,:,ib2))*hx_e**2
        psi = psi -ss*psi_out(:,:,ib2)
      end do
      ss = sum(psi**2)*hx_e**2
      psi = psi/sqrt(ss)
      if(icg == ncg)exit
      
      tpsi_e(-nx_e:nx_e,-nx_e:nx_e) = psi(-nx_e:nx_e,-nx_e:nx_e)
      call hpsi_e
      lambda = sum(psi*htpsi_e)*hx_e**2
      xi = lambda*psi - htpsi_e
      do ib2 = 1, ib-1
        ss = sum(xi*psi_out(:,:,ib2))*hx_e**2
        xi = xi -ss*psi_out(:,:,ib2)
      end do
      
      xixi = sum(xi**2)*hx_e**2

      if(xixi < res_c)then
        write(*,*)'break cg_bo'
        exit
      end if


    end do

    psi_out(:,:,ib) = psi
    eps_out(ib) = lambda


  end do


end subroutine cg_bo
!-----------------------------------------------------------------------------
subroutine reordering(psi_t, eps_t)
  use global_variables
  implicit none
  real(8),intent(inout) :: psi_t(-nx_e:nx_e,-nx_e:nx_e,n_bo)
  real(8),intent(inout) :: eps_t(n_bo)
  real(8) :: psi_tmp(-nx_e:nx_e,-nx_e:nx_e,n_bo)
  real(8) :: eps_tmp(n_bo)
  integer :: itable_ib(n_bo), itable_ib_tmp

  integer :: ib, ib2

  
  do ib = 1, n_bo
    itable_ib(ib) = ib
  end do

  do ib = 1, n_bo
    do ib2 = 2, n_bo
      if(eps_t(ib2) < eps_t(ib2-1))then
        itable_ib_tmp = itable_ib(ib2)
        itable_ib(ib2) = itable_ib(ib2-1)
        itable_ib(ib2-1) = itable_ib_tmp
      end if
    end do
  end do


  psi_tmp = psi_t
  eps_tmp = eps_t
  do ib = 1, n_bo
    psi_t(:,:,ib) = psi_tmp(:,:,itable_ib(ib))
    eps_t(ib) = eps_tmp(itable_ib(ib))
  end do

end subroutine reordering
!-----------------------------------------------------------------------------
subroutine calc_Etot_td(Etot_t,it)
  use global_variables
  implicit none
  real(8),intent(out) :: Etot_t
  integer,intent(in) :: it
  real(8) :: Et_t

  Et_t = Et(it)
  call calc_Hartree_pot_td

  ztpsi_e(-nx_e:nx_e,-nx_e:nx_e) = zwfn_e(-nx_e:nx_e,-nx_e:nx_e)
  call zhpsi_e(Et_t)

  ztpsi_p(-nx_p:nx_p) = zwfn_p(-nx_p:nx_p)
  call zhpsi_p(Et_t)

  Etot_t = sum(conjg(zwfn_e)*zhtpsi_e)*hx_e**2 &
         + sum(conjg(zwfn_p)*zhtpsi_p)*hx_p    &
         - sum(vh_p*rho_p)*hx_p

  

end subroutine calc_Etot_td
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine calc_bo_quantity(pop_bo,rho_ovlp_bo)
  use global_variables
  implicit none
  real(8),intent(out) :: pop_bo(n_bo), rho_ovlp_bo(n_bo,n_bo)
  complex(8) :: zchi_bo(-nx_p:nx_p,n_bo)
  integer :: ibo, jbo, ixp


  do ixp = -nx_p,nx_p

    do ibo = 1, n_bo
      zchi_bo(ixp,ibo) = sum(wfn_bo(:,:,ixp,ibo)*zwfn_e)*zwfn_p(ixp)*hx_e**2

    end do

  end do


! population
  do ibo = 1, n_bo
    pop_bo(ibo) = sum(abs(zchi_bo(:,ibo))**2)*hx_p
  end do

! coherence
  do ibo = 1, n_bo
    do jbo = ibo+1,n_bo
      rho_ovlp_bo(ibo,jbo) = sum(abs(zchi_bo(:,jbo)*zchi_bo(:,ibo))**2)*hx_p
      rho_ovlp_bo(jbo,ibo) = rho_ovlp_bo(ibo,jbo)
    end do
  end do

  

end subroutine calc_bo_quantity
