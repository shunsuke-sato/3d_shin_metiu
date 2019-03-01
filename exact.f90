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

  call reduced_density_matrix

end program main
!===============================================================================
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
!===============================================================================
subroutine initialization
  use global_variables
  implicit none
  integer :: ixe1, ixe2, ixp
  real(8) :: xx_e1, xx_e2, xx_p

  write(*,*)"Start: initialization"

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
          -erf_x(abs(xx_e2+0.5d0*Ldist_m)/Rl_m)/Rl_m &
          +erf_x(abs(xx_e2-xx_e1)/Re_m)/Re_m 


      end do
    end do
  end do

  write(*,*)"End: initialization"


end subroutine initialization
!===============================================================================
subroutine ground_state
  use global_variables
  implicit none
  integer :: ixp,ixe1,ixe2
  real(8) :: ss

  write(*,*)"Start: ground_state"

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

  open(20,file="density_1e.out")
  do ixe1 = -nx_e, nx_e
    write(20,"(999e26.16e3)")ixe1*hx_e,sum(wfn(ixe1,:,:)**2)*hx_e*hx_p
  end do
  close(20)

  open(20,file="density_1p.out")
  do ixp = -nx_p, nx_p
    write(20,"(999e26.16e3)")ixp*hx_p,sum(wfn(:,:,ixp)**2)*hx_e**2
  end do
  close(20)

  open(20,file="density_2pe.out")
  do ixp = -nx_p, nx_p
    do ixe1 = -nx_e, nx_e
      write(20,"(999e26.16e3)")ixp*hx_p,ixe1*hx_e,sum(wfn(:,ixe1,ixp)**2)*hx_e
    end do
    write(20,*)
  end do
  close(20)

  open(20,file="density_2ee.out")
  do ixe1 = -nx_e, nx_e
    do ixe2 = -nx_e, nx_e
      write(20,"(999e26.16e3)")ixe1*hx_e,ixe2*hx_e,sum(wfn(ixe1,ixe2,:)**2)*hx_p
    end do
    write(20,*)
  end do
  close(20)


  write(*,*)"End: ground_state"

end subroutine ground_state
!===============================================================================
subroutine conjugate_gradient
  use global_variables
  implicit none
  integer,parameter :: ncg = 300
  integer :: icg
  real(8) :: psi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
  real(8) :: xi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
  real(8) :: phi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
  real(8) :: phi_old(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
  real(8) :: lambda, xixi, xixi_old, gamma
  real(8) :: ss, aa, bb, theta 

  psi = wfn

  tpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p) = psi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
  call hpsi
  lambda = sum(psi*htpsi)*hx_p*hx_e**2
  xi = lambda*psi-htpsi

  xixi = sum(xi**2)*hx_p*hx_e**2
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
    ss = sum(phi*psi)*hx_p*hx_e**2
    phi = phi -ss*psi
    ss = sum(phi**2)*hx_p*hx_e**2
    phi = phi/sqrt(ss)

    bb = 2d0*sum(phi*htpsi)*hx_p*hx_e**2
    tpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p) = phi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
    call hpsi
    aa = sum(phi*htpsi)*hx_p*hx_e**2 -lambda
    aa = - aa
    if(aa /= 0d0)then
      theta = 0.5d0*atan(bb/aa)
    else
      theta = 0.25d0*pi
    end if
    psi = cos(theta)*psi + sin(theta)*phi

    ss = sum(psi**2)*hx_p*hx_e**2
    psi = psi/sqrt(ss)
    if(icg == ncg)exit

    tpsi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p) = psi(-nx_e:nx_e,-nx_e:nx_e,-nx_p:nx_p)
    call hpsi
    lambda = sum(psi*htpsi)*hx_p*hx_e**2
    xi = lambda*psi - htpsi 

    xixi = sum(xi**2)*hx_p*hx_e**2

    write(*,*)icg,lambda,xixi

  end do

   wfn = psi


end subroutine conjugate_gradient
!===============================================================================
subroutine reduced_density_matrix
  use global_variables
  implicit none
  integer :: n,ix,iy,i,j
  real(8),allocatable :: rho_dm(:,:), w(:)
  real(8),allocatable :: work(:)
  integer :: lwork, info


! nuclear part
  n = 2*nx_p + 1
  lwork = n**2
  allocate(rho_dm(n,n), w(n), work(lwork))

  do ix = -nx_p,nx_p
    i = ix + nx_p + 1
    do iy = -nx_p,nx_p
      j = iy + nx_p + 1 
      rho_dm(i,j) = sum(wfn(:,:,ix)*wfn(:,:,iy))*hx_e**2*hx_p
    end do
  end do

  open(20,file="density_1p_chk.out")
  do ix = -nx_p, nx_p
    i = ix + nx_p + 1
    write(20,"(999e26.16e3)")ix*hx_p,rho_dm(i,i)/hx_p
  end do
  close(20)

!  call dsyev('N', 'U', n, rho_dm, n, w, work, lwork, info)
  call dsyev('V', 'U', n, rho_dm, n, w, work, lwork, info)
  open(40,file="natural_occ_1rdm_p.out")
  write(40,"(A,2x,e26.16e3)")"# Total population =",sum(w)
  do i = 1, n
    write(40,"(I8,2x,999e26.16e3)")i,w(n-i+1)
  end do
  close(40)

  open(40,file="natural_orb_1rdm_p.out")
  do ix = -nx_p,nx_p
    i = ix + nx_p + 1
    write(40,"(999e26.16e3)")ix*hx_p,rho_dm(i,n),rho_dm(i,n-1),rho_dm(i,n-2)
  end do
  close(40)

  return
! electronic part
  n = 2*nx_e + 1
  lwork = n**2
  deallocate(rho_dm, w, work)
  allocate(rho_dm(n,n), w(n), work(lwork))

  do ix = -nx_e,nx_e
    i = ix + nx_e + 1
    do iy = -nx_e,nx_e
      j = iy + nx_e + 1 
      rho_dm(i,j) = sum(wfn(ix,:,:)*wfn(iy,:,:))*hx_e**2*hx_p
    end do
  end do

  open(20,file="density_1e_chk.out")
  do ix = -nx_e, nx_e
    i = ix + nx_e + 1
    write(20,"(999e26.16e3)")ix*hx_e,rho_dm(i,i)/hx_e
  end do
  close(20)

!  call dsyev('N', 'U', n, rho_dm, n, w, work, lwork, info)
  call dsyev('V', 'U', n, rho_dm, n, w, work, lwork, info)
  open(40,file="natural_occ_1rdm_e.out")
  write(40,"(A,2x,e26.16e3)")"# Total population =",sum(w)
  do i = 1, n
    write(40,"(I8,2x,999e26.16e3)")i,w(n-i+1)
  end do
  close(40)

  open(40,file="natural_orb_1rdm_e.out")
  do ix = -nx_e,nx_e
    i = ix + nx_e + 1
    write(40,"(999e26.16e3)")ix*hx_e,rho_dm(i,n),rho_dm(i,n-1),rho_dm(i,n-2)
  end do
  close(40)
  

end subroutine reduced_density_matrix
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

!$omp parallel
!$omp do private(ixp,ixe1,ixe2)
  do ixp = -nx_p,nx_p
    do ixe2 = -nx_e, nx_e
      do ixe1 = -nx_e, nx_e


        htpsi(ixe1,ixe2,ixp) = l0pe*tpsi(ixe1,ixe2,ixp) &
          +l1p*(tpsi(ixe1,ixe2,ixp+1)+tpsi(ixe1,ixe2,ixp-1)) &
          +l2p*(tpsi(ixe1,ixe2,ixp+2)+tpsi(ixe1,ixe2,ixp-2)) &
          +l1e*(tpsi(ixe1+1,ixe2,ixp)+tpsi(ixe1-1,ixe2,ixp) &
               +tpsi(ixe1,ixe2+1,ixp)+tpsi(ixe1,ixe2-1,ixp) ) &
          +l2e*(tpsi(ixe1+2,ixe2,ixp)+tpsi(ixe1-2,ixe2,ixp) &
               +tpsi(ixe1,ixe2+2,ixp)+tpsi(ixe1,ixe2-2,ixp) )

      end do
    end do
  end do
!$omp end do

!$omp do private(ixp,ixe1,ixe2)
  do ixp = -nx_p,nx_p
    do ixe2 = -nx_e, nx_e
      do ixe1 = -nx_e, nx_e
        htpsi(ixe1,ixe2,ixp) = htpsi(ixe1,ixe2,ixp) &
          + vpot(ixe1,ixe2,ixp)*tpsi(ixe1,ixe2,ixp)
      end do
    end do
  end do
!$omp end do

!$omp end parallel

end subroutine hpsi
!===============================================================================
!===============================================================================
!===============================================================================
