!! This module

module IC_mod
  use, intrinsic :: iso_fortran_env
  use k_Rosseland_mod, only : k_Ross_Freedman, k_Ross_Valencia, gam_Parmentier, Bond_Parmentier
  implicit none
  !! Precision variables
  integer, parameter :: dp = REAL64


  public :: IC_profile
  private :: adiabat_correction, Parmentier_IC, Guillot_IC, Tdeep_IC, Iso_IC

contains

  subroutine IC_profile(iIC,corr,nlay,p0,pl,pe,k_V,k_IR,Tint,mu,Tirr,grav,fl,Tl,prc)
    implicit none

    !! Input flags
    integer, intent(in) :: iIC
    logical, intent(in) :: corr

    !! Input quantities
    integer, intent(in) :: nlay
    real(kind=dp), intent(in) :: p0, Tint, mu, Tirr, grav, fl
    real(kind=dp), dimension(nlay), intent(in) :: pl
    real(kind=dp), dimension(nlay+1), intent(in) :: pe
    real(kind=dp), dimension(nlay), intent(in) :: k_V, k_IR

    !! Output quantities
    real(kind=dp), dimension(nlay), intent(out) :: Tl
    real(kind=dp), intent(out) :: prc

    select case(iIC)

    case(1)
      call Iso_IC(nlay,Tint,Tl)
    case(2)
      call Tdeep_IC(nlay,k_V(1),k_IR(1),Tirr,Tl)
    case(3)
      call Guillot_IC(nlay,p0,pl,k_V(1),k_IR(1),Tint,mu,Tirr,grav,fl,Tl)
    case(4)
      call Parmentier_IC(nlay,pl,pe,Tint,mu,Tirr,grav,Tl)
    case default
      print*, 'Invalid IC integer in IC_mod, stopping'
      stop

    end select

    if (corr .eqv. .True.) then
      call adiabat_correction(nlay,Tl,pl,prc)
    else
      prc = p0
    end if

  end subroutine IC_profile

  subroutine Iso_IC(nlay,Tint,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(kind=dp), intent(in) :: Tint

    !! Output quantities
    real(kind=dp), dimension(nlay), intent(out) :: Tl

    Tl(:) = Tint

  end subroutine Iso_IC

  subroutine Tdeep_IC(nlay,k_V,k_IR,Tirr,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(kind=dp), intent(in) :: k_V, k_IR, Tirr

    !! Output quantities
    real(kind=dp), dimension(nlay), intent(out) :: Tl

    !! Work variables
    real(kind=dp) :: gam

    gam = k_V/k_IR

    Tl(:) = Tirr * (3.0_dp/4.0_dp * (1.0_dp/gam + 2.0_dp/3.0_dp))**(1.0_dp/4.0_dp)

  end subroutine Tdeep_IC

  subroutine Guillot_IC(nlay,p0,pl,k_V,k_IR,Tint,mu,Tirr,grav,fl,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(kind=dp), intent(in) :: p0, k_V, k_IR, Tint, mu, Tirr, grav, fl
    real(kind=dp), dimension(nlay), intent(in) :: pl

    !! Output quantities
    real(kind=dp), dimension(nlay), intent(out) :: Tl

    !! Work variables
    real(kind=dp), dimension(nlay) :: tau_IRl
    real(kind=dp) :: gam, tau0

    gam = k_V/k_IR
    tau0 = k_IR/grav * p0

    tau_IRl(:) = fl * (pl(:)/p0 * tau0) + (1.0_dp - fl) * (pl(:)/p0 * tau0)**2

    Tl(:) = ((3.0_dp/4.0_dp) * Tint**4 * (tau_IRl(:) + 2.0_dp/3.0_dp))
    Tl(:) = Tl(:) + (mu * 3.0_dp * Tirr**4)/4.0_dp *  &
         & (2.0_dp/3.0_dp + mu/gam + ((gam/(3.0_dp*mu)) - mu/gam) * exp(-gam*tau_IRl(:)/mu))
    Tl(:) = Tl(:)**(1.0_dp/4.0_dp)

  end subroutine Guillot_IC

  !! This subroutine follows Parmentier & Guillot (2014, 2015) non-grey picket fence scheme
  subroutine Parmentier_IC(nlay,pl,pe,Tint,mu,Tirr,grav,Tl)
    implicit none

    integer, intent(in) :: nlay
    real(kind=dp), dimension(nlay), intent(in) :: pl
    real(kind=dp), dimension(nlay+1), intent(in) :: pe
    real(kind=dp), intent(in) :: Tint, mu, Tirr, grav


    real(kind=dp), dimension(nlay), intent(out) :: Tl

    integer, parameter :: table_num = 2
    real(kind=dp), parameter :: met = 0.0_dp

    real(kind=dp) :: Teff0, Teff, Tmu, Bond, Tskin
    real(kind=dp), dimension(3) :: gam_V, Beta_V
    real(kind=dp), dimension(2) :: Beta
    real(kind=dp) :: gam_1, gam_2, gam_P, tau_lim

    integer :: i, j
    real(kind=dp) :: a0, a1, b0, A, B, At1, At2, taul
    real(kind=dp), dimension(3) :: a2, a3, b1, b2, b3, Av1, Av2
    real(kind=dp), dimension(3) :: C, D, E
    real(kind=dp), dimension(nlay+1) :: tau
    real(kind=dp), dimension(nlay) :: kRoss

    !! Effective temperature parameter
    Tmu = (mu * Tirr**4)**(1.0_dp/4.0_dp)
    Teff0 = (Tint**4 + mu * Tirr**4)**(1.0_dp/4.0_dp)

    !! Find Bond albedo of planet
    call Bond_Parmentier(Teff0, grav, Bond)

    Teff = (Tint**4 + (1.0_dp - Bond) * mu * Tirr**4)**(1.0_dp/4.0_dp)

    Teff = 2100.0_dp

    !! Find the V band gamma, beta and IR gamma and beta ratios for this profile
    ! Passed mu, so make lat = acos(mu) and lon = 0
    call gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta, gam_1, gam_2, gam_P, tau_lim)


    gam_V(:) = gam_V(:) / mu

    !! Hard work starts here - first calculate all the required coefficents
    At1 = gam_1**2*log(1.0_dp + 1.0_dp/(tau_lim*gam_1))
    At2 = gam_2**2*log(1.0_dp + 1.0_dp/(tau_lim*gam_2))
    Av1(:) = gam_1**2*log(1.0_dp + gam_V(:)/gam_1)
    Av2(:) = gam_2**2*log(1.0_dp + gam_V(:)/gam_2)

    a0 = 1.0_dp/gam_1 + 1.0_dp/gam_2

    a1 = -1.0_dp/(3.0_dp * tau_lim**2) * (gam_P/(1.0_dp-gam_P) * (gam_1 + gam_2 - 2.0_dp)/(gam_1 + gam_2) &
    &  + (gam_1 + gam_2)*tau_lim - (At1 + At2)*tau_lim**2)

    a2(:) = tau_lim**2/(gam_P*gam_V(:)**2) &
    &  * ((3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*(gam_1+gam_2) &
    & - 3.0_dp*gam_V(:)*(6.0_dp*gam_1**2*gam_2**2-gam_V(:)**2*(gam_1**2+gam_2**2))) &
    & / (1.0_dp-gam_V(:)**2 * tau_lim**2)

    a3(:) = -tau_lim**2*(3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*(Av2(:)+Av1(:)) &
     &/(gam_P*gam_V(:)**3*(1.0_dp-gam_V(:)**2*tau_lim**2))

    b0 = 1.0_dp/(gam_1*gam_2/(gam_1-gam_2)*(At1-At2)/3.0_dp-(gam_1*gam_2)**2/sqrt(3.0_dp*gam_P)-(gam_1*gam_2)**3 &
    & / ((1.0_dp-gam_1)*(1.0_dp-gam_2)*(gam_1+gam_2)))

    b1(:) = gam_1*gam_2*(3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*tau_lim**2 &
    & / (gam_P*gam_V(:)**2*(gam_V(:)**2*tau_lim**2-1.0_dp))

    b2(:) = 3.0_dp*(gam_1+gam_2)*gam_V(:)**3/((3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2))

    b3(:) = (Av2(:)-Av1(:))/(gam_V(:)*(gam_1-gam_2))

    A = 1.0_dp/3.0_dp*(a0+a1*b0)
    B = -1.0_dp/3.0_dp*(gam_1*gam_2)**2/gam_P*b0
    C(:) = -1.0/3.0_dp*(b0*b1(:)*(1.0_dp+b2(:)+b3(:))*a1+a2(:)+a3(:))
    D(:) = 1.0/3.0_dp*(gam_1*gam_2)**2/gam_P*b0*b1(:)*(1.0_dp+b2(:)+b3(:))
    E(:) = (3.0_dp-(gam_V(:)/gam_1)**2)*(3.0_dp-(gam_V(:)/gam_2)**2)/(9.0_dp*gam_V(:)*((gam_V(:)*tau_lim)**2-1.0_dp))

    ! T-p structure calculation - we follow exactly V. Parmentier's method
    ! Estimate the skin temperature by setting tau = 0
    tau(1) = 0.0_dp
    Tskin = 3.0_dp*Tint**4/4.0_dp*(tau(1)+A+B*exp(-tau(1)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
    & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(1)/tau_lim)+E(:)*exp(-gam_V(:)*tau(1))))
    Tskin = Tskin**(1.0_dp/4.0_dp)
    ! Estimate the opacity TOA at the skin temperature - assume this is = first layer optacity
    call k_Ross_Freedman(Tskin, pl(1), met, kRoss(1))
    !call k_Ross_Valencia(Tskin, pe(1), met, kRoss(1))

    ! Recalculate the upmost tau with new kappa
    tau(1) = kRoss(1)/grav * pl(1)
    ! More accurate layer T at uppermost layer
    Tl(1) = 3.0_dp*Tint**4/4.0_dp*(tau(1)+A+B*exp(-tau(1)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
    & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(1)/tau_lim)+E(:)*exp(-gam_V(:)*tau(1))))
    Tl(1) = Tl(1)**(1.0_dp/4.0_dp)

    ! Now we can loop in optical depth space to find the T-p profile
    do i = 2, nlay
      ! Initial guess for layer
      call k_Ross_Freedman(Tl(i-1), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
      !call k_Ross_Valencia(Tl(i-1), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
      tau(i) = tau(i-1) + kRoss(i)/grav * (pl(i) - pl(i-1))
      Tl(i) = 3.0_dp*Tint**4/4.0_dp*(tau(i)+A+B*exp(-tau(i)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
      & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(i)/tau_lim)+E(:)*exp(-gam_V(:)*tau(i))))
      Tl(i) = Tl(i)**(1.0_dp/4.0_dp)
      ! Convergence loop
      do j = 1, 5
        call k_Ross_Freedman(sqrt(Tl(i-1)*Tl(i)), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
        !call k_Ross_Valencia(sqrt(Tl(i-1)*T(i)), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
        tau(i) = tau(i-1) + kRoss(i)/grav * (pl(i) - pl(i-1))
        Tl(i) = 3.0_dp*Tint**4/4.0_dp*(tau(i)+A+B*exp(-tau(i)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
        & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(i)/tau_lim)+E(:)*exp(-gam_V(:)*tau(i))))
        Tl(i) = Tl(i)**(1.0_dp/4.0_dp)

      end do

    end do

  end subroutine Parmentier_IC

  !! Subroutine that corrects for adiabatic region following Parmentier & Guillot (2015)
  subroutine adiabat_correction(nlay,Tl,pl,prc)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(kind=dp), dimension(nlay), intent(in) ::  pl

    !! Output quantities
    real(kind=dp), dimension(nlay), intent(inout) :: Tl
    real(kind=dp), intent(out) :: prc

    !! Work variables
    integer :: k, iRC, iRC1, taurc
    real(kind=dp), dimension(nlay) :: gradrad, gradad

    do k = 1, nlay-1
      gradrad(k) = (log10(Tl(k))-log10(Tl(k+1)))/(log10(pl(k))-log10(pl(k+1)))
      gradad(k) = 0.32_dp - 0.10_dp*Tl(k)/3000.0_dp
      !print*, k, gradrad(k), gradad(k)
    end do
    gradrad(nlay) = 0.0_dp
    gradad(nlay) = 0.0_dp

    iRC = nlay-1
    iRC1 = nlay-1

    do k = nlay-1, 1, -1
      if (IRC1 <= k+1) then
        if (gradrad(k) > 0.7_dp*gradad(k)) then
          iRC1 = k
        endif
        if (gradrad(k) > 0.98_dp*gradad(k)) then
         iRC = k
         prc = pl(iRC)
        endif
      end if
    end do

    if (iRC < nlay) then
      do k = iRC, nlay-1
        gradad(k) = 0.32_dp - 0.10_dp*Tl(k)/3000.0_dp
        if (gradad(k) < 0.0_dp) then
          gradad(k) = 0.0_dp
        end if
        !if (pl(k) > 1.0_dp*1e6_dp) then
        Tl(k+1)=Tl(k)*(pl(k+1)/pl(k))**gradad(k)
      !end if
      end do
    end if

  end subroutine adiabat_correction


end module IC_mod

program IC_mod_test
  use IC_mod, only : IC_profile
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: iIC, nlay, nlay1, k
  logical :: corr

  real(kind=dp) :: p0, mu, Tint, Tirr, fl, grav, prc, Teff
  real(kind=dp), allocatable, dimension(:) :: pe, pl, Tl, k_V_l, k_IR_l

  integer :: u, ua, ub
  character(len=20) :: a_sh, b_sh
  real(kind=dp), allocatable, dimension(:) :: a, b

  ! Number of layers and edges
  nlay = 50
  nlay1 = nlay + 1

  !! Read in sigma hybrid grid values
  a_sh = 'sig_hyb_HJ_51_a.txt'
  b_sh = 'sig_hyb_HJ_51_b.txt'

  open(newunit=ua,file=trim(a_sh),action='read')
  open(newunit=ub,file=trim(b_sh),action='read')

  allocate(a(nlay1),b(nlay1))
  do k = 1, nlay1
    read(ua,*) a(k)
    read(ub,*) b(k)
  end do

  ! Contruct pressure array in pa
  ! Surface pressure (pa)
  p0 = 1e8_dp
  allocate(pe(nlay1),pl(nlay))
  do k = 1, nlay1
    pe(k) = a(k) + b(k)*p0
    !print*, pe(k)/1e5_dp
  end do
  ! Pressure layers
  do k = 1, nlay
    pl(k) = (pe(k) + pe(k+1)) / 2.0_dp
  end do

  allocate(Tl(nlay))
  allocate(k_V_l(nlay),k_IR_l(nlay))

  ! sw Zenith angle
  mu = 1.0_dp / sqrt(3.0_dp)

  Tirr = 200.0_dp ! Irradiation temperature
  Tint = 100.0_dp ! Internal temperature

  Teff = (Tint**4 + mu*Tirr**4)**(1.0_dp/4.0_dp)

  k_IR_l(:) = 1e-3_dp
  k_V_l(:) = 6e-4_dp * sqrt(Tirr/2000.0_dp)

  grav = 25.0_dp
  fl = 1.0_dp!1.0_dp/2000.0_dp
  Tl = 0.0_dp

  iIC = 4
  corr = .True.

  call IC_profile(iIC,corr,nlay,p0,pl,pe,k_V_l,k_IR_l,Tint,mu,Tirr,grav,fl,Tl,prc)

  open(newunit=u,file='FMS_IC.out',action='readwrite')
  do k = 1, nlay
    write(u,*) k, pl(k), Tl(k), prc
  end do
  close(u)

  print*, Tirr, Tint, Teff, grav, mu,  (mu*Tirr**4)**(1.0_dp/4.0_dp), prc/1e5


end program IC_mod_test
