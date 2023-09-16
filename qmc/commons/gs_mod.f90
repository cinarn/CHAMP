module gs_mod !NC

  implicit none
  save
  
  ! gs_delta    small number to prevent singularities at pot. centers
  ! gs_npot     num. of gaussians
  ! gs_ndim     number of dimensions (must be equal to ndim)
  ! gs_scale    scale factor of the sums
  ! gs_vb       bias values of the sum
  ! gs_ncent    number of centers of each gaussian type                         (1:gs_npot)
  ! gs_v0, gs_rho, gs_s         shape parameters of each gaussian type          (1:gs_npot)
  ! gs_cent                     gaussian centers        (1:gs_ndim, 1:maxval(gs_ncent), 1:gs_npot)

  ! hsf_nf      number of heaviside step functions to enclose the system
  ! hsf_type    type of the enclosing (1: linear, 2:radial)         (1:hsf_nf)
  ! hsf_c       curvature of the step functions                     (1:hsf_nf)
  ! hsf_gargs   geometry arguments of the step functions,           (1:2, 1:hsf_nf)
  !             (radius for radial, the angles of crossings for linear)
  ! hsf_cent    center coordinates for the step functions           (1:gs_ndim, 1:hsf_nf) 

  double precision, parameter    :: gs_delta = 1.0d-15

  integer                        :: gs_npot, gs_ndim, hsf_nf
  integer, allocatable           :: gs_ncent(:), hsf_type(:)

  double precision               :: gs_scale, gs_vb
  double precision, allocatable  :: gs_v0(:), gs_rho(:), gs_s(:), gs_cent(:,:,:), hsf_cent(:,:), hsf_gargs(:,:), hsf_c(:)

  interface cosdeg
    module procedure cosdeg_fun
  end interface

  interface sindeg
    module procedure sindeg_fun
  end interface

  double precision, parameter, private :: pi = 4.d0*datan(1.d0)

contains
  function cosdeg_fun(t)
    double precision, intent(in) :: t
    double precision :: cosdeg_fun
    cosdeg_fun = cos(t*pi/180.d0)
  end function

  function sindeg_fun(t)
    double precision, intent(in) :: t
    double precision :: sindeg_fun
    sindeg_fun = sin(t*pi/180.d0)
  end function   
end module gs_mod

