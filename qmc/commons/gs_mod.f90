module gs_mod !NC

  implicit none
  save
  
  ! gs_delta    small number to prevent singularities at pot. centers
  ! gs_npot     num. of gaussians
  ! gs_ndim     number of dimensions (must be equal to ndim)
  ! gs_scale    scale factor of the sums
  ! gs_ncent    number of centers of each gaussian type                         (1:gs_npot)
  ! gs_vb                       bias value of the sum
  ! gs_v0, gs_rho, gs_s         shape parameters of each gaussian type          (1:gs_npot)
  ! gs_cent                     gaussian centers        (1:gs_ndim, 1:maxval(gs_ncent), 1:gs_npot)      

  double precision, parameter    :: gs_delta = 1.0d-15

  integer                        :: gs_npot, gs_ndim 
  integer, allocatable           :: gs_ncent(:)

  double precision               :: gs_scale, gs_vb
  double precision, allocatable  :: gs_v0(:), gs_rho(:), gs_s(:), gs_cent(:,:,:) 

end module gs_mod
